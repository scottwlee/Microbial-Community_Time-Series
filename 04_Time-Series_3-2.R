################################################################################
# [Phase 5 - Part 5 (b) (Master Version): Macro-ecological Phenology (Standardized)]
# 목적: 10년 치 전체 데이터를 월별(Month) 축으로 투영하여, 
#       GMM으로 분류된 세 그룹의 계절적 출현 패턴(Niche Partitioning)을 분석함.
# 특징: 
#   1) [True Zero-Included] Part 4, 7과 완벽히 동일한 결측치 0 복원 로직 적용.
#   2) [개별 ASV 표준화] Z-score 변환을 도입하여 극우점종 편향을 제거하고,
#      순수한 '계절적 타이밍(Phenological Synchrony)'만을 추출하여 시각화.
#   3) [enable_save_outputs] 원터치 마스터 스위치로 파일 저장 완벽 제어.
################################################################################

options(stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS 
#################################################
summary_method  <- "median"  
g_clusters      <- 3 

# 논문의 통일성을 위해 Zero-exclusion 스위치 제거 (강제 Zero-Included)
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/05_Seasonal_Phenology")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

log_file  <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_Standardized_Phenology_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_Standardized_Phenology_GAM.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # -------------------------------------------------------------------
  # Section 1. Data Load & Clustering (Synchronized)
  # -------------------------------------------------------------------
  log_msg("Step 1: Loading Data & Applying GMM Clustering...")
  
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # IS 대푯값은 상호작용이 존재하는 날(결측 아님)을 기준으로 추출
  df_merged_is <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  asv_summary <- df_merged_is %>%
    group_by(ASV) %>%
    summarise(
      agg_IS   = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE)
    ) %>% 
    filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  mapped_clusters <- factor(gmm_model$classification, 
                            levels = order(gmm_model$parameters$mean), 
                            labels = c("1_Negative", "2_Neutral", "3_Positive"))
  asv_clusters <- asv_summary %>%
    mutate(Cluster = mapped_clusters) %>%
    dplyr::select(ASV, Cluster)
  
  # -------------------------------------------------------------------
  # Section 2. True Zero-Included Monthly Aggregation
  # -------------------------------------------------------------------
  log_msg("Step 2: Aggregating Monthly Abundance with True Zero-Inclusion...")
  
  # 1) 월별 총 샘플링 횟수 도출
  month_counts <- df_ab_raw %>% 
    mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Date, Month) %>% 
    dplyr::distinct() %>% 
    group_by(Month) %>% 
    summarise(Num_Samples = n(), .groups = "drop")
  
  # 2) ASV별 Zero-included 평균 산출
  df_pheno_base <- df_ab_raw %>%
    mutate(Month = as.numeric(format(as.Date(Sample_Date), "%m"))) %>%
    group_by(ASV_ID, Month) %>%
    summarise(Sum_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(ASV_ID, Month = 1:12, fill = list(Sum_Abund = 0)) %>%
    # [수정된 부분] "AS" -> "ASV" 로 변경
    inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>% 
    left_join(month_counts, by = "Month") %>%
    mutate(Mean_Abund = Sum_Abund / Num_Samples)
  
  # -------------------------------------------------------------------
  # Section 3. Individual ASV Standardization (Z-score)
  # -------------------------------------------------------------------
  log_msg("Step 3: Standardizing abundance per ASV to extract pure phenology...")
  
  # 연구자님의 아이디어: 극단적인 절대량 차이를 소거하고 형태(Shape)만 비교
  df_standardized <- df_pheno_base %>%
    group_by(ASV_ID) %>%
    mutate(
      # 해당 ASV가 12개월 내내 0이거나 분산이 0일 경우 NaN 방지
      sd_abund = sd(Mean_Abund, na.rm = TRUE),
      Z_Score = ifelse(sd_abund == 0, 0, (Mean_Abund - mean(Mean_Abund, na.rm = TRUE)) / sd_abund)
    ) %>%
    ungroup() %>%
    filter(!is.na(Z_Score))
  
  # -------------------------------------------------------------------
  # Section 4. Visualization with Cyclic GAM Splines
  # -------------------------------------------------------------------
  log_msg("Step 4: Generating Standardized Seasonal Phenology Plot...")
  
  custom_colors <- c("1_Negative" = "#347433", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  
  caption_base <- sprintf(
    "Note: Clusters defined via 1D GMM (G=%d) on %s Temp_IS.\nData points represent Standardized (Z-score) monthly mean abundance per ASV (True Zero-Included).\nStandardization eliminates dominant taxa bias, revealing pure temporal niche partitioning.\nGAM lines fitted with Cyclic Cubic Splines (bs='cc').",
    g_clusters, tools::toTitleCase(summary_method)
  )
  
  p_phenology <- ggplot(df_standardized, aes(x = Month, y = Z_Score, color = Cluster, fill = Cluster)) +
    # 노이즈를 줄이기 위해 알파값 조정
    geom_jitter(alpha = 0.2, size = 1.2, width = 0.2) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), 
                alpha = 0.2, linewidth = 1.5) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(
      title = "Seasonal Phenology & Temporal Niche Partitioning",
      subtitle = "Standardized Population Dynamics (Z-Score of Monthly Means)",
      x = "Month of the Year",
      y = "Standardized Abundance (Z-Score)", 
      caption = caption_base
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
      plot.caption = element_text(hjust = 0, size = 10, color = "gray30", margin = margin(t = 15)),
      panel.grid.minor.x = element_blank()
    )
  
  # -------------------------------------------------------------------
  # Section 5. Export 
  # -------------------------------------------------------------------
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_phenology, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", file_plot))
  } else {
    log_msg("[SAFE MODE] Plot was generated in RStudio Viewer only (No file saved).")
  }
  
  print(p_phenology)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

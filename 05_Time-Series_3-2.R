# [Github: 04_Time-Series_3-2.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 2 (Master Version): Macro-ecological Phenology (Standardized)]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 5 (b)"에서 이관 및 재정비된 코드입니다.)
#
# 목적: 10년 치 전체 데이터를 월별(Month) 축으로 투영하여, 
#       GMM으로 분류된 세 그룹의 계절적 출현 패턴(Niche Partitioning)을 분석함.
# 특징: 
#   1) [True Zero-Included] 특정 종이 관찰되지 않은 샘플링 날짜를 0으로 강제 복원.
#   2) [개별 ASV 표준화] Z-score 변환을 도입하여 절대적 극우점종 편향을 제거하고,
#      순수한 '계절적 타이밍(Phenological Synchrony)'만을 추출하여 시각화.
#   3) [로깅 고도화] 플롯별 Title, Note 및 군집별 GAM 통계치를 로그 파일에 명시적 기록.
#   4) [다중 플롯 분리] 통합 오버레이 플롯과 개별 군집 플롯을 각각 분리하여 독립 저장.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
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
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택
#     - 옵션: "median", "mean"
#     - 권장: "median" (극단적인 이상치에 덜 민감함)
summary_method  <- "median"  

# [2] GMM 클러스터 개수
#     - 옵션: 양의 정수
#     - 권장: 3 (Negative, Neutral, Positive 생태학적 3분할)
g_clusters      <- 3 

# [3] 결과 저장 마스터 스위치
#     - 옵션: TRUE (지정된 폴더에 플롯과 로그 파일 저장), FALSE (RStudio 뷰어 및 콘솔에만 출력)
#     - 논문의 통일성을 위해 Zero-exclusion 스위치 제거 (강제 Zero-Included)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (자동 적용 - 폴더명 Part2 에 맞게 수정)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_V2_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/02_Seasonal_Phenology")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일 접두사 Part2_ 적용
log_file          <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Log.txt"))
file_plot_overlay <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Overlay.tiff"))
file_plot_neg     <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Negative.tiff"))
file_plot_neu     <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Neutral.tiff"))
file_plot_pos     <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Positive.tiff"))
file_plot_panels  <- file.path(out_dir, paste0("Part2_", toupper(summary_method), "_Standardized_Phenology_Combined_Panels.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Clustering (Synchronized)
  # ------------------------------------------------------------------- #
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
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE),
      .groups = "drop"
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
  
  # ------------------------------------------------------------------- #
  # Section 2. True Zero-Included Monthly Aggregation
  # ------------------------------------------------------------------- #
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
    inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>% 
    left_join(month_counts, by = "Month") %>%
    mutate(Mean_Abund = Sum_Abund / Num_Samples)
  
  # ------------------------------------------------------------------- #
  # Section 3. Individual ASV Standardization (Z-score)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Standardizing abundance per ASV to extract pure phenology...")
  
  # 극단적인 절대량 차이를 소거하고 형태(Shape)만 비교
  df_standardized <- df_pheno_base %>%
    group_by(ASV_ID) %>%
    mutate(
      sd_abund = sd(Mean_Abund, na.rm = TRUE),
      Z_Score = ifelse(sd_abund == 0, 0, (Mean_Abund - mean(Mean_Abund, na.rm = TRUE)) / sd_abund)
    ) %>%
    ungroup() %>%
    filter(!is.na(Z_Score))
  
  # ------------------------------------------------------------------- #
  # Section 4. Extracting GAM Statistics & Logging
  # ------------------------------------------------------------------- #
  log_msg("\n===================================================================")
  log_msg(" [Analysis Parameters & Clustering Info]")
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included"))
  log_msg(sprintf(" - Standardization  : Z-Score transformation per ASV"))
  log_msg(sprintf(" - Total ASVs Mapped: %d", nrow(asv_summary)))
  log_msg("===================================================================\n")
  
  main_title <- "Seasonal Phenology & Temporal Niche Partitioning"
  plot_note  <- sprintf("Note: Clusters defined via 1D GMM (G=%d) on %s Temp_IS. Data points represent Standardized (Z-score) monthly mean abundance per ASV (True Zero-Included). Standardization eliminates dominant taxa bias, revealing pure temporal niche partitioning. GAM lines fitted with Cyclic Cubic Splines (bs='cc', k=12).", g_clusters, tools::toTitleCase(summary_method))
  
  log_msg(sprintf(" [Plot Title] %s", main_title))
  log_msg(sprintf(" [Plot Note]  %s\n", plot_note))
  
  log_msg(" [GAM Statistics per Cluster (Z_Score ~ s(Month, bs='cc'))]")
  for (clust in levels(df_standardized$Cluster)) {
    sub_df <- df_standardized %>% filter(Cluster == clust)
    if (nrow(sub_df) > 12) {
      gam_mod <- mgcv::gam(Z_Score ~ s(Month, bs = "cc", k = 12), data = sub_df, method = "REML")
      gam_sum <- summary(gam_mod)
      pval <- gam_sum$s.table[1, "p-value"]
      pval_str <- ifelse(pval < 0.001, "p < 0.001", format.pval(pval, digits = 3))
      
      log_msg(sprintf("  -> %-12s: edf = %5.2f, F = %5.2f, %s, Deviance Explained = %4.1f%%", 
                      clust, gam_sum$s.table[1, "edf"], gam_sum$s.table[1, "F"], pval_str, gam_sum$dev.expl * 100))
    }
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 5. Visualization (Overlay & Individual Panels)
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Generating Overlay and Individual Panel Plots...")
  
  custom_colors <- c("1_Negative" = "#0065F8", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  
  base_theme <- theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
    panel.grid.minor.x = element_blank()
  )
  
  # 함수: 개별 또는 통합 플롯 생성기 (Note는 로그로 분리했으므로 Caption에서 제거)
  create_pheno_plot <- function(df, plot_title, show_legend = FALSE) {
    p <- ggplot(df, aes(x = Month, y = Z_Score, color = Cluster, fill = Cluster)) +
      geom_jitter(alpha = 0.2, size = 1.2, width = 0.2) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), alpha = 0.2, linewidth = 1.5) +
      scale_color_manual(values = custom_colors) +
      scale_fill_manual(values = custom_colors) +
      scale_x_continuous(breaks = 1:12, labels = month.abb) +
      labs(title = plot_title, subtitle = "Standardized Population Dynamics (Z-Score of Monthly Means)", x = "Month of the Year", y = "Standardized Abundance (Z-Score)") +
      base_theme
    
    if(show_legend) {
      p <- p + theme(legend.position = "top", legend.title = element_blank())
    }
    return(p)
  }
  
  # 1) 전체 통합 오버레이 플롯
  p_overlay <- create_pheno_plot(df_standardized, main_title, show_legend = TRUE)
  
  # 2) 군집별 개별 플롯
  p_neg <- create_pheno_plot(df_standardized %>% filter(Cluster == "1_Negative"), "1_Negative Phenology")
  p_neu <- create_pheno_plot(df_standardized %>% filter(Cluster == "2_Neutral"), "2_Neutral Phenology")
  p_pos <- create_pheno_plot(df_standardized %>% filter(Cluster == "3_Positive"), "3_Positive Phenology")
  
  # 3) 세로 다중 패널 결합 플롯
  p_combined_panels <- plot_grid(p_neg, p_neu, p_pos, ncol = 1, align = "v", labels = c("A", "B", "C"))
  
  # ------------------------------------------------------------------- #
  # Section 6. Export 
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    # 통합 플롯 저장
    ggsave(file_plot_overlay, plot = p_overlay, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    # 개별 군집 플롯 저장
    ggsave(file_plot_neg, plot = p_neg, device = "tiff", dpi = 600, width = 10, height = 6, compression = "lzw")
    ggsave(file_plot_neu, plot = p_neu, device = "tiff", dpi = 600, width = 10, height = 6, compression = "lzw")
    ggsave(file_plot_pos, plot = p_pos, device = "tiff", dpi = 600, width = 10, height = 6, compression = "lzw")
    # 3단 결합 플롯 저장
    ggsave(file_plot_panels, plot = p_combined_panels, device = "tiff", dpi = 600, width = 10, height = 15, compression = "lzw")
    
    log_msg("\n[SUCCESS] Overlay plot, Individual plots, and Combined panel plot saved to disk successfully.")
  } else {
    log_msg("\n[SAFE MODE] Plots were generated in RStudio Viewer only (No files saved).")
  }
  
  print(p_overlay)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

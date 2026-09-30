# [Github: 04_Time-Series_3-5.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 5 (Master Version): Ecological Response Time (Time-Lag) Spectrum]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 8"에서 이관 및 재정비된 코드입니다.)
#
# 목적: 각 ASV가 온도 변화를 겪은 후 실제 증식(Bloom)으로 반응하기까지 걸리는
#       고유 지연 시간(Best_TP, Optimal Time Lag)을 온도 민감도(IS) 군집별로 파악함.
# 특징: 
#   1) [True Zero-Included 동기화] 이전 파트들과 완벽히 동일한 결측치(0) 강제 
#      복원 엔진을 적용하여, GMM 클러스터링 배정 ASV 목록을 100% 일치시킴.
#   2) [통계 검정 자동화] 비모수 전역 검정(Kruskal-Wallis) 및 사후 검정 수행.
#   3) [로깅 고도화] 플롯 하단의 캡션을 제거하고 플롯 Title과 상세 통계/Note를 로그에 통합.
#   4) [색상 일관성] 논문 전반의 테마에 맞춰 1_Negative 그룹에 #0065F8 색상 적용.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "stats")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(mclust)
  library(stats)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] IS 통계량 요약 방식 
#     - 옵션: "median", "mean"
#     - 권장: "median" (극단적인 이상치에 덜 민감함)
summary_method  <- "median"  

# [2] GMM 클러스터(군집) 개수 
#     - 옵션: 양의 정수
#     - 권장: 3 (Negative, Neutral, Positive 생태학적 3분할)
g_clusters      <- 3 

# [3] 결과 저장 마스터 스위치 
#     - 옵션: TRUE (지정된 폴더에 플롯과 로그 파일 자동 저장), FALSE (RStudio 뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 갱신: Part 5)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_V2_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
out_dir    <- file.path(base_dir, "05_Phase5_Output/05_Time_Lag_Distribution")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

# 파일 접두사 Part5_ 적용
log_file  <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_Time_Lag_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_Time_Lag_Distribution.tiff"))
file_csv  <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_Time_Lag_Data.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & ID Mapping (Robust Selection)
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading Data and S-map Summary...")
  
  if(!file.exists(file_abundance)) stop("Error: Abundance Input file not found.")
  if(!file.exists(file_temp_is)) stop("Error: IS Input file not found.")
  if(!file.exists(file_smap_sum)) stop("Error: S-map Summary file not found.")
  
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE)
  
  asv_dict <- df_ab_raw %>% 
    dplyr::select(ASV_Hash, ASV = ASV_ID) %>% 
    dplyr::distinct()
  
  colnames(df_smap_raw)[1] <- "ASV_Hash"
  
  df_smap <- df_smap_raw %>%
    dplyr::inner_join(asv_dict, by = "ASV_Hash") %>%
    dplyr::select(ASV, Best_TP, Rho)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # ------------------------------------------------------------------- #
  # Section 2. Apply GMM Clustering (True Zero-Included Logic)
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Aggregating Data (True Zero-Included) & Applying GMM clustering...")
  
  asv_is_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
                     .groups = "drop")
  
  all_sample_dates <- unique(df_ab$Date)
  asv_ab_summary <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(asv_mean = mean(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  asv_summary <- dplyr::inner_join(asv_is_summary, asv_ab_summary, by = "ASV") %>%
    dplyr::filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_order <- order(gmm_model$parameters$mean) 
  
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)
  
  # ------------------------------------------------------------------- #
  # Section 3. Merge Clusters with Response Time (Best_TP)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Merging Response Time (Lag) data with GMM Clusters...")
  
  df_lag_analysis <- asv_clusters %>%
    dplyr::inner_join(df_smap, by = "ASV") %>%
    dplyr::mutate(Response_Lag_Weeks = abs(Best_TP)) %>%
    dplyr::filter(!is.na(Response_Lag_Weeks))
  
  # ------------------------------------------------------------------- #
  # Section 4. Statistical Testing (Kruskal-Wallis & Post-hoc)
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Performing Statistical Testing (Kruskal-Wallis & Pairwise Wilcoxon)...")
  
  kw_test <- kruskal.test(Response_Lag_Weeks ~ Cluster, data = df_lag_analysis)
  kw_pval <- kw_test$p.value
  kw_stat <- kw_test$statistic
  pval_text <- ifelse(kw_pval < 0.001, "p < 0.001", sprintf("p = %.3f", kw_pval))
  
  lag_summary <- df_lag_analysis %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(Mean_Lag = mean(Response_Lag_Weeks), 
                     Median_Lag = median(Response_Lag_Weeks),
                     Min_Lag = min(Response_Lag_Weeks),
                     Max_Lag = max(Response_Lag_Weeks),
                     .groups = "drop")
  
  plot_title <- "Ecological Response Time (Time-Lag) Spectrum"
  plot_note <- sprintf("'Response Delay' is the absolute value of Best_TP from CCM/S-map analysis. Clusters defined via 1D GMM (G=%d) on True Zero-Included %s Temp_IS.", g_clusters, tools::toTitleCase(summary_method))
  
  log_msg("\n===================================================================")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - GMM Clusters     : %d", g_clusters))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included"))
  log_msg("-------------------------------------------------------------------")
  log_msg(sprintf(" [Plot: %s]", plot_title))
  log_msg(sprintf("  -> Note           : %s", plot_note))
  log_msg(sprintf("  -> Kruskal-Wallis : \u03c7\u00b2 = %.2f, %s", kw_stat, pval_text))
  log_msg("===================================================================\n")
  
  log_msg(" -> Time-Lag Descriptive Statistics:")
  print(lag_summary)
  
  if(kw_pval < 0.05) {
    posthoc_res <- pairwise.wilcox.test(df_lag_analysis$Response_Lag_Weeks, 
                                        df_lag_analysis$Cluster, 
                                        p.adjust.method = "BH", 
                                        exact = FALSE)
    log_msg(" -> Post-hoc Pairwise Wilcoxon Test (FDR adjusted p-values):")
    print(posthoc_res)
  } else {
    log_msg(" -> No significant global difference; skipping post-hoc test.")
  }
  
  # ------------------------------------------------------------------- #
  # Section 5. Visualization (Violin + Jitter + Boxplot)
  # ------------------------------------------------------------------- #
  log_msg("Step 5: Generating Lag Distribution Plot...")
  
  # 논문 전반에 걸친 그룹별 공통 색상 팔레트 강제 적용
  custom_colors <- c("1_Negative" = "#0065F8", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  
  p_lag <- ggplot(df_lag_analysis, aes(x = Cluster, y = Response_Lag_Weeks, fill = Cluster, color = Cluster)) +
    geom_violin(alpha = 0.2, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = seq(0, max(df_lag_analysis$Response_Lag_Weeks), by = 2)) +
    labs(
      title = plot_title,
      subtitle = "Distribution of optimal response delay (Best_TP) across temperature sensitivity groups",
      x = paste0("GMM Cluster (", tools::toTitleCase(summary_method), " IS)"),
      y = "Ecological Response Delay (Weeks)"
      # Note 캡션은 로그로 이관됨 (Rule 4 적용)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
      axis.text.x = element_text(face = "bold", size = 12)
    )
  
  # ------------------------------------------------------------------- #
  # Section 6. Save Results & Export (단일 플롯 처리)
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_lag, device = "tiff", dpi = 600, width = 9, height = 7, compression = "lzw")
    write.csv(df_lag_analysis, file_csv, row.names = FALSE)
    log_msg(paste("[SUCCESS] Plot and CSV Data saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot generated in Viewer only (No files saved).")
  }
  
  print(p_lag)
  log_msg("Analysis complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

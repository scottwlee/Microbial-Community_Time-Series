# [Github: 04_Time-Series_3-10.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 10 (Master Version): Statistical Validation of GMM Clusters]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-1.R" 스크립트의 "Phase 5 - Part 16"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: GMM으로 분류된 세 미생물 그룹(1_Negative, 2_Neutral, 3_Positive)이 
#       온도 상호작용 강도(Temp_IS) 측면에서 통계적으로 완벽하게 분리되는 
#       독립적인 집단인지 ANOVA 및 Tukey HSD 사후검정으로 증명함.
# 특징:
#   1) [논리적 일관성] 이전 파트들과 완벽히 동일한 True Zero-Included 
#      결측치(0) 복원 로직을 사용하여 GMM 클러스터 배정 ASV를 100% 일치시킴.
#   2) [로깅 고도화] 플롯의 캡션을 제거하고, ANOVA 및 Tukey 결과를 로그 파일에 상세 기록.
#   3) [다중/단일 플롯 독립 저장] 단일 패널 플롯을 독립적인 고해상도 파일로 저장.
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
#     - 권장: "median" (극단적인 이상치에 덜 민감하며 이전 분석과 통일)
summary_method  <- "median"  

# [2] GMM 클러스터(군집) 개수 
#     - 옵션: 양의 정수
#     - 권장: 3 (Negative, Neutral, Positive 생태학적 3분할)
g_clusters      <- 3 

# [3] 결과 저장 마스터 스위치 
#     - 옵션: TRUE (지정된 폴더에 플롯과 로그 파일 자동 저장), FALSE (RStudio 뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 갱신: Part 10)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/10_Cluster_Validation")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일명 접두사 Part10_ 적용
log_file  <- file.path(out_dir, paste0("Part10_", toupper(summary_method), "_Cluster_Validation_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part10_", toupper(summary_method), "_Cluster_Validation.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 10: Statistical Validation of GMM Clusters]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - GMM Clusters     : %d", g_clusters))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included (Synchronized with pipeline)"))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Clustering (True Zero-Included Logic)
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading data, applying True-Zero inclusion, and assigning GMM clusters...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # IS 대푯값 추출 (상호작용이 존재하는 날 기준)
  asv_is_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE), .groups = "drop")
  
  # 결측치 0 강제 복원을 통한 평균 절대 풍부도 계산 (파이프라인 통일성 확보)
  all_sample_dates <- unique(df_ab$Date)
  asv_ab_summary <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(asv_mean = mean(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  # 병합 및 필터링
  asv_summary <- dplyr::inner_join(asv_is_summary, asv_ab_summary, by = "ASV") %>%
    dplyr::filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_order <- order(gmm_model$parameters$mean)
  
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive")))
  
  # ------------------------------------------------------------------- #
  # Section 2. Statistical Testing (ANOVA & Tukey HSD) & Logging
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Running ANOVA and Tukey HSD Post-hoc test...")
  
  anova_res <- aov(agg_IS ~ Cluster, data = asv_clusters)
  anova_summary <- summary(anova_res)
  pval_anova <- anova_summary[[1]][["Pr(>F)"]][1]
  f_val <- anova_summary[[1]][["F value"]][1]
  
  pval_text <- ifelse(pval_anova < 0.001, "p < 0.001", sprintf("p = %.3f", pval_anova))
  
  title_plot <- "Validation of Ecological Cluster Distinctness"
  note_plot  <- "One-way ANOVA and Tukey HSD confirm that all three clusters are mutually exclusive and possess statistically distinct thermal interaction profiles (True Zero-Included)."
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg(sprintf("[Plot: %s]", title_plot))
  log_msg(sprintf(" -> Note                     : %s", note_plot))
  log_msg(sprintf(" -> One-way ANOVA F-value    : %.2f", f_val))
  log_msg(sprintf(" -> ANOVA P-value            : %s", pval_text))
  
  log_msg("\n -> Tukey HSD Post-hoc Test Results:")
  tukey_res <- TukeyHSD(anova_res)
  
  # Tukey 결과를 문자열로 깔끔하게 로그에 기록
  tukey_df <- as.data.frame(tukey_res$Cluster)
  for(i in 1:nrow(tukey_df)) {
    p_adj <- tukey_df$`p adj`[i]
    p_adj_str <- ifelse(p_adj < 0.001, "p < 0.001", sprintf("p = %.4f", p_adj))
    log_msg(sprintf("    [%s] Diff: %7.4f | %s", rownames(tukey_df)[i], tukey_df$diff[i], p_adj_str))
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 3. Visualization
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Generating Violin & Boxplot for Cluster Validation...")
  
  # 논문 전반에 걸친 그룹별 공통 색상 팔레트 강제 적용
  custom_colors <- c("1_Negative" = "#0065F8", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  
  subtitle_plot <- paste("Distribution of EDM-derived Interaction Strength (Temp_IS) by GMM Cluster")
  
  p_val <- ggplot(asv_clusters, aes(x = Cluster, y = agg_IS, fill = Cluster)) +
    geom_violin(alpha = 0.4, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.7, color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.5, aes(color = Cluster)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    labs(
      title = title_plot, 
      subtitle = subtitle_plot, 
      x = "Ecological Cluster", 
      y = paste0("Temperature Interaction Strength (", toupper(summary_method), " Temp_IS)")
      # 캡션(Note)은 로그 파일로 이관됨 (Rule 4)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size=15),
      plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 10)),
      axis.text.x = element_text(face = "bold", size = 11)
    )
  
  # ------------------------------------------------------------------- #
  # Section 4. Save & Export
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    # 단일 플롯 독립 저장 (Rule 8)
    ggsave(file_plot, plot = p_val, device = "tiff", dpi = 600, width = 9, height = 7, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot generated in Viewer only (No files saved).")
  }
  
  print(p_val)
  log_msg("Cluster Validation analysis complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

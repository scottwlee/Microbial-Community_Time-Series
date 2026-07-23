# [Github: 04_Time-Series_3-1.R]==================================================#

# ------------------------------------------------------------------- # 
# [Phase 5 - Part 1 (Master Version): Macro-ecological Traits (True Zero-Included)]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 4 (b)"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: 10년 치 전체 데이터를 바탕으로 ASV의 온도 민감도(Temp_IS)에 따른
#       거시적 생태 특성(CV, SD, Mean)의 비선형적 관계를 분석함.
# 특징:
#   1) [핵심 혁신: True Zero-Included] tidyr::complete를 사용하여 특정 종이 
#      관찰되지 않은(결측된) 모든 샘플링 날짜를 0으로 강제 복원한 뒤 진짜 통계를 산출.
#   2) [논문 통일성 확보] 향후 분석과 완벽히 동일한 생태학적 기준(Zero-inclusion) 적용.
#   3) [클러스터링 검증] Silhouette Score, Uncertainty, ANOVA를 통한 3분할 타당성 입증.
#   4) [로깅 고도화] 플롯별 Title, Note 및 통계치(GAM, GMM)를 로그 파일에 명시적으로 기록.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
# [NEW] cluster 패키지 추가 (Silhouette Score 계산용)
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "cluster")
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
  library(cluster)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택 
#     - 옵션: "median", "mean"
#     - 권장: "median" (극단적인 이상치에 덜 민감함)
summary_method  <- "median"  

# [2] 플롯 시각화 색상 모드 선택
#     - 옵션: TRUE (군집별 3색 적용), FALSE (전체 단일 흑백 적용)
enable_gmm_colors <- TRUE 

# [3] GMM 클러스터 개수
#     - 옵션: 양의 정수
#     - 권장: 3 (Negative, Neutral, Positive 생태학적 3분할)
g_clusters      <- 3 

# [4] 이상치 필터링 파라미터 (IQR 기반)
#     - iqr_multiplier : 필터링 강도 (권장: 1.5 - 통상적인 박스플롯 수염 기준)
#     - iqr_filter_mode: 필터링 적용 축 (옵션: "x", "y", "both" | 권장: "y" - 종속변수 극단값만 제거)
iqr_multiplier  <- 1.5
iqr_filter_mode <- "y"

# [5] 결과 저장 마스터 스위치
#     - 옵션: TRUE (지정된 폴더에 플롯과 로그 파일 저장), FALSE (RStudio 뷰어 및 콘솔에만 출력)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 파일명 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/01_Macro_Ecological_Traits")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

color_mode_str <- ifelse(enable_gmm_colors, "Color", "Mono")

log_file       <- file.path(out_dir, paste0("Part1_", toupper(summary_method), "_TrueZero_", color_mode_str, "_Log.txt"))

file_plot_cv   <- file.path(out_dir, paste0("Part1_", toupper(summary_method), "_TrueZero_", color_mode_str, "_CV.tiff"))
file_plot_sd   <- file.path(out_dir, paste0("Part1_", toupper(summary_method), "_TrueZero_", color_mode_str, "_SD.tiff"))
file_plot_mean <- file.path(out_dir, paste0("Part1_", toupper(summary_method), "_TrueZero_", color_mode_str, "_Mean.tiff"))
file_plot_comb <- file.path(out_dir, paste0("Part1_", toupper(summary_method), "_TrueZero_", color_mode_str, "_Combined.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # ------------------------------------------------------------------- #
  # Section 1. Data Load 
  # ------------------------------------------------------------------- #
  log_msg("Loading Data...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # ------------------------------------------------------------------- #
  # Section 2. Advanced Aggregation (True Zero-Included Logic)
  # ------------------------------------------------------------------- #
  log_msg("Aggregating Data with True Zero-Included Logic...")
  
  df_merged_is <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  asv_is_summary <- df_merged_is %>%
    group_by(ASV) %>%
    summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
              .groups = "drop")
  
  all_sample_dates <- unique(df_ab$Date)
  
  asv_ab_summary <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    group_by(ASV) %>%
    summarise(
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE),
      asv_sd   = sd(Absolute_Abundance, na.rm = TRUE),
      agg_CV   = (asv_sd / asv_mean) * 100, 
      .groups = "drop"
    )
  
  asv_summary <- dplyr::inner_join(asv_is_summary, asv_ab_summary, by = "ASV") %>%
    dplyr::filter(!is.na(agg_IS) & !is.na(agg_CV) & !is.na(asv_sd) & asv_mean > 0)
  
  # ------------------------------------------------------------------- #
  # Section 3. Apply GMM Clustering & Validation
  # ------------------------------------------------------------------- #
  log_msg("Applying GMM Clustering and Validating Model...")
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_means <- gmm_model$parameters$mean
  cluster_order <- order(cluster_means) 
  mapped_clusters <- factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))
  asv_summary$Cluster <- mapped_clusters
  
  # 1) 실루엣 지수 (Silhouette Score) 계산
  dist_mat <- dist(asv_summary$agg_IS)
  sil_res <- cluster::silhouette(as.numeric(mapped_clusters), dist_mat)
  mean_sil_score <- mean(sil_res[, "sil_width"])
  
  # 2) 분류 불확실성 (Classification Uncertainty) 계산
  mean_uncertainty <- mean(gmm_model$uncertainty)
  
  # 3) 그룹 간 차이 검정 (ANOVA)
  aov_res <- aov(agg_IS ~ Cluster, data = asv_summary)
  pval_aov <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  pval_aov_str <- ifelse(pval_aov < 0.001, "p < 0.001", format.pval(pval_aov, digits = 3))
  
  log_msg("\n===================================================================")
  log_msg(" [Analysis Parameters & Clustering Info]")
  log_msg(sprintf(" - Plot Mode        : %s", ifelse(enable_gmm_colors, "Grouped (GMM Colors)", "Single Cluster (Mono Color)")))
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included (via tidyr::complete)"))
  log_msg(sprintf(" - Total ASVs Mapped: %d", nrow(asv_summary)))
  log_msg(sprintf(" - Outlier Filter   : %.1f IQR (Axis: %s)", iqr_multiplier, toupper(iqr_filter_mode)))
  log_msg("\n [GMM Clustering Validation (G = 3)]")
  log_msg(sprintf(" -> Mean Silhouette Score : %.3f (Target: > 0.5 for reasonable structure)", mean_sil_score))
  log_msg(sprintf(" -> Mean Uncertainty      : %.5f (Target: Close to 0 for clear boundaries)", mean_uncertainty))
  log_msg(sprintf(" -> ANOVA (Group Diff)    : F = %.2f, %s", summary(aov_res)[[1]][["F value"]][1], pval_aov_str))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 4. Independent IQR Filtering
  # ------------------------------------------------------------------- #
  filter_outliers_iqr_flexible <- function(df, x_var, y_var, k, mode) {
    df_filtered <- df
    if (mode %in% c("x", "both")) {
      q1x <- quantile(df[[x_var]], 0.25, na.rm = TRUE); q3x <- quantile(df[[x_var]], 0.75, na.rm = TRUE)
      lower_x <- q1x - k * (q3x - q1x); upper_x <- q3x + k * (q3x - q1x)
      df_filtered <- df_filtered %>% filter(.data[[x_var]] >= lower_x & .data[[x_var]] <= upper_x)
    }
    if (mode %in% c("y", "both")) {
      q1y <- quantile(df[[y_var]], 0.25, na.rm = TRUE); q3y <- quantile(df[[y_var]], 0.75, na.rm = TRUE)
      lower_y <- q1y - k * (q3y - q1y); upper_y <- q3y + k * (q3y - q1y)
      df_filtered <- df_filtered %>% filter(.data[[y_var]] >= lower_y & .data[[y_var]] <= upper_y)
    }
    return(df_filtered)
  }
  
  asv_filtered_cv   <- filter_outliers_iqr_flexible(asv_summary, "agg_IS", "agg_CV", iqr_multiplier, iqr_filter_mode)
  asv_filtered_sd   <- filter_outliers_iqr_flexible(asv_summary, "agg_IS", "asv_sd", iqr_multiplier, iqr_filter_mode)
  asv_filtered_mean <- filter_outliers_iqr_flexible(asv_summary, "agg_IS", "asv_mean", iqr_multiplier, iqr_filter_mode)
  
  # ------------------------------------------------------------------- #
  # Section 5. Fit Independent GAM Models & Log Statistics
  # ------------------------------------------------------------------- #
  log_msg("Fitting individual GAM models and extracting statistics...")
  
  extract_and_log_gam <- function(data, y_var, plot_title, plot_note) {
    fmla <- as.formula(paste(y_var, "~ s(agg_IS)"))
    gam_model <- mgcv::gam(fmla, data = data, method = "REML")
    gam_sum <- summary(gam_model)
    
    pval <- gam_sum$s.table[1, "p-value"]
    pval_str <- ifelse(pval < 0.001, "p < 0.001", format.pval(pval, digits = 3))
    
    log_msg(sprintf(" [%s]", plot_title))
    log_msg(sprintf("  -> Note: %s", plot_note))
    log_msg(sprintf("  -> GAM Stats: edf = %.2f, F = %.2f, %s, Deviance Explained = %.1f%%", 
                    gam_sum$s.table[1, "edf"], gam_sum$s.table[1, "F"], pval_str, gam_sum$dev.expl * 100))
    log_msg("")
    
    new_data <- data.frame(agg_IS = seq(min(data$agg_IS), max(data$agg_IS), length.out = 200))
    pred <- predict(gam_model, newdata = new_data, type = "response", se.fit = TRUE)
    new_data$fit <- pred$fit
    new_data$upr <- pred$fit + 1.96 * pred$se.fit
    new_data$lwr <- pred$fit - 1.96 * pred$se.fit
    
    return(new_data)
  }
  
  method_title <- tools::toTitleCase(summary_method)
  
  title_cv   <- paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year CV")
  title_sd   <- paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year SD")
  title_mean <- paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year Mean")
  
  note_cv   <- "Evaluates relative population volatility across the IS gradient (U-shape expected)."
  note_sd   <- "Evaluates absolute population volatility across the IS gradient."
  note_mean <- "Evaluates overall absolute abundance dominance across the IS gradient."
  
  log_msg("\n===================================================================")
  new_data_cv   <- extract_and_log_gam(asv_filtered_cv, "agg_CV", title_cv, note_cv)
  new_data_sd   <- extract_and_log_gam(asv_filtered_sd, "asv_sd", title_sd, note_sd)
  new_data_mean <- extract_and_log_gam(asv_filtered_mean, "asv_mean", title_mean, note_mean)
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 6. Visualization 
  # ------------------------------------------------------------------- #
  log_msg("Generating plot objects (Individual + Combined)...")
  
  trend_color   <- "#347433"  
  ci_color      <- "#347433"  
  
  base_theme <- theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "gray20", margin = margin(b = 10))
  )
  
  add_points <- function() {
    if (enable_gmm_colors) {
      custom_colors <- c("1_Negative" = "#0065F8", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
      list(
        geom_point(aes(color = Cluster), alpha = 0.8, size = 2),
        scale_color_manual(values = custom_colors)
      )
    } else {
      list(
        geom_point(color = "#000000", alpha = 0.8, size = 2)
      )
    }
  }
  
  p_cv <- ggplot(asv_filtered_cv, aes(x = agg_IS, y = agg_CV)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_cv, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_cv, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = title_cv, subtitle = note_cv, x = paste0(method_title, " Interaction Strength"), y = "CV of Absolute Abundance (%)") +
    theme(legend.position = ifelse(enable_gmm_colors, "top", "none")) + base_theme
  
  p_sd <- ggplot(asv_filtered_sd, aes(x = agg_IS, y = asv_sd)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_sd, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_sd, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = title_sd, subtitle = note_sd, x = paste0(method_title, " Interaction Strength"), y = "SD of Absolute Abundance") +
    theme(legend.position = "none") + base_theme
  
  p_mean <- ggplot(asv_filtered_mean, aes(x = agg_IS, y = asv_mean)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_mean, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_mean, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = title_mean, subtitle = note_mean, x = paste0(method_title, " Interaction Strength"), y = "Mean Absolute Abundance") +
    theme(legend.position = "none") + base_theme
  
  p_combined <- plot_grid(p_cv, p_sd, p_mean, ncol = 1, align = "v", rel_heights = c(1.1, 1, 1))
  
  # ------------------------------------------------------------------- #
  # Section 7. Export 
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    ggsave(file_plot_cv, plot = p_cv, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_sd, plot = p_sd, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_mean, plot = p_mean, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 10, height = 22, compression = "lzw")
    
    log_msg("\n[SUCCESS] 3 Individual plots and 1 Combined plot saved to disk successfully.")
  } else {
    log_msg("\n[SAFE MODE] Plots were generated in RStudio Viewer only (No files saved).")
  }
  
  print(p_combined)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

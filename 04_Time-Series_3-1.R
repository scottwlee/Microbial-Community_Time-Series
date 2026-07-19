################################################################################
# [Phase 5 - Part 4 (Version_2.1): Macro-ecological Traits (GMM Clustering - Expanded)]
# 목적: 10년 치 전체 데이터를 바탕으로 ASV의 온도 민감도(Temp_IS)에 따른
#       거시적 생태 특성(CV, SD, Mean)의 비선형적 관계를 분석함.
# 특징:
#   1) [자동 패키지 설치] 환경에 누락된 필수 라이브러리를 감지하고 자동 설치함.
#   2) Part 1, 2, 3의 분석을 통합하여 동일한 모집단(Core ASVs) 기반으로 통계적 오류 제거.
#   3) Y축 명칭에 'Absolute'를 명시하여 데이터의 정량적 본질을 학술적으로 엄밀히 표현.
#   4) USER SETTINGS의 [enable_gmm_colors] 스위치를 통해 3색 분리 플롯 생성.
#   5) [enable_save_outputs] 원터치 마스터 스위치로 파일 저장 여부 완벽 제어.
################################################################################

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 1) 필수 패키지 목록 정의 및 자동 설치
required_packages <- c("dplyr", "ggplot2", "mgcv", "cowplot", "mclust")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

# 2) 패키지 로드
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택 (옵션: "median", "mean")
summary_method  <- "median"  

# [2] 플롯 시각화 모드 선택 (옵션: TRUE, FALSE)
# - TRUE  : 호냉성, 범존종, 호열성 3가지 색상(#347433, #999999, #DC2525)으로 분리.
enable_gmm_colors <- TRUE 

# [3] GMM 클러스터 개수 (옵션: 정수, 기본값 3)
g_clusters      <- 3 

# [4] 제로 풍부도 배제 여부 (옵션: TRUE, FALSE)
exclude_zeros   <- FALSE     

# [5] 이상치 필터링 파라미터 (권장: "y"축 1.5 IQR 필터링)
iqr_multiplier  <- 1.5
iqr_filter_mode <- "y"

# [6] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 파일명 설정 (자동 적용)
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/04_GMM_Clustering")

# [원터치 제어] 마스터 스위치가 TRUE일 때만 디렉토리 생성
if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일명에 GMM 컬러 활성화 여부(Color/Mono)를 반영하여 덮어쓰기 방지
color_mode_str <- ifelse(enable_gmm_colors, "Color", "Mono")

log_file       <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", color_mode_str, "_Log.txt"))

file_plot_cv   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", color_mode_str, "_CV.tiff"))
file_plot_sd   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", color_mode_str, "_SD.tiff"))
file_plot_mean <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", color_mode_str, "_Mean.tiff"))
file_plot_comb <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", color_mode_str, "_Combined.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # -------------------------------------------------------------------
  # Section 1. Data Load & Merge 
  # -------------------------------------------------------------------
  log_msg("Loading and Merging Data...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  df_merged <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  if (exclude_zeros) {
    df_merged <- df_merged %>% filter(Absolute_Abundance > 0)
    zero_text <- "Zero-excluded"
  } else {
    zero_text <- "Zero-included"
  }
  
  # -------------------------------------------------------------------
  # Section 2. Aggregate Data & Apply GMM Clustering
  # -------------------------------------------------------------------
  log_msg("Aggregating data and running GMM...")
  asv_summary <- df_merged %>%
    group_by(ASV) %>%
    summarise(
      agg_IS   = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE),
      asv_sd   = sd(Absolute_Abundance, na.rm = TRUE),
      agg_CV   = (asv_sd / asv_mean) * 100
    ) %>% 
    filter(!is.na(agg_IS) & !is.na(agg_CV) & !is.na(asv_sd) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_means <- gmm_model$parameters$mean
  cluster_order <- order(cluster_means) 
  mapped_clusters <- factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))
  
  asv_summary$Cluster <- mapped_clusters
  
  log_msg("\n===================================================================")
  log_msg(" [Analysis Parameters & Clustering Info]")
  log_msg(sprintf(" - Plot Mode        : %s", ifelse(enable_gmm_colors, "Grouped (GMM Colors)", "Single Cluster (Mono Color)")))
  log_msg(sprintf(" - Output Save Mode : %s", ifelse(enable_save_outputs, "ENABLED (Writing to disk)", "DISABLED (Dry-run mode)")))
  log_msg(sprintf(" - GMM Clusters     : %d", g_clusters))
  log_msg(sprintf(" - Zero Abundance   : %s", zero_text))
  log_msg(sprintf(" - Outlier Filter   : %.1f IQR (Axis: %s)", iqr_multiplier, toupper(iqr_filter_mode)))
  log_msg("===================================================================\n")
  
  # -------------------------------------------------------------------
  # Section 3. Independent IQR Filtering
  # -------------------------------------------------------------------
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
  
  # -------------------------------------------------------------------
  # Section 4. Fit Independent GAM Models & Log Statistics
  # -------------------------------------------------------------------
  log_msg("Fitting individual GAM models and extracting statistics...")
  
  extract_and_log_gam <- function(data, y_var, name_str) {
    fmla <- as.formula(paste(y_var, "~ s(agg_IS)"))
    gam_model <- mgcv::gam(fmla, data = data, method = "REML")
    gam_sum <- summary(gam_model)
    
    pval <- gam_sum$s.table[1, "p-value"]
    pval_str <- ifelse(pval < 0.001, "p < 0.001", format.pval(pval, digits = 3))
    
    log_msg(sprintf(" [%s GAM Stats] edf = %.2f, F = %.2f, %s, Deviance Explained = %.1f%%", 
                    name_str, gam_sum$s.table[1, "edf"], gam_sum$s.table[1, "F"], pval_str, gam_sum$dev.expl * 100))
    
    new_data <- data.frame(agg_IS = seq(min(data$agg_IS), max(data$agg_IS), length.out = 200))
    pred <- predict(gam_model, newdata = new_data, type = "response", se.fit = TRUE)
    new_data$fit <- pred$fit
    new_data$upr <- pred$fit + 1.96 * pred$se.fit
    new_data$lwr <- pred$fit - 1.96 * pred$se.fit
    
    return(new_data)
  }
  
  log_msg("\n===================================================================")
  new_data_cv   <- extract_and_log_gam(asv_filtered_cv, "agg_CV", "CV")
  new_data_sd   <- extract_and_log_gam(asv_filtered_sd, "asv_sd", "SD")
  new_data_mean <- extract_and_log_gam(asv_filtered_mean, "asv_mean", "Mean")
  log_msg("===================================================================\n")
  
  # -------------------------------------------------------------------
  # Section 5. Visualization (Switchable Logic & Absolute Y-Axis Labels)
  # -------------------------------------------------------------------
  log_msg("Generating individual plot objects (p_cv, p_sd, p_mean)...")
  
  method_title <- tools::toTitleCase(summary_method)
  
  trend_color   <- "#0065F8"  
  ci_color      <- "#0065F8"  
  
  base_theme <- theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "gray20", margin = margin(b = 10))
  )
  
  add_points <- function() {
    if (enable_gmm_colors) {
      custom_colors <- c("1_Negative" = "#347433", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
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
  
  # 1) Plot: IS vs CV [Y축: CV of Absolute Abundance (%)]
  p_cv <- ggplot(asv_filtered_cv, aes(x = agg_IS, y = agg_CV)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_cv, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_cv, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year CV"),
         subtitle = "Relative Population Volatility",
         x = paste0(method_title, " Interaction Strength"), y = "CV of Absolute Abundance (%)") +
    theme(legend.position = ifelse(enable_gmm_colors, "top", "none")) + base_theme
  
  # 2) Plot: IS vs SD [Y축: SD of Absolute Abundance]
  p_sd <- ggplot(asv_filtered_sd, aes(x = agg_IS, y = asv_sd)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_sd, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_sd, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year SD"),
         subtitle = "Absolute Population Volatility",
         x = paste0(method_title, " Interaction Strength"), y = "SD of Absolute Abundance") +
    theme(legend.position = "none") + base_theme
  
  # 3) Plot: IS vs Mean [Y축: Mean Absolute Abundance]
  p_mean <- ggplot(asv_filtered_mean, aes(x = agg_IS, y = asv_mean)) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) + 
    add_points() +
    geom_ribbon(data = new_data_mean, aes(x = agg_IS, ymin = lwr, ymax = upr), alpha = 0.15, fill = ci_color, inherit.aes = FALSE) +
    geom_line(data = new_data_mean, aes(x = agg_IS, y = fit), color = trend_color, linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = paste0("Macro-ecological Traits: ", method_title, " IS vs 10-Year Mean"),
         subtitle = "Overall Abundance Dominance",
         x = paste0(method_title, " Interaction Strength"), y = "Mean Absolute Abundance") +
    theme(legend.position = "none") + base_theme
  
  # 세 플롯 상하 3단 결합 
  p_combined <- plot_grid(p_cv, p_sd, p_mean, ncol = 1, align = "v", rel_heights = c(1.1, 1, 1))
  
  # -------------------------------------------------------------------
  # Section 6. Export (원터치 제어)
  # -------------------------------------------------------------------
  if (enable_save_outputs) {
    ggsave(file_plot_cv, plot = p_cv, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_sd, plot = p_sd, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_mean, plot = p_mean, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 10, height = 22, compression = "lzw")
    log_msg("\n[SUCCESS] All plots saved to disk successfully.")
  } else {
    log_msg("\n[SAFE MODE] Plots were generated in RStudio Viewer only (No files saved).")
  }
  
  # RStudio Viewer 출력
  print(p_combined)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

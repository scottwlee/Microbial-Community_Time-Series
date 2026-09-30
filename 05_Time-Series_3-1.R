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
#   2) [논문 방어용 자동 통계 라우팅] 등분산성(Levene's test) 및 정규성(Shapiro-Wilk) 
#      가정 검정을 수행한 후, 조건에 따라 ANOVA 또는 Welch's ANOVA를 자동 선택.
#   3) [로깅 고도화] 포맷팅 에러 방지를 위해 paste0 기반의 직관적 로깅 적용.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "cluster", "car")
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
  library(car)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
summary_method  <- "median"  
enable_gmm_colors <- TRUE 
g_clusters      <- 3 
iqr_multiplier  <- 1.5
iqr_filter_mode <- "y"
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 파일명 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_V2_Output/01_Data_Integration")
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
  # Section 3. Apply GMM Clustering & Advanced Validation Routing
  # ------------------------------------------------------------------- #
  log_msg("Applying GMM Clustering and Validating Model...");
  set.seed(414);
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters);
  cluster_means <- gmm_model$parameters$mean;
  cluster_order <- order(cluster_means); 
  mapped_clusters <- factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"));
  asv_summary$Cluster <- mapped_clusters;
  
  dist_mat <- dist(asv_summary$agg_IS);
  sil_res <- cluster::silhouette(as.numeric(mapped_clusters), dist_mat);
  mean_sil_score <- mean(sil_res[, "sil_width"]);
  mean_uncertainty <- mean(gmm_model$uncertainty);      base_aov <- aov(agg_IS ~ Cluster, data = asv_summary);   aov_resid <- residuals(base_aov);   if(length(aov_resid) > 5000) { aov_resid <- sample(aov_resid, 5000) };   shapiro_res <- shapiro.test(aov_resid);      levene_res <- car::leveneTest(agg_IS ~ Cluster, data = asv_summary);   levene_pval <- levene_res$`Pr(>F)`[1];
  
  kw_res <- kruskal.test(agg_IS ~ Cluster, data = asv_summary);
  
  if (levene_pval >= 0.05) {
    selected_test_name <- "Standard One-way ANOVA";
    final_stat <- summary(base_aov)[[1]][["F value"]][1];
    final_pval <- summary(base_aov)[[1]][["Pr(>F)"]][1];
    stat_label <- "F";
  } else {
    selected_test_name <- "Welch's ANOVA (Variance Heterogeneity Adjusted)";
    welch_res <- oneway.test(agg_IS ~ Cluster, data = asv_summary, var.equal = FALSE);
    final_stat <- welch_res$statistic;
    final_pval <- welch_res$p.value;
    stat_label <- "F";
  };
  
  format_pval <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 4)));
  
  log_msg("\n===================================================================");
  log_msg(" [Analysis Parameters & Clustering Info]");
  log_msg(paste0(" - Plot Mode        : ", ifelse(enable_gmm_colors, "Grouped (GMM Colors)", "Single Cluster (Mono Color)")));
  log_msg(paste0(" - Summary Method   : ", toupper(summary_method)));
  log_msg(" - Zero Handling    : True Zero-Included (via tidyr::complete)");
  log_msg(paste0(" - Total ASVs Mapped: ", nrow(asv_summary)));
  log_msg(paste0(" - Outlier Filter   : ", iqr_multiplier, " IQR (Axis: ", toupper(iqr_filter_mode), ")"));
  
  log_msg("\n [GMM Clustering Validation (G = 3)]");
  log_msg(paste0(" -> Mean Silhouette Score : ", round(mean_sil_score, 3), " (Target: > 0.5 for reasonable structure)"));
  log_msg(paste0(" -> Mean Uncertainty      : ", round(mean_uncertainty, 5), " (Target: Close to 0 for clear boundaries)"));
  
  log_msg("\n [Statistical Assumption Tests]");
  log_msg(paste0(" -> Normality (Shapiro-Wilk)       : W = ", round(shapiro_res$statistic, 3), ", ", format_pval(shapiro_res$p.value)));   log_msg(paste0(" -> Homoscedasticity (Levene's)    : F = ", round(levene_res$`F value`[1], 3), ", ", format_pval(levene_pval)));
  
  log_msg("\n [Selected Group Difference Test]");
  log_msg(paste0(" -> Applied Test Model             : ", selected_test_name));
  log_msg(paste0(" -> Main Result                    : ", stat_label, " = ", round(final_stat, 2), ", ", format_pval(final_pval)));
  log_msg(paste0(" -> Non-parametric Reference (KW)  : Chi-squared = ", round(kw_res$statistic, 2), ", ", format_pval(kw_res$p.value)));
  log_msg("===================================================================\n");
  
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
    
    log_msg(paste0(" [", plot_title, "]"))
    log_msg(paste0("  -> Note: ", plot_note))
    log_msg(paste0("  -> GAM Stats: edf = ", round(gam_sum$s.table[1, "edf"], 2), ", F = ", round(gam_sum$s.table[1, "F"], 2), ", ", pval_str, ", Deviance Explained = ", round(gam_sum$dev.expl * 100, 1), "%"))
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

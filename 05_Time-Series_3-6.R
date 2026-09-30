# [Github: 04_Time-Series_3-6.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 6 (Master Version): Auto-Linked TOW & Thermophilic Expansion]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 11 & 13"에서 통합, 이관 및 재정비된 코드입니다.)
#
# 목적: Tukey's HSD를 통해 호열성 군집의 최전성기(TOW)를 추출하고, 이를 
#       시차 보정(Phase-Aligned) 팽창 분석으로 자동 연동함.
# 특징:
#   1) [생물학적 대조군] Part 4 (호냉성 붕괴)와 완벽한 대칭(거울상)을 이룸.
#   2) [논문 방어 검증] 정규성, 등분산성 검정 및 비모수 교차검증(KW) 자동 기록.
#   3) [정밀 True Zero] 시계열 왜곡을 막기 위한 결측치 0 복원 엔진 100% 가동.
#   4) [로깅 고도화] 포맷팅 에러 방지를 위해 세미콜론 및 paste0 기반 로깅 적용.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "stats", "agricolae", "car")
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
  library(stats)
  library(agricolae) 
  library(car)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
summary_method <- "median"  
g_clusters     <- 3 
target_cluster <- "3_Positive" 

peak_selection_method <- "tukey"  
peak_k_clusters       <- 3  
enable_save_outputs   <- TRUE

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정 (논문 전반 적용 대원칙)
# ------------------------------------------------------------------- #
pos_base_color <- "#DC2525" 
pos_fill       <- "#FFCCCC" 
trend_color    <- "#347433" 

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_V2_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
out_dir    <- file.path(base_dir, "05_Phase5_Output/06_Auto_Linked_TOW_Expansion")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

log_file       <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_", toupper(peak_selection_method), "_Integrated_Log.txt"))

f_plot_env     <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_1_Indiv_Env.tiff"))
f_plot_pheno   <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_2_Indiv_Phenology.tiff"))
f_plot_anova   <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_3_Indiv_ANOVA.tiff"))
f_plot_exp     <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_4_Indiv_Expansion.tiff"))
f_plot_therm   <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_5_Indiv_ThermalLimitation.tiff"))

f_plot_comb_tow <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_Combined_TOW.tiff"))
f_plot_comb_im  <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_Combined_Impact.tiff"))
f_csv_impact    <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_Impact_Data.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

format_pval <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 4)))

tryCatch({
  log_msg("\n===================================================================");
  log_msg(" [Phase 5 - Part 6: Auto-Linked TOW & Thermophilic Expansion Pipeline]");
  log_msg(" [Analysis Parameters & Settings]");
  log_msg(paste0(" - Target Cluster        : ", target_cluster, " (Biological Control)"));
  log_msg(paste0(" - Summary Method        : ", toupper(summary_method)));
  log_msg(paste0(" - Peak Selection Method : ", toupper(peak_selection_method)));
  log_msg(" - Zero Handling         : True Zero-Included");
  log_msg("===================================================================\n");
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Common GMM Clustering 
  # ------------------------------------------------------------------- #
  log_msg("Step 1. Loading Data & Applying True Zero-Included GMM...");
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE)
  
  colnames(df_smap_raw)[1] <- "ASV_Hash"
  asv_dict <- df_ab_raw %>% dplyr::select(ASV_Hash, ASV_ID) %>% dplyr::distinct()
  df_smap <- df_smap_raw %>% dplyr::inner_join(asv_dict, by = "ASV_Hash") %>% dplyr::select(ASV = ASV_ID, Best_TP)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  asv_is_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE), .groups = "drop")
  
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
  
  # 에러 방지 처리 완료 (백슬래시 제거)
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)

  # ------------------------------------------------------------------- #
  # Section 2. Peak Window (TOW) Definition & Assumption Tests
  # ------------------------------------------------------------------- #
  log_msg("\n>>> STARTING TOW DEFINITION <<<");
  
  df_env <- df_ab_raw %>% 
    dplyr::mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(!is.na(Temperature)) %>%
    dplyr::group_by(Date, Month) %>%
    dplyr::summarise(Temperature = mean(Temperature, na.rm = TRUE), .groups = "drop")
  
  sampled_ym <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Year, Month) %>% dplyr::distinct()
  
  df_cluster_abund <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    dplyr::filter(Cluster == target_cluster) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  df_cluster_monthly <- sampled_ym %>%
    dplyr::left_join(df_cluster_abund, by = c("Year", "Month")) %>%
    dplyr::mutate(
      Total_Abund = tidyr::replace_na(Total_Abund, 0),
      Log_Abundance = log10(Total_Abund + 1),
      Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
      Year_Fct = factor(Year)
    )
  
  target_months <- c()
  assumption_log <- ""
  
  if (peak_selection_method == "tukey") {
    aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
    
    # 1) 가정 검정 (Assumption Tests)
    aov_resid <- residuals(aov_model)
    if(length(aov_resid) > 5000) { aov_resid <- sample(aov_resid, 5000) }
    shapiro_res <- shapiro.test(aov_resid)
    
    levene_res <- car::leveneTest(Log_Abundance ~ Month_Fct, data = df_cluster_monthly)
    levene_pval <- levene_res$`Pr(>F)`[1]
    
    # 2) 비모수 교차검증 (Kruskal-Wallis)
    kw_res <- kruskal.test(Log_Abundance ~ Month_Fct, data = df_cluster_monthly)
    
    # Tukey HSD 진행
    tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
    top_group <- tukey_res$groups
    
    tow_months_str <- rownames(top_group)[grepl("a", top_group$groups)]
    target_months <- which(month.abb %in% tow_months_str)
    
    f_val <- summary(aov_model)[[1]][["F value"]][1]
    pval_anova <- summary(aov_model)[[1]][["Pr(>F)"]][1]
    
  } else if (peak_selection_method == "kmeans") {
    monthly_means <- df_cluster_monthly %>% dplyr::group_by(Month, Month_Fct) %>% dplyr::summarise(mean_abund = mean(Log_Abundance), .groups="drop")
    set.seed(414)
    km_res <- kmeans(monthly_means$mean_abund, centers = peak_k_clusters)
    monthly_means$KM_Cluster <- km_res$cluster
    target_months <- monthly_means %>% dplyr::filter(KM_Cluster == which.max(km_res$centers)) %>% dplyr::pull(Month) %>% sort()
    
    aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
    f_val <- summary(aov_model)[[1]][["F value"]][1]
    pval_anova <- summary(aov_model)[[1]][["Pr(>F)"]][1]
  }
  
  target_months <- sort(target_months)
  tow_months_str <- paste(month.abb[target_months], collapse = ", ")
  
  tow_label <- ifelse(length(target_months) == 1, paste0("TOW (", month.abb[target_months], ")"), 
                      ifelse(all(diff(target_months) == 1), paste0("TOW (", month.abb[min(target_months)], "-", month.abb[max(target_months)], ")"), 
                             paste0("TOW (", tow_months_str, ")")))
  
  df_cluster_monthly <- df_cluster_monthly %>% dplyr::mutate(Is_Peak = ifelse(Month %in% target_months, tow_label, "Non-Peak"))
  
  pval_text_tow <- format_pval(pval_anova)
  
  title_C <- paste0("C. Statistical Validation of ", tow_label);
  note_C  <- paste0("Thermophilic Optimal Window (TOW) defined via ", toupper(peak_selection_method), " on True Zero-Included data.");
  
  log_msg("\n-------------------------------------------------------------------");
  log_msg(paste0("[Plot C: ", title_C, "]"));
  log_msg(paste0(" -> Note                           : ", note_C));
  log_msg(paste0(" -> Two-way ANOVA (Block: Year)    : F = ", round(f_val, 2), ", ", pval_text_tow));
  
  if (peak_selection_method == "tukey") {
    log_msg(" [Statistical Assumption Tests]");
    log_msg(paste0(" -> Normality (Shapiro-Wilk)       : W = ", round(shapiro_res$statistic, 3), ", ", format_pval(shapiro_res$p.value)));     log_msg(paste0(" -> Homoscedasticity (Levene's)    : F = ", round(levene_res$`F value`[1], 3), ", ", format_pval(levene_pval)));
    log_msg(paste0(" -> Non-parametric Reference (KW)  : Chi-squared = ", round(kw_res$statistic, 2), ", ", format_pval(kw_res$p.value)));
  }
  log_msg("-------------------------------------------------------------------\n");
  
  # ------------------------------------------------------------------- #
  # Section 3. Plotting TOW Definition (Individual & Combined)
  # ------------------------------------------------------------------- #
  p11_temp <- ggplot(df_env, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = trend_color, fill = trend_color, alpha = 0.2, linewidth=1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(add = 0.5)) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p11_abund <- ggplot(df_cluster_monthly, aes(x = Month, y = Log_Abundance)) +
    geom_jitter(color = pos_base_color, alpha = 0.5, size = 1.5, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = trend_color, fill = trend_color, alpha = 0.2, linewidth = 1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(add = 0.5)) +
    labs(title = paste0("B. Biological Phenology (", target_cluster, ")"), x = "Month", y = "Log10(Cluster Abundance + 1)") + 
    theme(plot.title = element_text(face = "bold"))
  
  fill_palette <- c(pos_fill, "gray80"); names(fill_palette) <- c(tow_label, "Non-Peak")
  
  p11_anova <- ggplot(df_cluster_monthly, aes(x = Month_Fct, y = Log_Abundance, fill = Is_Peak)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
    geom_jitter(color = pos_base_color, width = 0.15, alpha = 0.6, size = 1.2) +
    scale_fill_manual(values = fill_palette) + 
    labs(title = title_C, subtitle = paste0("Two-way ANOVA (Month Effect): ", pval_text_tow), x = "Month of the Year", y = "Log10(Cluster Abundance + 1)") +
    theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(face = "bold")) 
  
  p11_combined <- plot_grid(plot_grid(p11_temp, p11_abund, ncol = 1, align = "v"), p11_anova, ncol = 2, rel_widths = c(1, 1.2))
  
  # ------------------------------------------------------------------- #
  # Section 4. Linkage & Lag-Adjusted Thermophilic Expansion
  # ------------------------------------------------------------------- #
  tow_months <- target_months 
  log_msg(paste0("\n>>> [LINKAGE BRIDGE] TOW (Months ", paste(tow_months, collapse=","), ") automatically passed to Expansion Section <<<"));
  log_msg(">>> STARTING LAG-ADJUSTED THERMOPHILIC EXPANSION <<<\n");
  
  df_env_daily <- df_ab_raw %>% dplyr::select(Date = Sample_Date, Temperature) %>% dplyr::distinct() %>% dplyr::filter(!is.na(Temperature)) %>% dplyr::arrange(Date)
  
  df_aligned <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::arrange(ASV, as.Date(Date)) %>% 
    dplyr::inner_join(asv_clusters, by = "ASV") %>% 
    dplyr::inner_join(df_smap, by = "ASV") %>%
    dplyr::mutate(Response_Lag_Weeks = abs(Best_TP)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::mutate(shift_n = as.integer(Response_Lag_Weeks[1]), Shifted_Abundance = dplyr::lead(Absolute_Abundance, n = shift_n[1])) %>%
    dplyr::ungroup() %>% dplyr::filter(!is.na(Shifted_Abundance)) %>%
    dplyr::inner_join(df_env_daily, by = "Date") %>% dplyr::mutate(Date = as.Date(Date))
  
  analysis_df <- df_aligned %>%
    dplyr::mutate(Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(Month %in% tow_months & Cluster == target_cluster) %>% 
    dplyr::group_by(Date, Year, Temperature) %>%
    dplyr::summarise(Total_Aligned_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Total_Aligned_Abund + 1), Year_Fct = as.factor(Year))
  
  lm_model <- lm(Log_Abund ~ Year, data = analysis_df)
  lm_summary <- summary(lm_model)
  pval_abund <- lm_summary$coefficients[2, 4]
  slope_abund <- lm_summary$coefficients[2, 1]
  
  gam_model <- mgcv::gam(Log_Abund ~ s(Temperature, k=6), data = analysis_df, method="REML")
  pval_direct <- summary(gam_model)$s.table[1, 4]
  
  title_exp   <- "D. Expansion of Thermophilic Engine";
  title_therm <- "E. High-Resolution TOW Thermal Limitation";
  note_impact <- paste0("Analyzed strictly within the Dynamically Linked TOW (Months ", paste(tow_months, collapse="-"), "). Abundances are phase-aligned using True Zero-Included Best_TP.");
  
  log_msg("-------------------------------------------------------------------");
  log_msg(paste0("[Plot D: ", title_exp, "]"));
  log_msg(paste0(" -> Note                     : ", note_impact));
  log_msg(paste0(" -> Linear Regression Slope  : ", round(slope_abund, 4)));
  log_msg(paste0(" -> P-value                  : ", format_pval(pval_abund)));
  log_msg("");
  log_msg(paste0("[Plot E: ", title_therm, "]"));
  log_msg(paste0(" -> Note                     : ", note_impact));
  log_msg(paste0(" -> GAM P-value              : ", format_pval(pval_direct)));
  log_msg("-------------------------------------------------------------------\n");
  
  # ------------------------------------------------------------------- #
  # Section 5. Plotting Climate Impact (Individual & Combined)
  # ------------------------------------------------------------------- #
  subtitle_exp <- paste0("Rise of Phase-Aligned Abundance (Slope = ", round(slope_abund, 3), ", ", format_pval(pval_abund), ")")
  
  p13_abund <- ggplot(analysis_df, aes(x = Year, y = Log_Abund)) +
    geom_boxplot(aes(group=Year_Fct), fill=pos_fill, alpha=0.3, color=pos_base_color, outlier.shape=NA) +
    geom_jitter(color=pos_base_color, alpha=0.6, width=0.15) + 
    geom_smooth(method="lm", color=trend_color, fill=trend_color, alpha=0.2, linewidth=1.2) +
    labs(title=title_exp, subtitle=subtitle_exp, x="Year", y="Log10(Lag-Adjusted TOW Abundance + 1)") +
    scale_x_continuous(breaks = min(analysis_df$Year):max(analysis_df$Year), expand = expansion(add = 0.5)) + 
    theme(plot.title = element_text(face="bold"))
  
  subtitle_therm <- paste0("GAM correlation between Trigger Temp & Abundance (", format_pval(pval_direct), ")")
  
  p13_direct <- ggplot(analysis_df, aes(x = Temperature, y = Log_Abund)) +
    geom_point(color=pos_base_color, size=2.5, alpha=0.6) + 
    geom_smooth(method="gam", formula=y~s(x, bs="cs", k=6), color=trend_color, fill=trend_color, alpha=0.2, linetype="solid", linewidth=1.2) +
    labs(title=title_therm, subtitle=subtitle_therm, x="TOW Trigger Temp (°C)", y="Log10(Lag-Adjusted TOW Abundance + 1)") +
    theme(plot.title = element_text(face="bold"))
  
  p13_combined <- plot_grid(p13_abund, p13_direct, ncol=2, align="h")
  
  # ------------------------------------------------------------------- #
  # Section 6. Export (All Variations)
  # ------------------------------------------------------------------- #
  if(enable_save_outputs) { 
    ggsave(f_plot_env,   plot = p11_temp,   device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_pheno, plot = p11_abund,  device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_anova, plot = p11_anova,  device = "tiff", dpi = 600, width = 8, height = 8, compression = "lzw")
    ggsave(f_plot_exp,   plot = p13_abund,  device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_therm, plot = p13_direct, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    
    ggsave(f_plot_comb_tow, plot = p11_combined, device = "tiff", dpi = 600, width = 16, height = 8, compression = "lzw") 
    ggsave(f_plot_comb_im,  plot = p13_combined, device = "tiff", dpi = 600, width = 13, height = 6, compression = "lzw") 
    
    write.csv(analysis_df, f_csv_impact, row.names = FALSE)
    
    log_msg("[SUCCESS] 5 Individual Plots, 2 Combined Plots, and CSV Data saved successfully.")
  } else {
    log_msg("[SAFE MODE] Plots generated in Viewer (No files saved).")
  }
  
  print(p11_combined)
  print(p13_combined)
  
}, error = function(e) { log_msg(paste0("ERROR: ", e$message)); stop(e) })

##### END. ######################################################################

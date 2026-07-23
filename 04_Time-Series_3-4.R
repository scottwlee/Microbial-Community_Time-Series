# [Github: 04_Time-Series_3-4.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 4 (Master Version): Auto-Linked POW & Climate Impact]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 7 & 9"에서 통합, 이관 및 재정비된 코드입니다.)
#
# 목적: Tukey's HSD 그룹핑(Letter Grouping)을 통해 통계적으로 완벽하게 정의된 
#       최전성기(POW)를 추출하고, 이를 시차 보정 억제 분석으로 자동 연동함.
# 특징:
#   1) [Tukey HSD 그룹핑] 'a' 그룹(최고 생물량과 통계적 차이가 없는 달) 자동 추출.
#   2) [Two-way ANOVA] 10년 시계열 특성을 반영하여 연도(Year)를 블록 효과로 통제.
#   3) [로깅 고도화] 단일 통합 로그파일에 분석 파라미터, 플롯 Title, Note, 통계치 기록.
#   4) [다중 플롯 분리] 통합 패널 플롯 2종 외에 5개의 개별 플롯 모두 독립적으로 저장.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "stats", "agricolae")
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
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################

# [1] 공통 GMM 클러스터 설정
# - summary_method: 온도 민감도(IS) 요약 대푯값 기준 (옵션: "median"(권장), "mean")
# - g_clusters    : 군집 분할 수 (권장: 3)
# - target_cluster: POW 및 기후 변화 억제 효과를 검증할 메인 타겟 군집 (옵션: "1_Negative", "2_Neutral", "3_Positive")
summary_method <- "median"  
g_clusters     <- 3 
target_cluster <- "1_Negative" 

# [2] 생태적 최전성기(Peak Window) 선별 방식 설정
# - peak_selection_method: 최전성기를 정의하는 통계 알고리즘.
#   * "tukey": 연도 블록 효과를 통제한 ANOVA 후 Tukey HSD 사후검정으로 최상위 그룹 추출 (논문 방어용 권장)
#   * "kmeans": 월별 궤적을 순수 기계학습(거리 기반)으로 분할
# - peak_k_clusters: (kmeans 선택 시) 분할할 계층 수 (권장: 3)
peak_selection_method <- "tukey"  
peak_k_clusters       <- 3  

# [3] 결과 저장 마스터 스위치
# - 옵션: TRUE (폴더에 플롯과 로그 자동 저장), FALSE (뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정
# ------------------------------------------------------------------- #
neg_base_color <- "#0065F8" # Negative 데이터 포인트 및 박스플롯 테두리 (Blue)
neg_fill       <- "#99C2FF" # Negative 박스플롯 내부 배경색 (Light Blue)
trend_color    <- "#347433" # 생물학적 & 물리적 추세선 및 CI (Green)

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 개편: Part 4)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")

# 통합 디렉토리 설정
out_dir    <- file.path(base_dir, "05_Phase5_Output/04_Auto_Linked_POW_Impact")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

# 파일명 접두사 Part4_ 적용 및 단일 통합 로그 생성
log_file <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_", toupper(peak_selection_method), "_Integrated_Log.txt"))

# 플롯 파일명 지정 (개별 및 결합)
f_plot_env     <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_1_Indiv_Env.tiff"))
f_plot_pheno   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_2_Indiv_Phenology.tiff"))
f_plot_anova   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_3_Indiv_ANOVA.tiff"))
f_plot_supp    <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_4_Indiv_Suppression.tiff"))
f_plot_therm   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_5_Indiv_ThermalLimitation.tiff"))

f_plot_comb_pw <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_Combined_PeakWindow.tiff"))
f_plot_comb_im <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_Combined_Impact.tiff"))
f_csv_impact   <- file.path(out_dir, paste0("Part4_", toupper(summary_method), "_Impact_Data.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 4: Auto-Linked POW & Climate Impact Pipeline]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Target Cluster        : %s", target_cluster))
  log_msg(sprintf(" - Summary Method        : %s", toupper(summary_method)))
  log_msg(sprintf(" - Peak Selection Method : %s", toupper(peak_selection_method)))
  log_msg(sprintf(" - Zero Handling         : True Zero-Included"))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Common GMM Clustering 
  # ------------------------------------------------------------------- #
  log_msg("Step 1. Loading Data & Applying True Zero-Included GMM...")
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
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)
  
  # ------------------------------------------------------------------- #
  # Section 2. Peak Window Definition 
  # ------------------------------------------------------------------- #
  log_msg("\n>>> STARTING PEAK WINDOW DEFINITION <<<")
  
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
  
  if (peak_selection_method == "tukey") {
    aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
    tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
    top_group <- tukey_res$groups
    
    peak_months_str <- rownames(top_group)[grepl("a", top_group$groups)]
    target_months <- which(month.abb %in% peak_months_str)
    
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
  peak_months_str <- paste(month.abb[target_months], collapse = ", ")
  
  peak_label <- ifelse(length(target_months) == 1, paste0("Peak Window (", month.abb[target_months], ")"), 
                       ifelse(all(diff(target_months) == 1), paste0("Peak Window (", month.abb[min(target_months)], "-", month.abb[max(target_months)], ")"), 
                              paste0("Peak Window (", peak_months_str, ")")))
  
  df_cluster_monthly <- df_cluster_monthly %>% dplyr::mutate(Is_Peak = ifelse(Month %in% target_months, peak_label, "Non-Peak"))
  
  pval_text_p7 <- ifelse(pval_anova < 0.001, "p < 0.001", sprintf("p = %.3f", pval_anova))
  
  title_C <- paste("C. Statistical Validation of", peak_label)
  note_C  <- sprintf("Peak Window defined via %s on True Zero-Included data.", toupper(peak_selection_method))
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg(sprintf("[Plot C: %s]", title_C))
  log_msg(sprintf(" -> Note                     : %s", note_C))
  log_msg(sprintf(" -> Two-way ANOVA (Block: Year) F-value : %.2f", f_val))
  log_msg(sprintf(" -> P-value                             : %s", pval_text_p7))
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 3. Plotting Peak Window Definition (Individual & Combined)
  # ------------------------------------------------------------------- #
  p7_temp <- ggplot(df_env, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = trend_color, fill = trend_color, alpha = 0.2, linewidth=1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p7_abund <- ggplot(df_cluster_monthly, aes(x = Month, y = Log_Abundance)) +
    geom_jitter(color = neg_base_color, alpha = 0.5, size = 1.5, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = trend_color, fill = trend_color, alpha = 0.2, linewidth = 1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = paste0("B. Biological Phenology (", target_cluster, ")"), x = "Month", y = "Log10(Cluster Abundance + 1)") + 
    theme(plot.title = element_text(face = "bold"))
  
  fill_palette <- c(neg_fill, "gray80"); names(fill_palette) <- c(peak_label, "Non-Peak")
  
  p7_anova <- ggplot(df_cluster_monthly, aes(x = Month_Fct, y = Log_Abundance, fill = Is_Peak)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
    geom_jitter(color = neg_base_color, width = 0.15, alpha = 0.6, size = 1.2) +
    scale_fill_manual(values = fill_palette) + 
    labs(title = title_C, subtitle = paste("Two-way ANOVA (Month Effect):", pval_text_p7), x = "Month of the Year", y = "Log10(Cluster Abundance + 1)") +
    theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(face = "bold")) 
  
  p7_combined <- plot_grid(plot_grid(p7_temp, p7_abund, ncol = 1, align = "v"), p7_anova, ncol = 2, rel_widths = c(1, 1.2))
  
  # ------------------------------------------------------------------- #
  # Section 4. Linkage & Lag-Adjusted Climate Impact
  # ------------------------------------------------------------------- #
  eco_window_months <- target_months 
  log_msg(sprintf("\n>>> [LINKAGE BRIDGE] POW (%s) automatically passed to Climate Impact Section <<<", paste(eco_window_months, collapse=",")))
  log_msg(">>> STARTING LAG-ADJUSTED CLIMATE IMPACT <<<\n")
  
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
    dplyr::filter(Month %in% eco_window_months & Cluster == target_cluster) %>% 
    dplyr::group_by(Date, Year, Temperature) %>%
    dplyr::summarise(Total_Aligned_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Total_Aligned_Abund + 1), Year_Fct = as.factor(Year))
  
  lm_model <- lm(Log_Abund ~ Year, data = analysis_df)
  lm_summary <- summary(lm_model)
  pval_abund <- lm_summary$coefficients[2, 4]
  slope_abund <- lm_summary$coefficients[2, 1]
  
  gam_model <- mgcv::gam(Log_Abund ~ s(Temperature, k=6), data = analysis_df, method="REML")
  pval_direct <- summary(gam_model)$s.table[1, 4]
  
  format_pval <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  
  title_supp <- "D. Suppression of Psychrophilic Engine"
  title_therm <- "E. High-Resolution POW Thermal Limitation"
  note_impact <- sprintf("Analyzed strictly within the Dynamically Linked POW (Months %s). Abundances are phase-aligned using True Zero-Included Best_TP.", paste(eco_window_months, collapse="-"))
  
  log_msg("-------------------------------------------------------------------")
  log_msg(sprintf("[Plot D: %s]", title_supp))
  log_msg(sprintf(" -> Note                     : %s", note_impact))
  log_msg(sprintf(" -> Linear Regression Slope  : %.4f", slope_abund))
  log_msg(sprintf(" -> P-value                  : %s", format_pval(pval_abund)))
  log_msg("")
  log_msg(sprintf("[Plot E: %s]", title_therm))
  log_msg(sprintf(" -> Note                     : %s", note_impact))
  log_msg(sprintf(" -> GAM P-value              : %s", format_pval(pval_direct)))
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 5. Plotting Climate Impact (Individual & Combined)
  # ------------------------------------------------------------------- #
  subtitle_supp <- sprintf("Decline of Phase-Aligned Abundance (Slope = %.3f, %s)", slope_abund, format_pval(pval_abund))
  
  p9_abund <- ggplot(analysis_df, aes(x = Year, y = Log_Abund)) +
    geom_boxplot(aes(group=Year_Fct), fill=neg_fill, alpha=0.3, color=neg_base_color, outlier.shape=NA) +
    geom_jitter(color=neg_base_color, alpha=0.6, width=0.15) + 
    geom_smooth(method="lm", color=trend_color, fill=trend_color, alpha=0.2, linewidth=1.2) +
    labs(title=title_supp, subtitle=subtitle_supp, x="Year", y="Log10(Lag-Adjusted POW Abundance + 1)") +
    scale_x_continuous(breaks = min(analysis_df$Year):max(analysis_df$Year)) + 
    theme(plot.title = element_text(face="bold"))
  
  subtitle_therm <- paste("GAM correlation between Trigger Temp & Abundance (", format_pval(pval_direct), ")")
  
  p9_direct <- ggplot(analysis_df, aes(x = Temperature, y = Log_Abund)) +
    geom_point(color=neg_base_color, size=2.5, alpha=0.6) + 
    geom_smooth(method="gam", formula=y~s(x, bs="cs", k=6), color=trend_color, fill=trend_color, alpha=0.2, linetype="solid", linewidth=1.2) +
    labs(title=title_therm, subtitle=subtitle_therm, x="POW Trigger Temp (°C)", y="Log10(Lag-Adjusted POW Abundance + 1)") +
    theme(plot.title = element_text(face="bold"))
  
  p9_combined <- plot_grid(p9_abund, p9_direct, ncol=2, align="h")
  
  # ------------------------------------------------------------------- #
  # Section 6. Export (All Variations)
  # ------------------------------------------------------------------- #
  if(enable_save_outputs) { 
    # 1. 5개의 개별 플롯 저장
    ggsave(f_plot_env,   plot = p7_temp,   device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_pheno, plot = p7_abund,  device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_anova, plot = p7_anova,  device = "tiff", dpi = 600, width = 8, height = 8, compression = "lzw")
    ggsave(f_plot_supp,  plot = p9_abund,  device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(f_plot_therm, plot = p9_direct, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    
    # 2. 통합 묶음 플롯 저장
    ggsave(f_plot_comb_pw, plot = p7_combined, device = "tiff", dpi = 600, width = 16, height = 8, compression = "lzw") 
    ggsave(f_plot_comb_im, plot = p9_combined, device = "tiff", dpi = 600, width = 13, height = 6, compression = "lzw") 
    
    # 3. CSV 저장
    write.csv(analysis_df, f_csv_impact, row.names = FALSE)
    
    log_msg("[SUCCESS] 5 Individual Plots, 2 Combined Plots, and CSV Data saved successfully.")
  } else {
    log_msg("[SAFE MODE] Plots generated in Viewer (No files saved).")
  }
  
  print(p7_combined)
  print(p9_combined)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

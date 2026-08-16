# [Github: 04_Time-Series_3-7.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 7 (Master Version): Neutral Baseline Diagnostic]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 14"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: 온도에 둔감한 범존종(2_Neutral) 그룹이 기후 변화 속에서도 생태계의 
#       총량을 방어하는 '거시적 안전판(Macro-Buffer)' 역할을 수행함을 증명함.
# 특징:
#   1) [Panel A] 특정 전성기(Peak Window) 부재 증명: 월별(Month) 평탄성 (ANOVA).
#   2) [Panel B] 10년 장기 안정성 증명: 1년 전체(Year-round) 생물량 유지 (Linear Regression).
#   3) [로깅 고도화] 플롯 내 캡션을 제거하고 플롯 Title과 통계치를 1:1 매칭하여 상세 로깅.
#   4) [다중 플롯 분리] 통합 2-Panel 플롯 및 개별 플롯 모두 분리하여 고해상도 저장.
#   5) [시각화 대원칙] Neutral 고유색(Dark Gray) 및 공통 추세선(Green) 적용.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "stats", "agricolae")
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
  library(agricolae) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################

# [1] 공통 GMM 클러스터 설정
# - summary_method: 온도 민감도(IS) 요약 대푯값 기준 (옵션: "median"(권장), "mean")
# - g_clusters    : 군집 분할 수 (권장: 3)
# - target_cluster: 거시적 안정성을 검증할 생태계 안전판(범존종) 군집 (고정 권장: "2_Neutral")
summary_method <- "median"  
g_clusters     <- 3 
target_cluster <- "2_Neutral" 

# [2] 결과 저장 마스터 스위치
# - 옵션: TRUE (폴더에 플롯과 로그 자동 저장), FALSE (뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정 (논문 전반 적용 대원칙)
# ------------------------------------------------------------------- #
neu_base_color <- "#737373" # 2_Neutral 데이터 포인트 테두리 (Dark Gray)
neu_fill       <- "#D9D9D9" # 2_Neutral 박스플롯 내부 배경색 (Light Gray)
trend_color    <- "#347433" # [대원칙] 통계적 추세선 고정 색상 (Green)

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 개편: Part 7)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/07_Neutral_Macro_Stability")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일명 접두사 Part7_ 적용 및 플롯 분리
log_file       <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_Macro_Stability_Log.txt"))

file_plot_A    <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_1_Indiv_Absence_Seasonal_Peak.tiff"))
file_plot_B    <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_2_Indiv_10Yr_Macro_Stability.tiff"))
file_plot_comb <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_Combined_Macro_Stability.tiff"))
file_csv       <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_Macro_Stability_Data.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 7: Macro-Stability of Generalist (Neutral) Group]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Target Cluster   : %s", target_cluster))
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - GMM Clusters     : %d", g_clusters))
  log_msg(sprintf(" - Scope            : Year-round (All months included)"))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included (via tidyr::complete)"))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Common GMM Clustering 
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading Data & Applying True Zero-Included GMM...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
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
  # Section 2. Year-Round Aggregation for Neutral Group
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Aggregating Year-Round Total Abundance for Neutral Group...")
  
  sampled_ym <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Year, Month) %>% dplyr::distinct()
  
  df_cluster_abund <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    dplyr::filter(Cluster == target_cluster) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  # 시차 보정(Time Machine)을 적용하지 않음: 거시적 기초 대사량(Baseline) 자체의 안정성을 보기 위함.
  df_cluster_monthly <- sampled_ym %>%
    dplyr::left_join(df_cluster_abund, by = c("Year", "Month")) %>%
    dplyr::mutate(
      Total_Abund = tidyr::replace_na(Total_Abund, 0),
      Log_Abundance = log10(Total_Abund + 1),
      Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
      Year_Fct = factor(Year)
    )
  
  # ------------------------------------------------------------------- #
  # Section 3. Statistical Testing (ANOVA for Seasonality, LM for Long-term trend)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Running Statistical Diagnostics...")
  
  # 1. Seasonality Test (ANOVA)
  aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
  f_val_anova <- summary(aov_model)[[1]][["F value"]][1]
  pval_anova <- summary(aov_model)[[1]][["Pr(>F)"]][1]
  
  tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
  top_group <- tukey_res$groups
  
  log_msg("\n[Tukey's HSD Letter Grouping (Expectation: Wide spread of 'a' group)]")
  for(i in 1:nrow(top_group)) {
    log_msg(sprintf(" -> %s: Group '%s' (Mean Log_Abund: %.2f)", rownames(top_group)[i], top_group$groups[i], top_group$Log_Abundance[i]))
  }
  
  # 2. Long-term Stability Test (Linear Regression)
  lm_model <- lm(Log_Abundance ~ Year, data = df_cluster_monthly)
  lm_summary <- summary(lm_model)
  slope_abund <- lm_summary$coefficients[2, 1]
  pval_lm <- lm_summary$coefficients[2, 4]
  
  format_pval <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  
  # --- Logging Title, Note, and Stats 1:1 Mapping ---
  title_a <- paste("A. Absence of Seasonal Peak (", target_cluster, ")", sep="")
  note_a  <- "Evaluates the presence of a specific seasonal peak. A flat distribution indicates year-round continuous presence."
  
  title_b <- paste("B. 10-Year Macro-Stability (", target_cluster, ")", sep="")
  note_b  <- "Analyzed using all year-round data without seasonal restriction. A non-significant slope implies macro-level stability."
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg(sprintf("[Plot A: %s]", title_a))
  log_msg(sprintf(" -> Note                     : %s", note_a))
  log_msg(sprintf(" -> Two-way ANOVA F-value    : %.2f", f_val_anova))
  log_msg(sprintf(" -> P-value                  : %s", format_pval(pval_anova)))
  log_msg("")
  log_msg(sprintf("[Plot B: %s]", title_b))
  log_msg(sprintf(" -> Note                     : %s", note_b))
  log_msg(sprintf(" -> Linear Regression Slope  : %.4f", slope_abund))
  log_msg(sprintf(" -> P-value                  : %s", format_pval(pval_lm)))
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 4. Visualization (Individual & Combined Panels)
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Generating Visualization...")
  
  # Panel A: Seasonality Boxplot
  subtitle_a <- paste("Monthly abundance distribution (Two-way ANOVA:", format_pval(pval_anova), ")")
  
  p_season <- ggplot(df_cluster_monthly, aes(x = Month_Fct, y = Log_Abundance)) +
    geom_boxplot(fill = neu_fill, alpha = 0.6, color = neu_base_color, outlier.shape = NA) + 
    geom_jitter(color = neu_base_color, width = 0.15, alpha = 0.6, size = 2.5) +
    labs(title = title_a, subtitle = subtitle_a, x = "Month of the Year", y = "Log10(Cluster Abundance + 1)") +
    theme(plot.title = element_text(face = "bold", size = 14)) 
  
  # Panel B: Long-term Stability Trend
  subtitle_b <- sprintf("Year-round baseline abundance trend (Slope = %.3f, %s)", slope_abund, format_pval(pval_lm))
  
  p_trend <- ggplot(df_cluster_monthly, aes(x = Year, y = Log_Abundance)) +
    geom_boxplot(aes(group=Year_Fct), fill=neu_fill, alpha=0.3, color=neu_base_color, outlier.shape=NA) +
    geom_jitter(color=neu_base_color, alpha=0.6, width=0.15, size=2.5) + 
    geom_smooth(method="lm", color=trend_color, fill=trend_color, alpha=0.2, linewidth=1.5) +
    labs(title = title_b, subtitle = subtitle_b, x = "Year", y = "Log10(Cluster Abundance + 1)") +
    scale_x_continuous(breaks = min(df_cluster_monthly$Year):max(df_cluster_monthly$Year)) + 
    theme(plot.title = element_text(face = "bold", size = 14))
  
  # Combine
  p_combined <- plot_grid(p_season, p_trend, ncol = 2, align = "h")
  
  # ------------------------------------------------------------------- #
  # Section 5. Save & Export
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    # 1. 개별 플롯 저장
    ggsave(file_plot_A, plot = p_season, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    ggsave(file_plot_B, plot = p_trend, device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
    
    # 2. 통합 묶음 플롯 저장
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 14, height = 6, compression = "lzw")
    
    # 3. CSV 데이터 저장
    write.csv(df_cluster_monthly, file_csv, row.names = FALSE)
    log_msg(paste("[SUCCESS] 2 Individual Plots, 1 Combined Plot, and CSV Data saved successfully."))
  } else {
    log_msg("[SAFE MODE] Plots generated in Viewer only (No files saved).")
  }
  
  print(p_combined)
  log_msg("Macro-Stability analysis complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

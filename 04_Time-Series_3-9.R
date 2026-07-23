# [Github: 04_Time-Series_3-9.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 9 (Master Version): Ecosystem Absolute Abundance Balance]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 16"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: 붕괴한 호냉성(Negative) 군집의 '손실된 절대 풍부도(Lost Abundance)'와
#       기회주의적 승리자(Winning ASVs)들의 '획득한 절대 풍부도(Gained Abundance)'를
#       비교하여, Niche 공백의 대체율(Compensation Rate)을 산출함.
# 특징:
#   1) [엄밀성 확보] Biomass, Capacity 등의 추론적 용어를 배제하고 Absolute Abundance로 통일.
#   2) [스케일 전환] Log 스케일의 통계적 필터링 후, 산술 계산은 Linear 스케일로 전환.
#   3) [보상률 계산] (Winners의 10년 총 증가량 / Negative의 10년 총 감소량) * 100 (%)
#   4) [로깅 고도화] 분석 파라미터, 플롯 Note, 생태학적 최종 결론을 로그 파일에 상세 기록.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)
options(scipen = 999) 

# ------------------------------------------------------------------- #
# Section 0. Environment Setup 
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "stats", "agricolae")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot); library(mclust); library(stats); library(agricolae) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################

# [1] 공통 GMM 클러스터 설정
# - summary_method: 온도 민감도(IS) 요약 대푯값 기준 (옵션: "median"(권장), "mean")
# - g_clusters    : 군집 분할 수 (권장: 3)
summary_method <- "median"  
g_clusters     <- 3 

# [2] 붕괴 및 팽창 추적 대상 군집 설정
# - pow_target_cluster: 붕괴를 겪는 최전성기(POW) 타겟 군집 (고정 권장: "1_Negative")
# - candidate_clusters: 빈자리를 차지할 후보 군집 (고정 권장: c("2_Neutral", "3_Positive"))
pow_target_cluster <- "1_Negative"
candidate_clusters <- c("2_Neutral", "3_Positive")

# [3] 통계적 유의성 기준
# - alpha_threshold: 승리자(Winning ASVs) 판별을 위한 선형 회귀 p-value 커트라인 (권장: 0.05)
alpha_threshold <- 0.05

# [4] 결과 저장 마스터 스위치
# - 옵션: TRUE (폴더에 플롯과 로그 자동 저장), FALSE (뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 로깅 설정 (넘버링 개편: Part 9)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
out_dir    <- file.path(base_dir, "05_Phase5_Output/09_Abundance_Balance")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

# 파일명 접두사 Part9_ 적용
log_file  <- file.path(out_dir, paste0("Part9_", toupper(summary_method), "_Abundance_Balance_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part9_", toupper(summary_method), "_Abundance_Balance_Plot.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 9: Absolute Abundance Balance (The Epilogue)]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Collapse Target Cluster : %s", pow_target_cluster))
  log_msg(sprintf(" - Winner Candidates       : %s", paste(candidate_clusters, collapse=", ")))
  log_msg(sprintf(" - Alpha Threshold (p-value): %.2f", alpha_threshold))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1 & 2: Data Load, GMM & POW Extraction
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading Data, applying GMM & Extracting Target Window...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE) %>% dplyr::rename(ASV_Hash = 1)
  
  asv_dict <- df_ab_raw %>% dplyr::select(ASV_Hash, ASV = ASV_ID) %>% dplyr::distinct()
  df_smap <- df_smap_raw %>% dplyr::inner_join(asv_dict, by = "ASV_Hash") %>% dplyr::select(ASV, Best_TP)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  asv_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(agg_IS = median(Temp_IS, na.rm = TRUE), asv_mean = mean(Absolute_Abundance, na.rm=TRUE), .groups = "drop") %>%
    dplyr::filter(asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = order(gmm_model$parameters$mean), labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)
  
  df_monthly <- df_ab_raw %>%
    dplyr::mutate(Year = as.numeric(format(as.Date(Sample_Date), "%Y")), Month = as.numeric(format(as.Date(Sample_Date), "%m"))) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV"))
  
  aov_data <- df_monthly %>% dplyr::filter(Cluster == pow_target_cluster) %>%
    dplyr::group_by(Year, Month) %>% dplyr::summarise(Tot = sum(Absolute_Abundance, na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(Log_Tot = log10(Tot + 1), Month_Fct = factor(Month, levels=1:12, labels=month.abb), Year_Fct = factor(Year))
  tukey_res <- HSD.test(aov(Log_Tot ~ Month_Fct + Year_Fct, data = aov_data), "Month_Fct", group = TRUE)
  collapse_months <- which(month.abb %in% rownames(tukey_res$groups)[grepl("a", tukey_res$groups$groups)]) %>% sort()
  
  # ------------------------------------------------------------------- #
  # Section 3: Phase Alignment & Identifying Winners
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Phase-Aligning Abundances and Identifying Winning ASVs...")
  df_aligned <- df_ab %>%
    tidyr::complete(ASV, Date = unique(df_ab$Date), fill = list(Absolute_Abundance = 0)) %>%
    dplyr::arrange(ASV, as.Date(Date)) %>% 
    dplyr::inner_join(asv_clusters, by = "ASV") %>% 
    dplyr::inner_join(df_smap, by = "ASV") %>%
    dplyr::group_by(ASV) %>%
    dplyr::mutate(shift_n = as.integer(abs(Best_TP[1])), Shifted_Abundance = dplyr::lead(Absolute_Abundance, n = shift_n[1])) %>%
    dplyr::ungroup() %>% dplyr::filter(!is.na(Shifted_Abundance)) %>%
    dplyr::mutate(Year = as.numeric(format(as.Date(Date), "%Y")), Month = as.numeric(format(as.Date(Date), "%m")))
  
  winner_stats <- df_aligned %>%
    dplyr::filter(Month %in% collapse_months & Cluster %in% candidate_clusters) %>%
    dplyr::group_by(ASV, Cluster, Year) %>%
    dplyr::summarise(Yearly_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Yearly_Abund + 1)) %>%
    dplyr::group_by(ASV, Cluster) %>%
    dplyr::summarise(var_abund = var(Log_Abund),
                     Slope = if(var_abund>0) coef(lm(Log_Abund ~ Year))[2] else 0,
                     P_value = if(var_abund>0) summary(lm(Log_Abund ~ Year))$coefficients[2,4] else 1, .groups = "drop") %>%
    dplyr::filter(P_value < alpha_threshold & Slope > 0)
  
  winner_asvs <- winner_stats$ASV
  
  # ------------------------------------------------------------------- #
  # Section 4: Absolute Abundance Balance Calculation (Linear Scale)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Calculating Absolute Linear Trends for Balance Sheet...")
  neg_yearly <- df_aligned %>%
    dplyr::filter(Month %in% collapse_months & Cluster == "1_Negative") %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(Total_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop")
  
  lm_neg <- lm(Total_Abund ~ Year, data = neg_yearly)
  start_year <- min(neg_yearly$Year)
  end_year <- max(neg_yearly$Year)
  neg_start_val <- predict(lm_neg, newdata = data.frame(Year = start_year))
  neg_end_val <- predict(lm_neg, newdata = data.frame(Year = end_year))
  total_lost_abundance <- abs(neg_start_val - neg_end_val) 
  
  if (length(winner_asvs) > 0) {
    win_yearly <- df_aligned %>%
      dplyr::filter(Month %in% collapse_months & ASV %in% winner_asvs) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(Total_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop")
    
    lm_win <- lm(Total_Abund ~ Year, data = win_yearly)
    win_start_val <- predict(lm_win, newdata = data.frame(Year = start_year))
    win_end_val <- predict(lm_win, newdata = data.frame(Year = end_year))
    total_gained_abundance <- abs(win_end_val - win_start_val)
  } else {
    total_gained_abundance <- 0
  }
  
  compensation_rate <- (total_gained_abundance / total_lost_abundance) * 100
  uncompensated_deficit <- total_lost_abundance - total_gained_abundance
  
  # --- Logging Stats and Conclusion ---
  log_msg("\n-------------------------------------------------------------------")
  log_msg("[Ecosystem Abundance Balance Report (Linear Scale)]")
  log_msg(sprintf(" -> Evaluated Window          : Months %s (%d - %d)", paste(collapse_months, collapse="-"), start_year, end_year))
  log_msg(sprintf(" -> Total Lost by 'Negative'  : %12.2f (100.0%%)", total_lost_abundance))
  log_msg(sprintf(" -> Total Gained by 'Winners' : %12.2f ( %5.1f%%)", total_gained_abundance, compensation_rate))
  
  if(uncompensated_deficit > 0) {
    log_msg(sprintf(" -> Uncompensated Deficit     : %12.2f ( %5.1f%%)", uncompensated_deficit, (uncompensated_deficit/total_lost_abundance)*100))
    log_msg("\n [Conclusion] ECOLOGICAL DEGRADATION:")
    log_msg(" Winners failed to fully compensate for the collapse of psychrophilic species.")
    log_msg(" The ecosystem's carrying capacity within the targeted window has been permanently reduced.")
  } else {
    log_msg(sprintf(" -> Abundance Surplus         : %12.2f ( +%.1f%%)", abs(uncompensated_deficit), abs((uncompensated_deficit/total_lost_abundance)*100)))
    log_msg("\n [Conclusion] ECOLOGICAL RESILIENCE & SHIFT:")
    log_msg(" Winners completely filled the niche void left by psychrophilic species.")
    log_msg(" Total abundance was maintained, proving high functional redundancy despite species turnover.")
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 5: Visualization (Balance Bar Chart)
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Generating Mass Balance Visualization...")
  
  title_plot <- "Absolute Abundance Balance (Linear Scale)"
  note_plot  <- "Assesses the absolute linear change in abundance to determine if the expansion of winning ASVs compensated for the loss of the psychrophilic cluster."
  
  log_msg(sprintf("[Plot: %s]", title_plot))
  log_msg(sprintf(" -> Note: %s\n", note_plot))
  
  plot_data <- data.frame(
    Category = factor(c("1. Lost Abundance\n(Negative Cluster)", "2. Compensated Abundance\n(Winning ASVs)", "3. Uncompensated Deficit\n(Niche Vacancy)"),
                      levels = c("1. Lost Abundance\n(Negative Cluster)", "2. Compensated Abundance\n(Winning ASVs)", "3. Uncompensated Deficit\n(Niche Vacancy)")),
    Value = c(total_lost_abundance, total_gained_abundance, max(0, uncompensated_deficit)),
    FillColor = c("#1B4A7E", "#DC2525", "#808080") 
  )
  
  if(uncompensated_deficit < 0) {
    plot_data$Category[3] <- "3. Abundance Surplus\n(Over-compensation)"
    plot_data$Value[3] <- abs(uncompensated_deficit)
    plot_data$FillColor[3] <- "#347433" 
  }
  
  p_balance <- ggplot(plot_data, aes(x = Category, y = Value, fill = FillColor)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", (Value/total_lost_abundance)*100)), vjust = -0.8, size = 5, fontface = "bold") +
    scale_fill_identity() +
    labs(
      title = title_plot,
      subtitle = sprintf("Assessing Niche Compensation Rate During Target Window (Compensation = %.1f%%)", compensation_rate),
      x = "", y = "Absolute Abundance Change (Predicted Linear Trend)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y = element_text(size = 11)
    )
  
  if (enable_save_outputs) {
    # 단일 플롯이므로 하나의 파일로만 저장
    ggsave(file_plot, plot = p_balance, device = "tiff", dpi = 600, width = 10, height = 7, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot generated in Viewer only (No files saved).")
  }
  
  print(p_balance)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

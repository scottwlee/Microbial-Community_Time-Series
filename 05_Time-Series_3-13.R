# [Github: 05_Time-Series_3-13.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 13 (Master Version): Anatomy of Suppression - Cold-favored Losers]
# 
# 목적: 10년간 저온 선호 군집의 붕괴(Suppression)를 개별 ASV 단위로 쪼개어, 
#       실제 어떤 종들이 통계적으로 유의미하게 멸종/감소하고 있는지 규명함.
# 특징:
#   1) [독립적 POW 자동 추출] 저온 선호 군집의 최전성기(POW)를 동적으로 재추출.
#   2) [개별 회귀 분석] POW 내 각 ASV의 시차 보정된 연간 붕괴율(Slope, p-value) 계산.
#   3) [Volcano & Trajectory] 거시적 붕괴 분포(Volcano)와 상위 희생종의 궤적(Facet) 동시 시각화.
#   4) [마스터 로깅 시스템] 통계적으로 명확히 감소한 상위 ASV들의 상세 통계치
#      (Slope, R², APA 형식 p-value, Taxonomy)를 플롯 패널과 1:1 매칭하여 로그 기록.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "stats", "agricolae", "stringr", "broom")
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
  library(stringr)
  library(broom)
})

#################################################
# USER SETTINGS (스위치 및 중요 파라미터 제어)
#################################################

# [1] 결과 저장 마스터 스위치 
enable_save_outputs <- TRUE

# [2] 시각화 대상 및 레이아웃 설정
top_n_losers      <- 6         # 시각화할 가장 가파르게 붕괴한 ASV 개수 (권장: 6)
plot_columns      <- 2         # 통합 플롯(Facet Grid)의 열(Column) 개수
target_taxa_level <- "Class"   # Taxonomy 표시 수준 (Domain ~ Species)

# [3] GMM 및 통계 파라미터 (Part 12와 구조적 일관성 유지)
summary_method       <- "median"  
g_clusters           <- 3 
pow_target_cluster   <- "1_Negative" # 붕괴를 분석할 대상 (저온 선호)
alpha_threshold      <- 0.05         # 유의수준 (통계적으로 명확한 감소 기준)

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정
# ------------------------------------------------------------------- #
neg_base_color <- "#99C2FF" # 일반적인 Cold-favored 붕괴 궤적 (Light Blue)
neg_sig_color  <- "#0065F8" # 유의미하게 붕괴하는 핵심 희생종 (Deep Blue)
trend_color    <- "#347B34" # 추세선 (Green)
trend_fill     <- "#347B34" # 신뢰구간 음영

# ------------------------------------------------------------------- #
# 경로 및 독립적 통합 로깅(Logging) 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
# [수정됨] 04_Phase4_Output -> 04_Phase4_V2_Output
input_dir  <- file.path(base_dir, "04_Phase4_V2_Output/01_Data_Integration") 
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
out_dir    <- file.path(base_dir, "05_Phase5_Output/13_Anatomy_Of_Suppression")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

file_plot_volcano <- file.path(out_dir, paste0("Part13_", toupper(summary_method), "_Volcano_Plot.tiff"))
file_plot_facet   <- file.path(out_dir, paste0("Part13_", toupper(summary_method), "_Top_", top_n_losers, "_Losers_Facet.tiff"))
file_plot_comb    <- file.path(out_dir, paste0("Part13_", toupper(summary_method), "_Combined_Anatomy.tiff"))
file_csv_stats    <- file.path(out_dir, paste0("Part13_", toupper(summary_method), "_Decline_Stats.csv"))

log_file <- file.path(out_dir, paste0("Part13_Anatomy_Log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 13: Anatomy of Suppression (Cold-favored ASVs)]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Target Cluster to Analyze : %s", pow_target_cluster))
  log_msg(sprintf(" - Target Top Losers to Plot : %d", top_n_losers))
  log_msg(sprintf(" - Display Taxonomy Level    : %s", target_taxa_level))
  log_msg(sprintf(" - Alpha (p-value) Threshold : %.2f", alpha_threshold))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Load Data & Dynamic POW Extraction
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading Data and Dynamically Extracting Peak Window (POW)...")
  
  df_ab_raw   <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw   <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE)
  
  colnames(df_smap_raw)[1] <- "ASV_Hash"
  asv_dict <- df_ab_raw %>% dplyr::select(ASV_Hash, ASV_ID) %>% dplyr::distinct()
  df_smap  <- df_smap_raw %>% dplyr::inner_join(asv_dict, by = "ASV_Hash") %>% dplyr::select(ASV = ASV_ID, Best_TP)
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
  
  sampled_ym <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Year, Month) %>% dplyr::distinct()
  
  df_cluster_abund <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    dplyr::filter(Cluster == pow_target_cluster) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  df_cluster_monthly <- sampled_ym %>%
    dplyr::left_join(df_cluster_abund, by = c("Year", "Month")) %>%
    dplyr::mutate(Total_Abund = tidyr::replace_na(Total_Abund, 0),
                  Log_Abundance = log10(Total_Abund + 1),
                  Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
                  Year_Fct = factor(Year))
  
  aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
  tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
  peak_months_str <- rownames(tukey_res$groups)[grepl("a", tukey_res$groups$groups)]
  pow_months <- which(month.abb %in% peak_months_str) %>% sort()
  
  log_msg(sprintf(" -> Dynamically Extracted POW Months: %s", paste(pow_months, collapse="-")))
  
  # ------------------------------------------------------------------- #
  # Section 2. Lag-Alignment and ASV-Level Linear Regression
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Phase-aligning Abundances and Running ASV-level Regressions...")
  
  df_aligned <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::arrange(ASV, as.Date(Date)) %>% 
    dplyr::inner_join(asv_clusters, by = "ASV") %>% 
    dplyr::inner_join(df_smap, by = "ASV") %>%
    dplyr::mutate(Response_Lag_Weeks = abs(Best_TP)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::mutate(shift_n = as.integer(Response_Lag_Weeks[1]), 
                  Shifted_Abundance = dplyr::lead(Absolute_Abundance, n = shift_n[1])) %>%
    dplyr::ungroup() %>% dplyr::filter(!is.na(Shifted_Abundance))
  
  # POW 기간 내의 ASV별 연간 평균 풍부도 계산
  asv_yearly_df <- df_aligned %>%
    dplyr::mutate(Date = as.Date(Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(Month %in% pow_months & Cluster == pow_target_cluster) %>%
    dplyr::group_by(ASV, Year) %>%
    dplyr::summarise(Mean_Shifted_Abund = mean(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Mean_Shifted_Abund + 1))
  
  # 개별 ASV 회귀분석 및 R2 추출
  asv_suppression_stats <- asv_yearly_df %>%
    dplyr::group_by(ASV) %>%
    dplyr::filter(n() >= 5) %>% # 최소 5년 이상의 데이터가 있는 ASV만 대상
    dplyr::group_modify(~ {
      fit <- lm(Log_Abund ~ Year, data = .x)
      tidy_res <- broom::tidy(fit) %>% dplyr::filter(term == "Year")
      r2_val <- summary(fit)$r.squared
      data.frame(slope = tidy_res$estimate, p_value = tidy_res$p.value, R2 = r2_val)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(neg_log_p = -log10(p_value))
  
  # 버블 크기를 위한 10년 전체 평균
  asv_mean_abund <- asv_yearly_df %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(Overall_Mean = mean(Mean_Shifted_Abund, na.rm = TRUE), .groups = "drop")
  
  # ------------------------------------------------------------------- #
  # Section 3. Taxonomy Mapping & Selection of Top Losers (수정됨)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Taxonomy Mapping and Identifying Top Suppressed ASVs...")
  
  taxa_col_mapping <- c("Domain"="L1", "Phylum"="L2", "Class"="L3", "Order"="L4", "Family"="L5", "Genus"="L6", "Species"="L7")
  actual_col_name <- taxa_col_mapping[target_taxa_level]
  
  # 플롯용 Target Taxa와 함께 전체 Taxonomy(L1~L7) 모두 추출
  taxa_extract <- df_ab_raw %>%
    dplyr::select(ASV = ASV_ID, L1, L2, L3, L4, L5, L6, L7) %>%
    dplyr::distinct() %>%
    dplyr::group_by(ASV) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Taxa_Val = !!sym(actual_col_name))
  
  asv_stats_final <- asv_suppression_stats %>%
    dplyr::left_join(asv_mean_abund, by = "ASV") %>%
    dplyr::left_join(taxa_extract, by = "ASV") %>%
    dplyr::mutate(
      Target_Taxa = ifelse(is.na(Taxa_Val) | Taxa_Val == "" | Taxa_Val == "Unassigned", paste0("Unclassified ", target_taxa_level), Taxa_Val),
      Significance = ifelse(p_value < alpha_threshold & slope < 0, "Significant Decline", "Non-significant")
    ) %>%
    dplyr::arrange(slope)
  
  # 가장 붕괴가 심한 ASV 추출
  top_losers <- asv_stats_final %>%
    dplyr::filter(Significance == "Significant Decline") %>%
    head(top_n_losers)
  
  actual_n <- nrow(top_losers)
  
  format_apa_pval <- function(p) {
    if (p < 0.001) return("< .001")
    if (p < 0.01) return("< .01")
    if (p < 0.05) return("< .05")
    return(paste0("= ", sub("^0", "", sprintf("%.3f", p))))
  }
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg("[LOGGED STATS: TOP COLD-FAVORED LOSERS (SEVERE SUPPRESSION)]")
  log_msg(sprintf(" -> Total Cold-favored ASVs analyzed: %d", nrow(asv_stats_final)))
  log_msg(sprintf(" -> Significant Decliners (p < 0.05): %d", sum(asv_stats_final$Significance == "Significant Decline")))
  log_msg("※ Detailed Statistics matching the Trajectory Plots:")
  
  for(i in 1:actual_n) {
    log_msg(sprintf(" [Rank %d] ASV ID: %s | Plot Label (%s): %s", i, top_losers$ASV[i], target_taxa_level, top_losers$Target_Taxa[i]))
    log_msg(sprintf("   - Full Taxa : D:%s | P:%s | C:%s | O:%s | F:%s | G:%s | S:%s", 
                    top_losers$L1[i], top_losers$L2[i], top_losers$L3[i], top_losers$L4[i], top_losers$L5[i], top_losers$L6[i], top_losers$L7[i]))
    log_msg(sprintf("   - Slope     : %.4f", top_losers$slope[i]))
    log_msg(sprintf("   - R-squared : %s", sub("^0", "", sprintf("%.3f", top_losers$R2[i]))))
    log_msg(sprintf("   - p-value   : p %s", format_apa_pval(top_losers$p_value[i])))
    log_msg("   -------------------------------------------------")
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 4. Plotting
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Generating Volcano and Trajectory Plots...")
  
  # [Plot A] Volcano Plot
  p_volcano <- ggplot(asv_stats_final, aes(x = slope, y = neg_log_p)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_hline(yintercept = -log10(alpha_threshold), color = "black", linetype = "dashed", linewidth = 1) +
    geom_point(aes(color = Significance), size = 3, alpha = 0.7) +  # size 맵핑 해제 및 고정값(3) 부여
    scale_color_manual(values = c("Significant Decline" = neg_sig_color, "Non-significant" = neg_base_color)) +
    # scale_size_continuous(range = c(2, 8), guide = "none") +
    labs(title = "A. Anatomy of Cold-favored Suppression",
         x = "Rate of Decline (Slope)", y = "-Log10(p-value)", color = "Status") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0), legend.position = "bottom")
  
  # [Plot B] Facet Trajectories (Top Losers)
  top_losers <- top_losers %>%
    dplyr::mutate(Facet_Label = factor(sprintf("%s\n(%s)", ASV, Target_Taxa), 
                                       levels = sprintf("%s\n(%s)", ASV, Target_Taxa)))
  
  df_top_losers_plot <- asv_yearly_df %>%
    dplyr::inner_join(top_losers %>% dplyr::select(ASV, Facet_Label), by = "ASV")
  
  years_seq <- min(df_top_losers_plot$Year):max(df_top_losers_plot$Year)
  shading_ranges <- data.frame(xmin = years_seq - 0.5, xmax = years_seq + 0.5, ymin = -Inf, ymax = Inf, year = years_seq) %>% dplyr::filter(year %% 2 == 0)
  
  p_trajectories <- ggplot(df_top_losers_plot, aes(x = Year, y = Log_Abund)) +
    geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_point(color = neg_sig_color, size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = trend_color, fill = trend_fill, alpha = 0.15, linetype = "solid", se = TRUE) +
    facet_wrap(~ Facet_Label, scales = "free_y", ncol = plot_columns) +
    scale_x_continuous(breaks = years_seq, expand = expansion(add = 0.5)) +
    labs(title = sprintf("B. Trajectories of Top %d Suppressed Taxa", actual_n),
         x = "Year", y = "Log10(Phase-Aligned Yearly Abund + 1)") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
          strip.background = element_rect(fill = "grey90"), 
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", hjust = 0),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  
  p_combined <- plot_grid(p_volcano, p_trajectories, rel_widths = c(1, 1.3), ncol = 2)
  print(p_combined)
  
  # ------------------------------------------------------------------- #
  # Section 5. Export 
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    ggsave(file_plot_volcano, plot = p_volcano, device = "tiff", dpi = 600, width = 7, height = 6, compression = "lzw")
    ggsave(file_plot_facet, plot = p_trajectories, device = "tiff", dpi = 600, width = 8, height = 8, compression = "lzw")
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 14, height = 7, compression = "lzw")
    write.csv(asv_stats_final, file_csv_stats, row.names = FALSE)
    log_msg("[SUCCESS] All files and logs have been successfully saved to /13_Anatomy_Of_Suppression.")
  } else {
    log_msg("[SAFE MODE] Output generated in Viewer (No files saved).")
  }
  
  log_msg("=== Phase 5 - Part 13 Analysis Complete ===\n")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

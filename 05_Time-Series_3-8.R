# [Github: 04_Time-Series_3-8.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 8 (Master Version): Micro-Niche Replacement Tracking]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 15"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: 호냉성(Negative) 군집이 붕괴하는 최전성기(POW) 기간 동안, 빈 생태적 
#       지위를 차지하며 팽창한 '기회주의적 승리자(Winning ASVs)'를 색출함.
# 특징:
#   1) [시각화 최적화] 롤리팝 플롯(Panel B)의 Y축을 ASV ID로 제한하여 가독성 확보.
#   2) [로깅 고도화] 전체 Taxonomy를 ';' 기호로 병합하여 통계치, 플롯 Note와 함께 로그에 상세 기록.
#   3) [다중 플롯 분리] 통합 2-Panel 플롯 및 개별 플롯(Volcano, Lollipop) 모두 분리 저장.
#   4) [동적 연동] Part 4의 통계 로직을 내장하여 붕괴 기간(POW)을 자동으로 추출해 적용.
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
summary_method <- "median"  
g_clusters     <- 3 

# [2] 붕괴 및 팽창 추적 대상 군집 설정
# - pow_target_cluster: 붕괴를 겪는 최전성기(POW) 타겟 군집 (고정 권장: "1_Negative")
# - peak_selection_method: 붕괴 기간 추출을 위한 통계 알고리즘 (옵션: "tukey"(권장), "kmeans")
# - candidate_clusters: 빈자리를 차지할 후보 군집 (고정 권장: c("2_Neutral", "3_Positive"))
pow_target_cluster    <- "1_Negative" 
peak_selection_method <- "tukey"
candidate_clusters    <- c("2_Neutral", "3_Positive")

# [3] 통계적 유의성 기준
# - alpha_threshold: 승리자(Winning ASVs) 판별을 위한 선형 회귀 p-value 커트라인 (권장: 0.05)
alpha_threshold <- 0.05

# [4] 결과 저장 마스터 스위치
# - 옵션: TRUE (폴더에 플롯과 로그 자동 저장), FALSE (뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정
# ------------------------------------------------------------------- #
neu_color <- "#737373" # 2_Neutral 고유색
pos_color <- "#DC2525" # 3_Positive 고유색
sig_line_color <- "#000000" # 통계적 유의성 기준선
non_sig_color <- "#E0E0E0"  # 유지/정체 ASV 색상

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 개편: Part 8)
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
out_dir    <- file.path(base_dir, "05_Phase5_Output/08_Micro_Niche_Replacement")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

# 파일명 접두사 Part8_ 적용
log_file       <- file.path(out_dir, paste0("Part8_", toupper(summary_method), "_Volcano_Log.txt"))

file_plot_A    <- file.path(out_dir, paste0("Part8_", toupper(summary_method), "_1_Indiv_Volcano.tiff"))
file_plot_B    <- file.path(out_dir, paste0("Part8_", toupper(summary_method), "_2_Indiv_Lollipop.tiff"))
file_plot_comb <- file.path(out_dir, paste0("Part8_", toupper(summary_method), "_Combined_Replacement.tiff"))
file_csv       <- file.path(out_dir, paste0("Part8_", toupper(summary_method), "_Micro_Replacement_Data.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 8: ASV-Level Micro-Niche Replacement Tracking]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Collapse Window defined by  : %s", pow_target_cluster))
  log_msg(sprintf(" - Replacement Candidates      : %s", paste(candidate_clusters, collapse=", ")))
  log_msg(sprintf(" - Alpha Threshold (p-value)   : %.2f", alpha_threshold))
  log_msg(sprintf(" - Summary Method              : %s", toupper(summary_method)))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load, GMM & Complete Taxonomy Extraction
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading Data, Building Taxonomy Dictionary & Applying GMM...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE)
  
  colnames(df_smap_raw)[1] <- "ASV_Hash"
  asv_dict <- df_ab_raw %>% dplyr::select(ASV_Hash, ASV_ID) %>% dplyr::distinct()
  df_smap <- df_smap_raw %>% dplyr::inner_join(asv_dict, by = "ASV_Hash") %>% dplyr::select(ASV = ASV_ID, Best_TP)
  
  # --- Taxonomy 문자열을 ';' 기호로 연결하여 딕셔너리 생성 ---
  tax_cols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxonomy")
  avail_tax_cols <- intersect(colnames(df_ab_raw), tax_cols)
  
  if(length(avail_tax_cols) > 0) {
    tax_dict <- df_ab_raw %>%
      dplyr::select(ASV = ASV_ID, dplyr::all_of(avail_tax_cols)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(ASV) %>% dplyr::slice(1) %>% dplyr::ungroup()
    
    tax_dict$Full_Taxonomy <- apply(tax_dict[avail_tax_cols], 1, function(x) {
      x[is.na(x) | x == "" | x == "NA"] <- "Unassigned"
      paste(x, collapse = ";")
    })
    tax_dict <- tax_dict %>% dplyr::select(ASV, Full_Taxonomy)
  } else {
    tax_dict <- data.frame(ASV = unique(df_ab_raw$ASV_ID)) %>% dplyr::mutate(Full_Taxonomy = "Taxonomy_Not_Available")
  }
  
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
  # Section 2. Dynamic POW Extraction (Auto-Linkage)
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Dynamically extracting Collapse Window (POW) from 1_Negative...")
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
    dplyr::mutate(
      Total_Abund = tidyr::replace_na(Total_Abund, 0),
      Log_Abundance = log10(Total_Abund + 1),
      Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
      Year_Fct = factor(Year)
    )
  
  aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
  tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
  peak_months_str <- rownames(tukey_res$groups)[grepl("a", tukey_res$groups$groups)]
  collapse_months <- which(month.abb %in% peak_months_str) %>% sort()
  
  log_msg(sprintf(" -> Dynamically Extracted Target Window: Months %s", paste(collapse_months, collapse=", ")))
  
  # ------------------------------------------------------------------- #
  # Section 3. Time Machine Phase-Alignment & ASV-Level Aggregation
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Applying Phase-Alignment (Time-Lag adjustment) for all Candidate ASVs...")
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
  
  asv_yearly_abund <- df_aligned %>%
    dplyr::mutate(Date = as.Date(Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(Month %in% collapse_months & Cluster %in% candidate_clusters) %>%
    dplyr::group_by(ASV, Cluster, Year) %>%
    dplyr::summarise(Yearly_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Yearly_Abund + 1))
  
  # ------------------------------------------------------------------- #
  # Section 4. ASV-Level Linear Regression & Taxonomy Logging
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Running independent 10-year linear regressions and logging Taxonomy...")
  
  format_pval <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p))
  
  asv_stats <- asv_yearly_abund %>%
    dplyr::group_by(ASV, Cluster) %>%
    dplyr::summarise(
      n_years = n(),
      var_abund = var(Log_Abund),
      Slope = if(var_abund > 0) coef(lm(Log_Abund ~ Year))[2] else 0,
      P_value = if(var_abund > 0) summary(lm(Log_Abund ~ Year))$coefficients[2,4] else 1,
      .groups = "drop"
    ) %>%
    dplyr::left_join(tax_dict, by = "ASV") %>% 
    dplyr::mutate(
      NegLog10P = -log10(P_value),
      Significance = case_when(
        P_value < alpha_threshold & Slope > 0 ~ "Winner (Significant Expansion)",
        P_value < alpha_threshold & Slope < 0 ~ "Loser (Significant Decline)",
        TRUE ~ "Stable (No Change)"
      ),
      Plot_Color = case_when(
        Significance == "Winner (Significant Expansion)" & Cluster == "2_Neutral" ~ neu_color,
        Significance == "Winner (Significant Expansion)" & Cluster == "3_Positive" ~ pos_color,
        Significance == "Loser (Significant Decline)" ~ "#666666", 
        TRUE ~ non_sig_color 
      )
    )
  
  # --- Logging Title, Note, and Stats 1:1 Mapping ---
  title_a <- "A. Niche Replacement Tracking (Volcano Plot)"
  note_a  <- "Evaluates the expansion rate (Slope) and significance (-Log10 P-value) of Neutral and Positive ASVs during the Negative cluster's collapse window. True Zero-Included phase-aligned abundances were used."
  
  title_b <- "B. Expansion Rates of Winning ASVs"
  note_b  <- "Lollipop chart displaying the specific opportunistic ASVs that significantly expanded into the vacant ecological niches."
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg(sprintf("[Plot A: %s]", title_a))
  log_msg(sprintf(" -> Note                     : %s", note_a))
  log_msg("")
  log_msg(sprintf("[Plot B: %s]", title_b))
  log_msg(sprintf(" -> Note                     : %s", note_b))
  log_msg("-------------------------------------------------------------------\n")
  
  # --- Logging Detailed Winners Info with ';' Separated Taxonomy ---
  winners_df <- asv_stats %>% 
    dplyr::filter(Significance == "Winner (Significant Expansion)") %>% 
    dplyr::arrange(desc(Slope))
  
  num_winners_neu <- sum(winners_df$Cluster == "2_Neutral")
  num_winners_pos <- sum(winners_df$Cluster == "3_Positive")
  total_neu <- sum(asv_stats$Cluster == "2_Neutral")
  total_pos <- sum(asv_stats$Cluster == "3_Positive")
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg("[Micro-Niche Replacement Summary]")
  log_msg(sprintf(" -> Total Neutral ASVs Evaluated : %d | Winners: %d (%.1f%%)", total_neu, num_winners_neu, (num_winners_neu/total_neu)*100))
  log_msg(sprintf(" -> Total Positive ASVs Evaluated: %d | Winners: %d (%.1f%%)", total_pos, num_winners_pos, (num_winners_pos/total_pos)*100))
  
  log_msg("\n[List of Winning ASVs (Taxonomic Identity & Growth Rate)]")
  if(nrow(winners_df) > 0) {
    for(i in 1:nrow(winners_df)) {
      # 로그에 ASV ID, Slope, P-value, ';'로 분리된 전체 Taxonomy 기록
      log_msg(sprintf(" %2d. %s [%s] | Slope: %.4f | p-value: %s | Taxonomy: %s", 
                      i, winners_df$ASV[i], winners_df$Cluster[i], winners_df$Slope[i], format_pval(winners_df$P_value[i]), winners_df$Full_Taxonomy[i]))
    }
  } else {
    log_msg(" -> No significant winners found. Niche replacement failed completely.")
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 5. Visualization (Individual & Combined Panels)
  # ------------------------------------------------------------------- #
  log_msg("Step 5: Generating Plot Objects (Volcano + Clean Lollipop)...")
  
  custom_colors <- c("2_Neutral" = neu_color, "3_Positive" = pos_color)
  
  subtitle_a <- sprintf("ASV-level expansion during Collapse Window (Months %s)", paste(collapse_months, collapse="-"))
  
  # 화산 플롯 생성
  p_volcano <- ggplot(asv_stats, aes(x = Slope, y = NegLog10P, fill = Cluster, color = Cluster)) +
    geom_hline(yintercept = -log10(alpha_threshold), linetype = "dashed", color = sig_line_color, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.6, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    labs(title = title_a, subtitle = subtitle_a, x = "Rate of Expansion (Linear Regression Slope)", y = "-Log10(P-value)") +
    theme(plot.title = element_text(face = "bold", size = 14), legend.position = "none")
  
  # 롤리팝 플롯 생성 (긴 Taxonomy 제거, 깔끔한 ASV ID 사용)
  if(nrow(winners_df) > 0) {
    p_lollipop <- ggplot(winners_df, aes(x = reorder(ASV, Slope), y = Slope, color = Cluster)) +
      geom_segment(aes(x = reorder(ASV, Slope), xend = reorder(ASV, Slope), y = 0, yend = Slope), color = "gray80", linewidth = 1) +
      geom_point(size = 3.5, alpha = 0.8) +
      scale_color_manual(values = custom_colors) +
      coord_flip() +
      labs(title = title_b, subtitle = "Significant replacers/intruders (p < 0.05, Slope > 0)", x = "Opportunistic Species (ASV ID)", y = "Rate of Expansion (Slope)") +
      theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom", legend.title = element_blank())
  } else {
    p_lollipop <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No Significant Winners Found") +
      theme_void() +
      labs(title = title_b) + theme(plot.title = element_text(face = "bold", size = 14))
  }
  
  # 플롯 병합 시 비율을 동등하게(1:1) 맞추어 공간 낭비 최소화
  p_combined <- plot_grid(p_volcano, p_lollipop, ncol = 2, align = "h", rel_widths = c(1, 1))
  
  # ------------------------------------------------------------------- #
  # Section 6. Save & Export
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    # 1. 개별 플롯 저장
    ggsave(file_plot_A, plot = p_volcano, device = "tiff", dpi = 600, width = 8, height = 7, compression = "lzw")
    ggsave(file_plot_B, plot = p_lollipop, device = "tiff", dpi = 600, width = 8, height = 7, compression = "lzw")
    
    # 2. 통합 묶음 플롯 저장
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 14, height = 7, compression = "lzw")
    
    # 3. CSV 데이터 저장
    write.csv(asv_stats %>% arrange(desc(Slope)), file_csv, row.names = FALSE)
    
    log_msg(paste("[SUCCESS] 2 Individual Plots, 1 Combined Plot, and CSV Data saved successfully."))
  } else {
    log_msg("[SAFE MODE] Plots generated in Viewer only (No files saved).")
  }
  
  print(p_combined)
  log_msg("Micro-Niche Replacement analysis complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

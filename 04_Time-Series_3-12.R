# [Github: 05_Time-Series_5-12.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 12 (Master Version): Top Winning ASVs Time-Series Visualization]
# 
# 목적: Part 8에서 식별된 '기회주의적 승리자(Winning ASVs)' 중 팽창 속도(Slope)가 
#       가장 높은 상위 N개의 개별 궤적을 시계열 선형 플롯(Linear Plot)으로 시각화함.
# 특징:
#   1) [생태적 지위 대체 검증] Cold-favored taxa가 붕괴하는 특정 시기(Collapse Window)
#      내에서의 연도별 팽창률만을 정밀 타격하여 경쟁적 해방 현상을 증명함.
#   2) [Taxonomy 다이렉트 매칭] 원시 데이터의 L1~L7 컬럼 구조와 직관적인 분류군 명칭
#      (Phylum, Class 등)을 동적으로 매핑하여 정확한 분류군 정보를 출력함. (버전 호환성 강화)
#   3) [통계값 표시 스위치] 플롯 패널 우측 하단에 Slope 및 P-value 삽입(ON/OFF) 기능.
#   4) [디자인 원칙] 추세선 초록색 실선(#347B34) 및 2x2 배치 통합/개별 플롯 동시 산출 지원.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "stats", "agricolae", "stringr")
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
})

#################################################
# USER SETTINGS (스위치 및 중요 파라미터 제어)
#################################################

# [1] 결과 저장 마스터 스위치 
# TRUE: 지정된 폴더에 플롯(.tiff)과 로그(.txt) 물리적 저장 / FALSE: 뷰어 출력만 실행
enable_save_outputs <- TRUE

# [2] 플롯 내 통계값 표시 스위치
# TRUE: 우측 하단에 통계(Slope, p-value) 표시 / FALSE: 텍스트 없는 Clean Plot
display_stats_on_plot <- TRUE

# [3] 시각화 대상 및 분류군 표기 설정 (L1~L7 매핑)
top_n_asvs        <- 20        # [수정됨] 상위 20개 ASV로 고정
plot_columns      <- 2         # [수정됨] 통합 플롯(Facet Grid)의 열을 2개로 설정 (2x2)
target_taxa_level <- "Class"   # 명시할 수준 입력 ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# [4] 이전 Part 8 분석 파라미터 유지 (정합성 보장용 - 변경 금지 권장)
summary_method        <- "median"  
g_clusters            <- 3 
pow_target_cluster    <- "1_Negative"
candidate_clusters    <- c("2_Neutral", "3_Positive")
alpha_threshold       <- 0.05

# ------------------------------------------------------------------- #
# 시각화 공통 색상 테마 설정
# ------------------------------------------------------------------- #
neu_color     <- "#737373" # Eurythermal taxa 색상
pos_color     <- "#DC2525" # Warm-favored taxa 색상
trend_color   <- "#347B34" # [합의 반영] 추세선 초록색 실선
trend_fill    <- "#347B34" # [합의 반영] 신뢰구간(CI) 초록색 음영

# ------------------------------------------------------------------- #
# 경로 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")
part8_dir  <- file.path(base_dir, "05_Phase5_Output/08_Micro_Niche_Replacement")
out_dir    <- file.path(base_dir, "05_Phase5_Output/12_Top_Winners_TimeSeries")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")
file_part8_csv <- file.path(part8_dir, paste0("Part8_", toupper(summary_method), "_Micro_Replacement_Data.csv"))

log_file       <- file.path(out_dir, paste0("Part12_TopWinners_Plot_Log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
file_plot_comb <- file.path(out_dir, paste0("Part12_Top_", top_n_asvs, "_Winners_Combined_Facet.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 12: Top Winning ASVs Linear Time-Series Visualization]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Target Top ASVs to Plot   : %d", top_n_asvs))
  log_msg(sprintf(" - Display Taxonomy Level    : %s", target_taxa_level))
  log_msg(sprintf(" - Show Stats on Plot        : %s", display_stats_on_plot))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Load Part 8 Results & Filter Top Winners
  # ------------------------------------------------------------------- #
  if(!file.exists(file_part8_csv)) stop("Part 8 CSV file not found. Please run Part 8 first.")
  df_part8 <- read.csv(file_part8_csv, stringsAsFactors = FALSE)
  
  top_winners <- df_part8 %>%
    dplyr::filter(Significance == "Winner (Significant Expansion)") %>%
    dplyr::arrange(desc(Slope)) %>%
    head(top_n_asvs)
  
  if(nrow(top_winners) == 0) stop("No significant winning ASVs found in Part 8 results.")
  
  actual_n <- nrow(top_winners)
  log_msg(sprintf("Step 1: Successfully extracted Top %d ASVs (Target was %d).", actual_n, top_n_asvs))
  target_asvs <- top_winners$ASV
  
  # ------------------------------------------------------------------- #
  # Section 2. Dynamic Data Prep (Linear: Yearly Aggregation)
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Preparing phase-aligned yearly abundances for precise trendline matching...")
  
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
  
  # POW (Collapse Months) 동적 추출
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
  
  log_msg(sprintf(" -> Dynamically Extracted Target Window: Months %s", paste(collapse_months, collapse="-")))
  
  # Time-lag 적용 및 특정 기간 연도별 합산
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
  
  plot_data <- df_aligned %>%
    dplyr::mutate(Date = as.Date(Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(Month %in% collapse_months & ASV %in% target_asvs) %>%
    dplyr::group_by(ASV, Cluster, Year) %>%
    dplyr::summarise(Yearly_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Yearly_Abund + 1))
  
  # ------------------------------------------------------------------- #
  # Section 3. Formatting Taxonomy Labels & Logging Statistics (L1~L7 동적 매핑)
  # ------------------------------------------------------------------- #
  log_msg("\n-------------------------------------------------------------------")
  log_msg("[LOGGED PLOT NOTES: TOP WINNING ASVS TIME-SERIES]")
  
  base_info <- top_winners %>% dplyr::select(ASV, Slope, P_value)
  
  # 데이터셋의 L1~L7 컬럼과 직관적인 분류군 명칭 매핑
  taxa_col_mapping <- c(
    "Domain"  = "L1",
    "Phylum"  = "L2",
    "Class"   = "L3",
    "Order"   = "L4",
    "Family"  = "L5",
    "Genus"   = "L6",
    "Species" = "L7"
  )
  
  actual_col_name <- taxa_col_mapping[target_taxa_level]
  log_msg(sprintf("-> Target Taxonomy Level Displayed: %s (Mapped to column '%s')", target_taxa_level, actual_col_name))
  
  # Base R 인덱싱을 통해 안전하게 추출
  if(!is.na(actual_col_name) && actual_col_name %in% colnames(df_ab_raw)) {
    taxa_extract <- df_ab_raw[, c("ASV_ID", actual_col_name)]
    colnames(taxa_extract) <- c("ASV", "Taxa_Val")
    
    taxa_extract <- taxa_extract %>%
      dplyr::distinct() %>%
      dplyr::group_by(ASV) %>% dplyr::slice(1) %>% dplyr::ungroup()
    
    base_info <- base_info %>%
      dplyr::left_join(taxa_extract, by="ASV") %>%
      dplyr::mutate(Target_Taxa = ifelse(is.na(Taxa_Val) | Taxa_Val == "" | Taxa_Val == "NA" | Taxa_Val == "Unassigned", 
                                         paste0("Unclassified ", target_taxa_level), 
                                         Taxa_Val)) %>%
      dplyr::select(-Taxa_Val)
  } else {
    log_msg(sprintf(" ! WARNING: Column mapping failed. Defaulting to Unclassified."))
    base_info$Target_Taxa <- paste0("Unclassified ", target_taxa_level)
  }
  
  # 패널 라벨 및 통계 라벨 생성
  base_info <- base_info %>%
    dplyr::mutate(
      Facet_Label = sprintf("%s\n(%s)", ASV, Target_Taxa),
      Stat_Label  = sprintf("Slope: +%.4f\np = %.2e", Slope, P_value) 
    ) %>%
    dplyr::arrange(desc(Slope))
  
  ordered_labels <- base_info$Facet_Label
  
  plot_data <- plot_data %>% 
    dplyr::left_join(base_info %>% dplyr::select(ASV, Facet_Label, Stat_Label), by = "ASV") %>%
    dplyr::mutate(
      Facet_Label = factor(Facet_Label, levels = ordered_labels),
      Eco_Group = ifelse(Cluster == "3_Positive", "Warm-favored taxa", "Eurythermal taxa")
    )
  
  log_msg("\n[Top ASVs Detailed Statistics (Mapped perfectly to Plot Panel Titles)]")
  for(i in 1:nrow(base_info)) {
    grp <- ifelse(top_winners$Cluster[top_winners$ASV == base_info$ASV[i]] == "3_Positive", "Warm-favored", "Eurythermal")
    log_msg(sprintf(" [Panel %d] %s", i, gsub("\n", " ", base_info$Facet_Label[i])))
    log_msg(sprintf("   - Group     : %s taxa", grp))
    log_msg(sprintf("   - Slope     : +%.4f", base_info$Slope[i]))
    log_msg(sprintf("   - P-value   : %.2e", base_info$P_value[i]))
  }
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 4. Generate Plots (Combined & Individual)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Generating Combined Facet Plot and Individual Plots...")
  
  custom_colors <- c("Eurythermal taxa" = neu_color, "Warm-favored taxa" = pos_color)
  years_seq <- min(plot_data$Year):max(plot_data$Year)
  shading_ranges <- data.frame(
    xmin = years_seq - 0.5, xmax = years_seq + 0.5, ymin = -Inf, ymax = Inf, year = years_seq
  ) %>% dplyr::filter(year %% 2 == 0)
  
  # 4-1. 공통 테마 설정
  plot_theme <- theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), 
          axis.text.y = element_text(color = "black"),
          strip.background = element_rect(fill = "grey90"), 
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "bottom")
  
  # 4-2. 통합 플롯 (Facet Grid) 생성 (2x2 배치 적용)
  p_combined <- ggplot(plot_data, aes(x = Year, y = Log_Abund)) +
    geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_point(aes(color = Eco_Group), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = trend_color, fill = trend_fill, alpha = 0.15, linetype = "solid", se = TRUE) +
    facet_wrap(~ Facet_Label, scales = "free_y", ncol = plot_columns) +
    scale_color_manual(values = custom_colors) +
    scale_x_continuous(breaks = years_seq, expand = c(0, 0)) +
    labs(title = sprintf("Linear Expansion Trajectories During Collapse Window - Top %d Taxa", actual_n),
         x = "Year", y = "Log10(Phase-Aligned Yearly Abund + 1)", color = "Ecological Group") +
    plot_theme
  
  if (display_stats_on_plot) {
    p_combined <- p_combined + 
      geom_text(aes(x = Inf, y = -Inf, label = Stat_Label), 
                hjust = 1.05, vjust = -0.3, size = 3.5, color = "black", fontface = "italic", 
                inherit.aes = FALSE, check_overlap = TRUE)
  }
  
  print(p_combined)
  
  if (enable_save_outputs) {
    ggsave(file_plot_comb, plot = p_combined, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    log_msg(paste("[SUCCESS] Combined Facet Plot saved to:", file_plot_comb))
  }
  
  # 4-3. 개별 ASV 플롯 분리 생성 및 저장 루프
  log_msg("Step 4: Extracting and saving individual ASV plots...")
  
  for(i in 1:nrow(base_info)) {
    target_asv_id <- base_info$ASV[i]
    indiv_label <- base_info$Facet_Label[i]
    indiv_data <- plot_data %>% dplyr::filter(ASV == target_asv_id)
    
    p_indiv <- ggplot(indiv_data, aes(x = Year, y = Log_Abund)) +
      geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
      geom_point(aes(color = Eco_Group), size = 4, alpha = 0.8) +
      geom_smooth(method = "lm", color = trend_color, fill = trend_fill, alpha = 0.15, linetype = "solid", se = TRUE) +
      scale_color_manual(values = custom_colors) +
      scale_x_continuous(breaks = years_seq, expand = c(0, 0)) +
      labs(title = indiv_label,
           x = "Year", y = "Log10(Phase-Aligned Yearly Abund + 1)", color = "Ecological Group") +
      plot_theme
    
    if (display_stats_on_plot) {
      p_indiv <- p_indiv + 
        geom_text(aes(x = Inf, y = -Inf, label = Stat_Label), 
                  hjust = 1.05, vjust = -0.3, size = 4, color = "black", fontface = "italic", 
                  inherit.aes = FALSE, check_overlap = TRUE)
    }
    
    if (enable_save_outputs) {
      file_plot_indiv <- file.path(out_dir, sprintf("Part12_Rank%02d_%s_Linear.tiff", i, target_asv_id))
      ggsave(file_plot_indiv, plot = p_indiv, device = "tiff", dpi = 600, width = 6, height = 5, compression = "lzw")
      log_msg(sprintf(" -> Saved Individual Plot: Rank %d (%s)", i, target_asv_id))
    }
  }
  
  log_msg("=== Phase 5 - Part 12 Analysis Complete ===")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

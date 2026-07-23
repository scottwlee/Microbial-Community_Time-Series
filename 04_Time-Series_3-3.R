# [Github: 04_Time-Series_3-3.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 3 (Master Version): Unified Taxonomic Composition]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-2.R" 스크립트의 "Phase 5 - Part 19 (b)"에서 이관 및 재정비된 코드입니다.)
#
# 목적: ASV Count(생물다양성)와 Mean Abundance(생태적 우점도) 기반 분류군 조성을 비교함.
# 특징:
#   1) [True Zero-Included] 특정 종이 관찰되지 않은 샘플링 날짜를 0으로 강제 복원.
#   2) [동적 시각화] 절대값(Absolute)/상대비율(Relative) 전환 기능을 제공.
#   3) [로깅 고도화] 플롯별 Title, Note 및 카이제곱(Chi-square) 통계치를 로그에 기록.
#   4) [다중 플롯 분리] 통합 2단 플롯(Combined) 및 개별 플롯(Individual) 동시 자동 저장.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "RColorBrewer", "stats")
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
  library(RColorBrewer) 
  library(stats) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 분석할 분류군 레벨 
#     - 예: "L2" (Phylum), "L3" (Class), "L6" (Genus)
tax_level <- "L3"          

# [2] 'Others' 묶음 처리 임계치 (%)
#     - 군집 내 조성 비율이 이 값 미만인 분류군은 'Others'로 병합 (가독성 향상)
abundance_threshold <- 3.0 

# [3] 플롯 A (다양성 / Taxonomic Richness) Y축 값 표기 방식
#     - 옵션: "absolute" (절대 ASV 개수), "relative" (비율 %)
plot_A_value <- "absolute"  

# [4] 플롯 B (우점도 / Ecological Dominance) Y축 값 표기 방식
#     - 옵션: "absolute" (절대 생물량 평균), "relative" (비율 %)
plot_B_value <- "absolute"  

# [5] 요약 통계량 및 GMM 클러스터 설정
#     - summary_method: "median" (권장), "mean"
#     - g_clusters: 3 (Negative, Neutral, Positive 3분할 권장)
summary_method <- "median"
g_clusters <- 3 

# [6] 시각화 색상 팔레트 및 라벨 교정 딕셔너리
#     - 과거 분류군 명칭을 최신 명칭으로 업데이트
taxa_rename_dict <- c(
  "Proteobacteria" = "Pseudomonadota",
  "Bacteroidetes"  = "Bacteroidota"
)
taxa_palette <- "Spectral" 

# [7] 결과 저장 마스터 스위치
#     - 옵션: TRUE (폴더에 플롯과 로그 자동 저장), FALSE (뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정
# ------------------------------------------------------------------- #
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/03_Unified_Taxonomic_Composition") # 넘버링 갱신

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

value_str <- paste0("A-", plot_A_value, "_B-", plot_B_value)

# 파일명 접두사 Part3_ 적용
log_file          <- file.path(out_dir, paste0("Part3_Composition_", tax_level, "_", value_str, "_TrueZero_Log.txt"))
file_plot_A       <- file.path(out_dir, paste0("Part3_Indiv_A_", tax_level, "_", plot_A_value, "_TrueZero.tiff"))
file_plot_B       <- file.path(out_dir, paste0("Part3_Indiv_B_", tax_level, "_", plot_B_value, "_TrueZero.tiff"))
file_plot_comb    <- file.path(out_dir, paste0("Part3_Combined_", tax_level, "_", value_str, "_TrueZero.tiff"))
file_csv_div      <- file.path(out_dir, paste0("Part3_Data_Diversity_", tax_level, "_TrueZero.csv"))
file_csv_dom      <- file.path(out_dir, paste0("Part3_Data_Dominance_", tax_level, "_TrueZero.csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Taxonomic Renaming
  # ------------------------------------------------------------------- #
  log_msg("Step 1. Loading data and formatting taxonomies...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_taxa <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Taxon = all_of(tax_level)) %>% distinct()
  if (!is.null(taxa_rename_dict) && length(taxa_rename_dict) > 0) {
    for (old_name in names(taxa_rename_dict)) {
      df_taxa$Taxon <- ifelse(df_taxa$Taxon == old_name, taxa_rename_dict[[old_name]], df_taxa$Taxon)
    }
  }
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # ------------------------------------------------------------------- #
  # Section 2. True Zero-Included Aggregation & GMM Clustering
  # ------------------------------------------------------------------- #
  log_msg("Step 2. Aggregating data (True Zero-Included) and applying GMM...")
  
  # [엔진 1] IS(온도 민감도) 대푯값 추출
  asv_is_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance)) %>%
    group_by(ASV) %>%
    summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
              .groups = "drop")
  
  # [엔진 2] 평균 절대 풍부도 추출 (결측치 0 강제 복원)
  all_sample_dates <- unique(df_ab$Date)
  
  asv_ab_summary <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    group_by(ASV) %>%
    summarise(asv_mean = mean(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  # 두 통계 엔진 병합
  asv_summary <- dplyr::inner_join(asv_is_summary, asv_ab_summary, by = "ASV") %>%
    dplyr::filter(!is.na(agg_IS) & asv_mean > 0)
  
  # GMM 클러스터링 적용
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_order <- order(gmm_model$parameters$mean)
  
  asv_clusters <- asv_summary %>%
    mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive")))
  
  df_master <- asv_clusters %>% inner_join(df_taxa, by = "ASV")
  
  # ------------------------------------------------------------------- #
  # Section 3. Calculate Compositions (Absolute vs Relative)
  # ------------------------------------------------------------------- #
  log_msg("Step 3. Calculating Compositions (Diversity vs. Dominance)...")
  
  # (A) 다양성 (Diversity / Richness) - Count 계산
  tax_comp_count <- df_master %>%
    group_by(Cluster, Taxon) %>%
    summarise(ASV_Count = n(), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(Percentage = (ASV_Count / sum(ASV_Count)) * 100) %>% ungroup() %>%
    mutate(Taxon_Clean = ifelse(Percentage < abundance_threshold, "Others", Taxon)) %>%
    group_by(Cluster, Taxon_Clean) %>%
    summarise(ASV_Count = sum(ASV_Count), Percentage = sum(Percentage), .groups = "drop")
  
  # (B) 우점도 (Dominance) - 평균 절대 풍부도(asv_mean) 합산
  tax_comp_abundance <- df_master %>%
    group_by(Cluster, Taxon) %>%
    summarise(Taxon_Abundance = sum(asv_mean), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(Percentage = (Taxon_Abundance / sum(Taxon_Abundance)) * 100) %>% ungroup() %>%
    mutate(Taxon_Clean = ifelse(Percentage < abundance_threshold, "Others", Taxon)) %>%
    group_by(Cluster, Taxon_Clean) %>%
    summarise(Taxon_Abundance = sum(Taxon_Abundance), Percentage = sum(Percentage), .groups = "drop")
  
  # ------------------------------------------------------------------- #
  # Section 4. Statistical Testing & Advanced Logging
  # ------------------------------------------------------------------- #
  log_msg("Step 4. Performing Statistical Tests and Logging...")
  set.seed(414)
  
  # 다양성 분포 차이 통계 검정 (Chi-square)
  contingency_count <- xtabs(ASV_Count ~ Cluster + Taxon_Clean, data = tax_comp_count)
  chi_count <- chisq.test(contingency_count, simulate.p.value = TRUE, B = 2000)
  pval_count <- ifelse(chi_count$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", chi_count$p.value))
  
  # 동적 라벨링 설정
  x_label_dynamic <- paste0("GMM Cluster (", tools::toTitleCase(summary_method), " IS)")
  
  y_var_A   <- ifelse(plot_A_value == "absolute", "ASV_Count", "Percentage")
  title_A   <- ifelse(plot_A_value == "absolute", "(A) Taxonomic Richness", "(A) Taxonomic Diversity")
  sub_A     <- ifelse(plot_A_value == "absolute", "Absolute number of unique ASVs (Richness)", "Relative proportion of unique ASVs (%)")
  ylab_A    <- ifelse(plot_A_value == "absolute", "Total ASV Count (Richness)", "Relative Proportion (%)")
  
  y_var_B   <- ifelse(plot_B_value == "absolute", "Taxon_Abundance", "Percentage")
  title_B   <- ifelse(plot_B_value == "absolute", "(B) Ecological Dominance (Abundance)", "(B) Ecological Dominance (Relative %)")
  sub_B     <- ifelse(plot_B_value == "absolute", "Absolute mean abundance (True Zero-Included)", "Abundance-weighted relative proportion (%)")
  ylab_B    <- ifelse(plot_B_value == "absolute", "Absolute Mean Abundance", "Relative Proportion (%)")
  
  note_shared <- sprintf("Taxa with < %.1f%% abundance are grouped into 'Others' (True Zero-Included). Colors are globally synchronized.", abundance_threshold)
  stat_A      <- sprintf("Pearson's Chi-squared test: X-squared = %.2f, %s", chi_count$statistic, pval_count)
  
  log_msg("\n===================================================================")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Tax Level          : %s", tax_level))
  log_msg(sprintf(" - Others Threshold   : %.1f%%", abundance_threshold))
  log_msg(sprintf(" - Zero Handling      : True Zero-Included (via tidyr::complete)"))
  log_msg(sprintf(" - Total ASVs Mapped  : %d", nrow(asv_summary)))
  log_msg("-------------------------------------------------------------------")
  log_msg(" [Plot A: Taxonomic Diversity/Richness]")
  log_msg(sprintf("  -> Title: %s", title_A))
  log_msg(sprintf("  -> Note : %s", note_shared))
  log_msg(sprintf("  -> Stats: %s", stat_A))
  log_msg("")
  log_msg(" [Plot B: Ecological Dominance]")
  log_msg(sprintf("  -> Title: %s", title_B))
  log_msg(sprintf("  -> Note : %s", note_shared))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 5. Ordered Legend & Dynamic Color Palette
  # ------------------------------------------------------------------- #
  log_msg("Step 5. Generating selected palette and ordering legends...")
  
  all_taxa <- unique(c(tax_comp_count$Taxon_Clean, tax_comp_abundance$Taxon_Clean))
  special_taxa <- c("Unassigned", "Others")
  regular_taxa <- sort(setdiff(all_taxa, special_taxa))
  existing_special <- special_taxa[special_taxa %in% all_taxa]
  universal_levels <- c(regular_taxa, existing_special)
  
  tax_comp_count$Taxon_Clean <- factor(tax_comp_count$Taxon_Clean, levels = universal_levels)
  tax_comp_abundance$Taxon_Clean <- factor(tax_comp_abundance$Taxon_Clean, levels = universal_levels)
  
  color_count <- length(regular_taxa)
  max_pal_colors <- brewer.pal.info[taxa_palette, "maxcolors"]
  
  if (color_count <= max_pal_colors) {
    regular_colors <- brewer.pal(max_pal_colors, taxa_palette)[1:color_count]
  } else {
    get_palette <- colorRampPalette(brewer.pal(max_pal_colors, taxa_palette))
    regular_colors <- get_palette(color_count)
  }
  names(regular_colors) <- regular_taxa
  special_colors <- c("Unassigned" = "#A6A6A6", "Others" = "#E0E0E0") 
  universal_colors <- c(regular_colors, special_colors[existing_special])
  
  # ------------------------------------------------------------------- #
  # Section 6. Plotting (Individual and Combined Generation)
  # ------------------------------------------------------------------- #
  log_msg("Step 6. Generating Plot objects...")
  
  plot_theme <- theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray20", margin = margin(b = 10)),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",          
    legend.justification = "center",     
    legend.box.just = "center",          
    legend.title = element_blank(),      
    legend.text = element_text(size = 10, margin = margin(r = 15, l = 3)), 
    legend.spacing.x = unit(0.4, "cm"),                                    
    legend.spacing.y = unit(0.3, "cm"),                                    
    legend.key.size = unit(0.5, "cm"),                                     
    legend.margin = margin(t = 15, b = 5)                                  
  )
  
  # (A) 다양성 플롯 생성
  p_A <- ggplot(tax_comp_count, aes(x = Cluster, y = .data[[y_var_A]], fill = Taxon_Clean)) +
    geom_bar(stat = "identity", width = 0.6, color = "gray30", linewidth = 0.3) +
    scale_fill_manual(values = universal_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_A, subtitle = sub_A, x = x_label_dynamic, y = ylab_A) +
    plot_theme + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  
  # (B) 우점도 플롯 생성
  p_B <- ggplot(tax_comp_abundance, aes(x = Cluster, y = .data[[y_var_B]], fill = Taxon_Clean)) +
    geom_bar(stat = "identity", width = 0.6, color = "gray30", linewidth = 0.3) +
    scale_fill_manual(values = universal_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_B, subtitle = sub_B, x = x_label_dynamic, y = ylab_B) +
    plot_theme + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  
  # (C) 2단 결합 플롯 생성 (범례는 하단 공유)
  shared_legend <- get_legend(p_A)
  p_row <- plot_grid(
    p_A + theme(legend.position = "none"), 
    p_B + theme(legend.position = "none"), 
    ncol = 2, align = "h"
  )
  final_plot_combined <- plot_grid(p_row, shared_legend, ncol = 1, rel_heights = c(1, 0.2))
  
  # ------------------------------------------------------------------- #
  # Section 7. Export (All Variations)
  # ------------------------------------------------------------------- #
  if(enable_save_outputs) {
    # 1. 2단 결합 플롯 저장
    ggsave(file_plot_comb, plot = final_plot_combined, device = "tiff", dpi = 600, width = 14, height = 9, compression = "lzw")
    # 2. 개별 플롯 A, B 저장
    ggsave(file_plot_A, plot = p_A, device = "tiff", dpi = 600, width = 8, height = 9, compression = "lzw")
    ggsave(file_plot_B, plot = p_B, device = "tiff", dpi = 600, width = 8, height = 9, compression = "lzw")
    # 3. CSV 데이터 저장
    write.csv(tax_comp_count, file_csv_div, row.names = FALSE)
    write.csv(tax_comp_abundance, file_csv_dom, row.names = FALSE)
    
    log_msg("[SUCCESS] Combined Plot, Individual Plots (A & B), and CSV data saved successfully.")
  } else {
    log_msg("[SAFE MODE] Plots generated in Viewer (No files saved).")
  }
  
  print(final_plot_combined)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

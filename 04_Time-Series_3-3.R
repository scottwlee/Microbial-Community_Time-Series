################################################################################
# [Phase 5 - Part 19: Unified Taxonomic Composition (Dynamic Values & Legend Fixed)]
# 목적: ASV Count(생물다양성)와 Mean Abundance(생태적 우점도) 기반 분류군 조성을 비교함.
# 특징:
#   1) [자동 패키지 설치] 환경에 누락된 필수 라이브러리를 감지하고 자동 설치함.
#   2) [plot_A_value], [plot_B_value] 스위치를 통해 플롯 A와 B를 각각 
#      'absolute'(절대값) 또는 'relative'(비율, 100%)로 자유롭게 전환 가능.
#   3) (B) Absolute 계산 시 'IS vs Mean' 플롯과의 완벽한 수학적 일관성을 위해 
#      개별 ASV의 '10-year Mean Abundance'를 합산하여 군집의 절대 생물량을 도출.
#   4) 고해상도 TIFF 추출 시 하단 범례 겹침(Overlap) 방지 여백 최적화 적용.
#   5) [enable_save_outputs] 원터치 마스터 스위치로 파일/로그/CSV 저장 여부 완벽 제어.
################################################################################

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 1) 필수 패키지 목록 정의 및 자동 설치
required_packages <- c("dplyr", "ggplot2", "cowplot", "mclust", "RColorBrewer", "stats")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

# 2) 패키지 로드
suppressPackageStartupMessages({
  library(dplyr)
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
# [1] 분류군 해상도 (옵션: "L2" (Phylum), "L3" (Class), "L4" (Order))
# - 플롯에 표시할 분류학적 계급 수준.
tax_level <- "L3"          

# [2] 통합 임계값 (단위: %, 상대 비율 기준. 권장: 2.0 ~ 5.0)
# - 이 수치 미만의 비율을 차지하는 미소 분류군은 'Others'로 묶어 가독성 향상.
abundance_threshold <- 3.0 

# [3] 플롯 값 표현 방식 선택 (옵션: "absolute", "relative")
# - "absolute": (A) 실제 ASV 개수 합산(Richness), (B) 평균 절대 풍부도 합산.
# - "relative": (A) ASV 개수 비율(100%), (B) 풍부도 비율(100%).
plot_A_value <- "absolute"  # 권장: absolute (종 풍부도 비교 증명용)
plot_B_value <- "absolute"  # 권장: absolute (IS vs Mean 플롯과 논리적 연계용)

# [4] 플롯 출력 모드 및 캡션 설정
# - plot_mode 옵션: "combined" (하나의 패널로 묶음), "individual" (개별 플롯 생성)
plot_mode <- "combined"     
include_caption <- TRUE     

# [5] 분석 파라미터 일관성 유지 (Part 4, 5, 6과 동일)
# - summary_method (옵션: "median", "mean") : IS 대푯값 산정 기준
# - g_clusters (고정: 3) : 생태적 지위 그룹 개수
# - exclude_zeros (옵션: TRUE, FALSE) : 휴면기 데이터 배제 여부
summary_method <- "median"
g_clusters <- 3 
exclude_zeros <- FALSE

# [6] 특정 분류군 이름 변경 딕셔너리 (옵션)
# - 최신 분류학 명칭으로 업데이트하기 위한 딕셔너리 (불필요시 비워둠).
taxa_rename_dict <- c(
  "Proteobacteria" = "Pseudomonadota",
  "Bacteroidetes"  = "Bacteroidota"
)

# [7] 분류군 시각화 컬러 팔레트 (옵션: "Paired", "Set3", "Spectral" 등)
taxa_palette <- "Spectral" 

# [8] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
# - TRUE일 경우 디렉토리 생성, 로그 기록, TIFF 및 CSV 파일 저장을 실제로 수행함.
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/17_Unified_Taxonomic_Composition")

# [원터치 제어] 디렉토리 생성
if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일명에 세팅 상태 동적 반영
cap_str   <- ifelse(include_caption, "CapON", "CapOFF")
zero_str  <- ifelse(exclude_zeros, "ZeroExc", "ZeroInc")
value_str <- paste0("A-", plot_A_value, "_B-", plot_B_value)

log_file <- file.path(out_dir, paste0("Part19_Composition_", tax_level, "_", value_str, "_", cap_str, "_Log.txt"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  # [원터치 제어] 텍스트 파일에 로그 기록
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # -------------------------------------------------------------------
  # Section 1. Data Load & Taxonomic Renaming
  # -------------------------------------------------------------------
  log_msg("Step 1. Loading data and formatting taxonomies...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_taxa <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Taxon = all_of(tax_level)) %>% distinct()
  if (!is.null(taxa_rename_dict) && length(taxa_rename_dict) > 0) {
    for (old_name in names(taxa_rename_dict)) {
      df_taxa$Taxon <- ifelse(df_taxa$Taxon == old_name, taxa_rename_dict[[old_name]], df_taxa$Taxon)
    }
  }
  
  df_merged <- dplyr::inner_join(
    df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance),
    df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength),
    by = c("ASV", "Date")
  ) %>% dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  if (exclude_zeros) {
    df_merged <- df_merged %>% filter(Absolute_Abundance > 0)
    zero_text <- "Zero-excluded"
  } else {
    zero_text <- "Zero-included"
  }
  
  # -------------------------------------------------------------------
  # Section 2. GMM Clustering & Master Data Frame
  # -------------------------------------------------------------------
  log_msg("Step 2. Aggregating data and applying GMM clustering...")
  
  asv_summary <- df_merged %>%
    group_by(ASV) %>%
    summarise(
      agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE), 
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE)
    ) %>% 
    filter(asv_mean > 0 & !is.na(agg_IS))
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_order <- order(gmm_model$parameters$mean)
  
  asv_clusters <- asv_summary %>%
    mutate(Cluster = factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive")))
  
  df_master <- asv_clusters %>% inner_join(df_taxa, by = "ASV")
  
  # -------------------------------------------------------------------
  # Section 3. Calculate Compositions (Absolute vs Relative)
  # -------------------------------------------------------------------
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
  
  # -------------------------------------------------------------------
  # Section 4. Statistical Testing (Chi-square for Count ONLY)
  # -------------------------------------------------------------------
  log_msg("Step 4. Performing Statistical Tests (Applied to Count data only)...")
  set.seed(414)
  
  contingency_count <- xtabs(ASV_Count ~ Cluster + Taxon_Clean, data = tax_comp_count)
  chi_count <- chisq.test(contingency_count, simulate.p.value = TRUE, B = 2000)
  pval_count <- ifelse(chi_count$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", chi_count$p.value))
  
  # -------------------------------------------------------------------
  # Section 5. Ordered Legend & Dynamic Color Palette
  # -------------------------------------------------------------------
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
  
  # -------------------------------------------------------------------
  # Section 6. Plotting Parameters Setup (Dynamic Labels & Theme)
  # -------------------------------------------------------------------
  log_msg("Step 6. Assembling plots dynamically based on settings...")
  
  x_label_dynamic <- paste0("GMM Cluster (", tools::toTitleCase(summary_method), " IS)")
  
  # (A) 플롯 동적 라벨링
  y_var_A   <- ifelse(plot_A_value == "absolute", "ASV_Count", "Percentage")
  title_A   <- ifelse(plot_A_value == "absolute", "(A) Taxonomic Richness", "(A) Taxonomic Diversity")
  sub_A     <- ifelse(plot_A_value == "absolute", "Absolute number of unique ASVs (Richness)", "Relative proportion of unique ASVs (%)")
  ylab_A    <- ifelse(plot_A_value == "absolute", "Total ASV Count (Richness)", "Relative Proportion (%)")
  
  # (B) 플롯 동적 라벨링 
  y_var_B   <- ifelse(plot_B_value == "absolute", "Taxon_Abundance", "Percentage")
  title_B   <- ifelse(plot_B_value == "absolute", "(B) Ecological Dominance (Abundance)", "(B) Ecological Dominance (Relative %)")
  sub_B     <- ifelse(plot_B_value == "absolute", "Absolute mean abundance (Linked to 'IS vs Mean' plot)", "Abundance-weighted relative proportion (%)")
  ylab_B    <- ifelse(plot_B_value == "absolute", "Absolute Mean Abundance", "Relative Proportion (%)")
  
  cap_line1 <- sprintf("Note: Taxa with < %.1f%% abundance are grouped into 'Others' (%s). Colors are globally synchronized.", abundance_threshold, zero_text)
  cap_line2 <- sprintf("[Plot A] %s | Pearson's Chi-squared test: X-squared = %.2f, %s", title_A, chi_count$statistic, pval_count)
  cap_line3 <- sprintf("[Plot B] %s | %s", title_B, sub_B)
  full_caption <- paste(cap_line1, cap_line2, cap_line3, sep = "\n")
  
  # 범례 겹침 완벽 방지를 위한 여백 세밀화 테마
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
  
  # -------------------------------------------------------------------
  # Section 7. Render & Export (원터치 제어)
  # -------------------------------------------------------------------
  p_A <- ggplot(tax_comp_count, aes(x = Cluster, y = .data[[y_var_A]], fill = Taxon_Clean)) +
    geom_bar(stat = "identity", width = 0.6, color = "gray30", linewidth = 0.3) +
    scale_fill_manual(values = universal_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_A, subtitle = sub_A, x = x_label_dynamic, y = ylab_A) +
    plot_theme + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  
  p_B <- ggplot(tax_comp_abundance, aes(x = Cluster, y = .data[[y_var_B]], fill = Taxon_Clean)) +
    geom_bar(stat = "identity", width = 0.6, color = "gray30", linewidth = 0.3) +
    scale_fill_manual(values = universal_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_B, subtitle = sub_B, x = x_label_dynamic, y = ylab_B) +
    plot_theme + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  
  if (plot_mode == "combined") {
    shared_legend <- get_legend(p_A)
    p_row <- plot_grid(
      p_A + theme(legend.position = "none"), 
      p_B + theme(legend.position = "none"), 
      ncol = 2, align = "h"
    )
    
    if (include_caption) {
      caption_grob <- ggdraw() + draw_label(full_caption, x = 0.01, y = 0.8, hjust = 0, vjust = 1, fontface = "italic", size = 10, color = "gray30", lineheight = 1.4)
      final_plot <- plot_grid(p_row, shared_legend, caption_grob, ncol = 1, rel_heights = c(1, 0.2, 0.15))
    } else {
      final_plot <- plot_grid(p_row, shared_legend, ncol = 1, rel_heights = c(1, 0.2))
    }
    
    file_plot <- file.path(out_dir, paste0("Part19_Combined_", tax_level, "_", value_str, "_", zero_str, "_", cap_str, ".tiff"))
    file_csv_div <- file.path(out_dir, paste0("Part19_Data_Diversity_", tax_level, "_", zero_str, ".csv"))
    file_csv_dom <- file.path(out_dir, paste0("Part19_Data_Dominance_", tax_level, "_", zero_str, ".csv"))
    
    if(enable_save_outputs) {
      ggsave(file_plot, plot = final_plot, device = "tiff", dpi = 600, width = 14, height = 9, compression = "lzw")
      write.csv(tax_comp_count, file_csv_div, row.names = FALSE)
      write.csv(tax_comp_abundance, file_csv_dom, row.names = FALSE)
      log_msg("[SUCCESS] Combined Plot and CSV data saved successfully.")
    } else {
      log_msg("[SAFE MODE] Combined Plot generated in Viewer (No files saved).")
    }
    print(final_plot)
    
  } else if (plot_mode == "individual") {
    if (include_caption) {
      p_A_final <- p_A + labs(caption = paste(cap_line1, cap_line2, sep="\n")) + theme(plot.caption = element_text(hjust=0, color="gray30", face="italic", margin=margin(t=15)))
      p_B_final <- p_B + labs(caption = paste(cap_line1, cap_line3, sep="\n")) + theme(plot.caption = element_text(hjust=0, color="gray30", face="italic", margin=margin(t=15)))
    } else {
      p_A_final <- p_A
      p_B_final <- p_B
    }
    
    file_plot_A <- file.path(out_dir, paste0("Part19_Indiv_A_", tax_level, "_", plot_A_value, "_", zero_str, "_", cap_str, ".tiff"))
    file_plot_B <- file.path(out_dir, paste0("Part19_Indiv_B_", tax_level, "_", plot_B_value, "_", zero_str, "_", cap_str, ".tiff"))
    file_csv_div <- file.path(out_dir, paste0("Part19_Data_Diversity_", tax_level, "_", zero_str, ".csv"))
    file_csv_dom <- file.path(out_dir, paste0("Part19_Data_Dominance_", tax_level, "_", zero_str, ".csv"))
    
    if(enable_save_outputs) {
      ggsave(file_plot_A, plot = p_A_final, device = "tiff", dpi = 600, width = 8, height = 10, compression = "lzw")
      ggsave(file_plot_B, plot = p_B_final, device = "tiff", dpi = 600, width = 8, height = 10, compression = "lzw")
      write.csv(tax_comp_count, file_csv_div, row.names = FALSE)
      write.csv(tax_comp_abundance, file_csv_dom, row.names = FALSE)
      log_msg("[SUCCESS] Individual Plots and CSV data saved successfully.")
    } else {
      log_msg("[SAFE MODE] Individual Plots generated in Viewer (No files saved).")
    }
    print(p_A_final)
    print(p_B_final)
  }
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

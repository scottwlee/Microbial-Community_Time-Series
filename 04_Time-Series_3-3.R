################################################################################
# [Phase 5 - Part 19 (b) (Master Version): Unified Taxonomic Composition]
# 목적: ASV Count(생물다양성)와 Mean Abundance(생태적 우점도) 기반 분류군 조성을 비교함.
# 특징:
#   1) [핵심 혁신: True Zero-Included] Part 4, 5, 7과 완벽하게 동일한 
#      tidyr::complete 로직을 적용하여 극우점종/희귀종의 생물량 과대평가 오류 차단.
#   2) [plot_A_value], [plot_B_value] 스위치를 통해 절대값/비율 전환 가능.
#   3) 고해상도 TIFF 추출 시 하단 범례 겹침 방지 여백 최적화 적용.
################################################################################

options(stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
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
# USER SETTINGS 
#################################################
tax_level <- "L3"          
abundance_threshold <- 3.0 

plot_A_value <- "absolute"  
plot_B_value <- "absolute"  

plot_mode <- "combined"     
include_caption <- TRUE     

summary_method <- "median"
g_clusters <- 3 
# [수정] exclude_zeros 옵션 삭제 (강제 Zero-Included 적용으로 논문 통일성 확보)

taxa_rename_dict <- c(
  "Proteobacteria" = "Pseudomonadota",
  "Bacteroidetes"  = "Bacteroidota"
)
taxa_palette <- "Spectral" 
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/17_Unified_Taxonomic_Composition")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

cap_str   <- ifelse(include_caption, "CapON", "CapOFF")
value_str <- paste0("A-", plot_A_value, "_B-", plot_B_value)

# 파일명에 TrueZero 명시
log_file <- file.path(out_dir, paste0("Part19_Composition_", tax_level, "_", value_str, "_TrueZero_", cap_str, "_Log.txt"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
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
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # -------------------------------------------------------------------
  # Section 2. True Zero-Included Aggregation & GMM Clustering
  # -------------------------------------------------------------------
  log_msg("Step 2. Aggregating data (True Zero-Included) and applying GMM...")
  
  # [엔진 1] IS(온도 민감도) 대푯값 추출 (상호작용이 존재하는 날 기준)
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
  sub_B     <- ifelse(plot_B_value == "absolute", "Absolute mean abundance (True Zero-Included)", "Abundance-weighted relative proportion (%)")
  ylab_B    <- ifelse(plot_B_value == "absolute", "Absolute Mean Abundance", "Relative Proportion (%)")
  
  cap_line1 <- sprintf("Note: Taxa with < %.1f%% abundance are grouped into 'Others' (True Zero-Included). Colors are globally synchronized.", abundance_threshold)
  cap_line2 <- sprintf("[Plot A] %s | Pearson's Chi-squared test: X-squared = %.2f, %s", title_A, chi_count$statistic, pval_count)
  cap_line3 <- sprintf("[Plot B] %s | %s", title_B, sub_B)
  full_caption <- paste(cap_line1, cap_line2, cap_line3, sep = "\n")
  
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
    
    file_plot <- file.path(out_dir, paste0("Part19_Combined_", tax_level, "_", value_str, "_TrueZero_", cap_str, ".tiff"))
    file_csv_div <- file.path(out_dir, paste0("Part19_Data_Diversity_", tax_level, "_TrueZero.csv"))
    file_csv_dom <- file.path(out_dir, paste0("Part19_Data_Dominance_", tax_level, "_TrueZero.csv"))
    
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
    
    file_plot_A <- file.path(out_dir, paste0("Part19_Indiv_A_", tax_level, "_", plot_A_value, "_TrueZero_", cap_str, ".tiff"))
    file_plot_B <- file.path(out_dir, paste0("Part19_Indiv_B_", tax_level, "_", plot_B_value, "_TrueZero_", cap_str, ".tiff"))
    file_csv_div <- file.path(out_dir, paste0("Part19_Data_Diversity_", tax_level, "_TrueZero.csv"))
    file_csv_dom <- file.path(out_dir, paste0("Part19_Data_Dominance_", tax_level, "_TrueZero.csv"))
    
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

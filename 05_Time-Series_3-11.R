# [Github: 04_Time-Series_3-11.R]==================================================#

# ------------------------------------------------------------------- #
# [Phase 5 - Part 11 (Master Version): Phylogenetic Mapping & Mantel Test]
# (※ 기록: 이 스크립트는 과거 "04_Time-Series_3-1.R" 스크립트의 "Phase 5 - Part 17"에서 이관 및 재정비된 코드입니다.)
# 
# 목적: 1) 극단적으로 뭉친 계통수 문제를 해결하기 위해 균일한 분기도(Cladogram)로 시각화.
#       2) Mantel Test를 통해 온도 반응성(Temp_IS)의 계통적 보존성을 통계적으로 검증함.
# 특징:
#   1) [논리적 일관성] 파이프라인 공통의 True Zero-Included 결측치 복원 로직을 
#      강제 적용하여 GMM 클러스터 배정 ASV 목록을 100% 일치시킴.
#   2) [로깅 고도화] 플롯 캡션을 제거하고, Mantel Test 통계치 및 결론을 로그에 매칭 기록.
#   3) [다중 플롯 분리] 단일 ggtree 객체를 독립적인 고해상도 TIFF 파일로 자동 저장.
# ------------------------------------------------------------------- #

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------- #
# Section 0. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------- #
required_packages <- c("dplyr", "tidyr", "ggplot2", "cowplot", "mclust", "phyloseq", "ape", "ggtree", "treeio", "DECIPHER", "phangorn")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  # Bioconductor 패키지가 포함되어 있으므로 설치 로직 분기 필요 시 주의
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_packages)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(mclust)
  library(phyloseq)
  library(ape)      
  library(ggtree)   
  library(treeio)
  library(DECIPHER) 
  library(phangorn) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] IS 통계량 요약 방식 
#     - 옵션: "median", "mean"
#     - 권장: "median" (이상치 방어 및 파이프라인 통일성 유지)
summary_method  <- "median"  

# [2] GMM 클러스터(군집) 개수 
#     - 권장: 3 (Negative, Neutral, Positive 생태학적 3분할)
g_clusters      <- 3 

# [3] 결과 저장 마스터 스위치 
#     - 옵션: TRUE (지정된 폴더에 플롯과 로그 파일 자동 저장), FALSE (RStudio 뷰어 출력만)
enable_save_outputs <- TRUE

# ------------------------------------------------------------------- #
# 경로 및 동적 파일명 설정 (넘버링 갱신: Part 11)
# ------------------------------------------------------------------- #
base_dir <- "/home/scott/EDM_16SV4_PA"
file_phyloseq <- file.path(base_dir, "02_Phase2_Output", "Filtered_Phyloseq_Updated.rds") 

input_dir <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir   <- file.path(base_dir, "05_Phase5_Output/11_Phylogenetic_Mapping")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

# 파일명 접두사 Part11_ 적용
log_file  <- file.path(out_dir, paste0("Part11_", toupper(summary_method), "_Phylogenetic_Mapping_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part11_", toupper(summary_method), "_Phylogenetic_Niche_Cladogram.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Phase 5 - Part 11: Phylogenetic Mapping & Mantel Test]")
  log_msg(" [Analysis Parameters & Settings]")
  log_msg(sprintf(" - Summary Method   : %s", toupper(summary_method)))
  log_msg(sprintf(" - GMM Clusters     : %d", g_clusters))
  log_msg(sprintf(" - Zero Handling    : True Zero-Included (Synchronized)"))
  log_msg("===================================================================\n")
  
  # ------------------------------------------------------------------- #
  # Section 1. Data Load & Clustering (True Zero-Included Logic)
  # ------------------------------------------------------------------- #
  log_msg("Step 1: Loading data and assigning synchronized GMM clusters...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  # IS 대푯값 추출 (결측 아님)
  asv_is_summary <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE), .groups = "drop")
  
  # 결측치 0 강제 복원 (True Zero)
  all_sample_dates <- unique(df_ab$Date)
  asv_ab_summary <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::summarise(asv_mean = mean(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  # 병합 및 필터링
  asv_summary <- dplyr::inner_join(asv_is_summary, asv_ab_summary, by = "ASV") %>%
    dplyr::filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  asv_clusters <- asv_summary %>%
    dplyr::mutate(Cluster = factor(gmm_model$classification, levels = order(gmm_model$parameters$mean), labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster, agg_IS)
  
  target_asvs <- unique(asv_clusters$ASV)
  
  # ------------------------------------------------------------------- #
  # Section 2. Auto-Mapping & Extract Sequences
  # ------------------------------------------------------------------- #
  log_msg("Step 2: Loading phyloseq object and matching ID sequences...")
  
  if(!file.exists(file_phyloseq)) stop("Error: Phyloseq object file not found.")
  ps_obj <- readRDS(file_phyloseq)
  
  tax_matrix <- as(phyloseq::tax_table(ps_obj), "matrix")
  tax_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  
  hash_to_asv <- NULL
  for (col in colnames(tax_df)) {
    if (any(target_asvs %in% tax_df[[col]])) {
      hash_to_asv <- tax_df[[col]]
      names(hash_to_asv) <- rownames(tax_df)
      break
    }
  }
  
  if(is.null(hash_to_asv)) stop("Error: Target ASVs could not be mapped to Phyloseq taxa.")
  
  target_hashes <- names(hash_to_asv)[hash_to_asv %in% target_asvs]
  target_seqs <- phyloseq::refseq(ps_obj)[target_hashes]
  names(target_seqs) <- hash_to_asv[target_hashes]
  
  # ------------------------------------------------------------------- #
  # Section 3. Build Tree (MSA -> NJ Tree)
  # ------------------------------------------------------------------- #
  log_msg("Step 3: Performing Multiple Sequence Alignment (MSA)...")
  alignment <- DECIPHER::AlignSeqs(target_seqs, anchor = NA, verbose = FALSE)
  
  log_msg("Building Phylogenetic Tree using phangorn (NJ)...")
  phang_align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dna_dist <- phangorn::dist.ml(phang_align)
  treeNJ <- phangorn::NJ(dna_dist)
  
  treeNJ$tip.label <- as.character(treeNJ$tip.label)
  
  # ------------------------------------------------------------------- #
  # Section 4. Statistical Testing: Mantel Test & Logging
  # ------------------------------------------------------------------- #
  log_msg("Step 4: Running Mantel Test for Phylogenetic Conservatism...")
  
  phy_dist_mat <- ape::cophenetic.phylo(treeNJ)
  
  asv_to_is <- setNames(asv_clusters$agg_IS, asv_clusters$ASV)
  is_vec <- asv_to_is[treeNJ$tip.label] 
  trait_dist_mat <- as.matrix(dist(is_vec))
  
  set.seed(414)
  mantel_res <- ape::mantel.test(phy_dist_mat, trait_dist_mat, nperm = 999)
  
  # 통계 결과 정리
  if(mantel_res$p < 0.05) {
    conc_str <- "SIGNIFICANT Phylogenetic Conservatism detected. (Close relatives share highly similar temperature interactions.)"
    sig_text <- sprintf("Mantel Test p = %.3f (Significant Conservatism)", mantel_res$p)
  } else {
    conc_str <- "NO Significant Phylogenetic Conservatism. (Thermal traits are randomly distributed or convergent.)"
    sig_text <- sprintf("Mantel Test p = %.3f (No Conservatism)", mantel_res$p)
  }
  
  # 플롯 메타데이터 (로그 1:1 매칭용)
  title_plot <- "Phylogenetic Conservatism of Thermal Niches"
  note_plot  <- sprintf("Tips represent individual ASVs. Branch lengths are ignored to emphasize topology (Cladogram). True Zero-Included GMM mapping. %s", sig_text)
  
  log_msg("\n-------------------------------------------------------------------")
  log_msg(sprintf("[Plot: %s]", title_plot))
  log_msg(sprintf(" -> Note                     : %s", note_plot))
  log_msg("\n[Mantel Test Statistics]")
  log_msg(sprintf(" -> Z-statistic              : %f", mantel_res$z.stat))
  log_msg(sprintf(" -> P-value                  : %.3f (based on 999 permutations)", mantel_res$p))
  log_msg(sprintf(" -> Conclusion               : %s", conc_str))
  log_msg("-------------------------------------------------------------------\n")
  
  # ------------------------------------------------------------------- #
  # Section 5. FANCY Visualization (Circular Cladogram)
  # ------------------------------------------------------------------- #
  log_msg("Step 5: Generating Fancy Circular Cladogram...")
  
  # 논문 표준 색상 테마 적용
  custom_colors <- c("1_Negative" = "#0065F8", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  
  meta_data <- asv_clusters %>% dplyr::select(ASV, Cluster, agg_IS)
  
  p_tree <- ggtree(treeNJ, layout = "fan", open.angle = 15, branch.length = "none", linewidth = 0.4) %<+% meta_data +
    geom_tippoint(aes(color = Cluster), size = 3.5, alpha = 0.9, stroke = 0.3) +
    scale_color_manual(values = custom_colors, name = "Ecological Cluster") +
    labs(
      title = title_plot,
      subtitle = "Topological Cladogram of 16S rRNA mapped with GMM Clusters"
      # Note(Caption)는 로그로 이관됨.
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  # ------------------------------------------------------------------- #
  # Section 6. Save & Export
  # ------------------------------------------------------------------- #
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_tree, device = "tiff", dpi = 600, width = 11, height = 10, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot generated in Viewer only (No files saved).")
  }
  
  print(p_tree)
  log_msg("Phylogenetic mapping and statistical testing complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

##### END. ######################################################################

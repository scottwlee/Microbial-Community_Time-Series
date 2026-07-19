################################################################################
# [Phase 5 - Part 6: Taxonomic Composition of Niche Strategies]
# 목적: GMM으로 분류된 온도 민감도 3개 그룹(Negative, Neutral, Positive)의 
#       분류군 조성(Taxonomic Composition)을 비교 분석함.
# 특징: 
#   1) [자동 패키지 설치] 환경에 누락된 필수 라이브러리를 감지하고 자동 설치함.
#   2) 각 생태적 지위(Temporal Niche)를 차지하는 미생물들이 특정 분류군(Phylum 등)에 
#      진화적으로 보존되어 있는지 확인하기 위한 100% Stacked Bar 플롯 생성.
#   3) [enable_save_outputs] 원터치 마스터 스위치로 파일/로그 저장 여부 완벽 제어.
################################################################################

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 1) 필수 패키지 목록 정의 및 자동 설치
required_packages <- c("dplyr", "ggplot2", "cowplot", "mclust", "RColorBrewer")
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
})
theme_set(theme_cowplot())

################################################################################
# [Pre-check] 사용 가능한 Taxonomy Level(분류군 계급) 확인용 진단 코드
################################################################################
# 아래 세 줄만 드래그해서 먼저 실행해보시면, 데이터셋에 존재하는 
# 분류군 컬럼명(예: "L2", "L3", "L4", "L5", "L6", "Taxonomy")을 확인할 수 있습니다.
# temp_df <- read.csv(file.path("/home/scott/EDM_16SV4_PA/04_Phase4_Output/01_Data_Integration", "Target_ASVs_Absolute_Abundance_Calculated.csv"), nrows = 5)
# print(colnames(temp_df)[grep("^L[0-9]|Taxon", colnames(temp_df))])

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택 (옵션: "median", "mean")
# - 의미: ASV의 10년 치 온도 상호작용 강도(IS)를 대변할 대푯값. (권장: median)
summary_method  <- "median"  

# [2] GMM 클러스터 강제 할당 (옵션: 정수, 기본값: 3)
# - 의미: 분석의 논리적 연속성을 위해 [Negative, Neutral, Positive] 3그룹 고정.
g_clusters      <- 3 

# [3] 분류군 시각화 해상도 (옵션: "L2", "L3", "L4", "L5", "L6")
# - 의미: 플롯에 표시할 분류군 계급 (L2=Phylum, L3=Class, L4=Order 등).
tax_level <- "L2" 

# [4] 시각적 통합 임계값 (단위: %, 권장: 2.0 ~ 5.0)
# - 의미: 이 비율 미만의 미소 분류군들은 모두 "Others"(회색)로 통폐합하여 가독성 향상.
abundance_threshold <- 2.0 

# [5] 0(Zero) 풍부도 제외 여부 (옵션: TRUE, FALSE)
# - 의미: 휴면기(Abundance = 0) 데이터를 평균 계산에 포함할지 여부.
exclude_zeros <- FALSE

# [6] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
# - 의미: TRUE일 경우 디렉토리 생성, 로그 기록, TIFF/CSV 파일 저장을 실제로 수행함.
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/06_Taxonomic_Composition")

# [원터치 제어] 마스터 스위치가 TRUE일 때만 디렉토리 생성
if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

zero_str  <- ifelse(exclude_zeros, "ZeroExc", "ZeroInc")
log_file  <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_", tax_level, "_", zero_str, "_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_Composition_", tax_level, "_", zero_str, ".tiff"))
file_csv  <- file.path(out_dir, paste0("Part6_", toupper(summary_method), "_Composition_", tax_level, "_", zero_str, ".csv"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  # [원터치 제어] 마스터 스위치가 TRUE일 때만 텍스트 파일에 로그 기록
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # -------------------------------------------------------------------
  # Section 1. Data Load & Merge (Robust Selection)
  # -------------------------------------------------------------------
  log_msg("Loading and Merging Data with Taxonomy...")
  
  if(!file.exists(file_abundance)) stop("Error: Abundance Input file not found.")
  if(!file.exists(file_temp_is)) stop("Error: IS Input file not found.")
  
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  # 이전 데이터 추출에 Taxonomy 정보(tax_level) 추가 확보
  df_ab <- df_ab_raw %>% 
    dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance, Taxon = all_of(tax_level))
  df_is <- df_is_raw %>% 
    dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  df_merged <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  if (exclude_zeros) {
    df_merged <- df_merged %>% filter(Absolute_Abundance > 0)
    zero_text <- "Zero-excluded"
  } else {
    zero_text <- "Zero-included"
  }
  
  # ASV별 고유 Taxonomy 딕셔너리 생성
  asv_tax_dict <- df_ab %>% dplyr::select(ASV, Taxon) %>% distinct()
  
  # -------------------------------------------------------------------
  # Section 2. Apply GMM Clustering (Inherited exactly from Part 4)
  # -------------------------------------------------------------------
  log_msg("Applying GMM clustering to map ASVs to logical groups...")
  
  asv_summary <- df_merged %>%
    group_by(ASV) %>%
    summarise(
      agg_IS   = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE)
    ) %>% 
    filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  
  cluster_means <- gmm_model$parameters$mean
  cluster_order <- order(cluster_means) 
  mapped_clusters <- factor(gmm_model$classification, 
                            levels = cluster_order, 
                            labels = c("1_Negative", "2_Neutral", "3_Positive"))
  
  asv_clusters <- asv_summary %>%
    mutate(Cluster = mapped_clusters) %>%
    dplyr::select(ASV, Cluster) %>%
    inner_join(asv_tax_dict, by = "ASV")
  
  # -------------------------------------------------------------------
  # Section 3. Calculate Taxonomic Composition per Cluster
  # -------------------------------------------------------------------
  log_msg(paste0("Calculating taxonomic composition at [", tax_level, "] level..."))
  
  # 각 클러스터 내에서 해당 분류군이 차지하는 ASV 개수(비율) 집계
  tax_composition <- asv_clusters %>%
    group_by(Cluster, Taxon) %>%
    summarise(ASV_Count = n(), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(Percentage = ASV_Count / sum(ASV_Count) * 100) %>%
    ungroup()
  
  # 임계값 미만 분류군 'Others'로 통합
  tax_composition_clean <- tax_composition %>%
    mutate(Taxon_Clean = ifelse(Percentage < abundance_threshold, "Others", Taxon)) %>%
    group_by(Cluster, Taxon_Clean) %>%
    summarise(Percentage = sum(Percentage), .groups = "drop")
  
  # 'Others'가 범례의 가장 마지막에 오도록 요인(Factor) 레벨 재정렬
  tax_levels <- unique(tax_composition_clean$Taxon_Clean)
  tax_levels <- c(sort(tax_levels[tax_levels != "Others"]), "Others")
  tax_composition_clean$Taxon_Clean <- factor(tax_composition_clean$Taxon_Clean, levels = tax_levels)
  
  # -------------------------------------------------------------------
  # Section 4. Visualization (100% Stacked Bar Chart)
  # -------------------------------------------------------------------
  log_msg("Generating Stacked Bar Plot...")
  
  # 분류군 개수에 따른 색상 팔레트 자동 생성
  color_count <- length(tax_levels)
  get_palette <- colorRampPalette(brewer.pal(8, "Set2"))
  my_colors <- get_palette(color_count)
  names(my_colors) <- tax_levels
  my_colors["Others"] <- "#d9d9d9" 
  
  caption_base <- sprintf(
    "Note: Clusters defined via 1D GMM (G=%d) on %s_IS (%s).\nTaxa with < %.1f%% abundance are grouped into 'Others'.", 
    g_clusters, tools::toTitleCase(summary_method), zero_text, abundance_threshold
  )
  
  p_tax <- ggplot(tax_composition_clean, aes(x = Cluster, y = Percentage, fill = Taxon_Clean)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = my_colors, name = paste0("Taxonomy (", tax_level, ")")) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = "Taxonomic Composition of Temperature Sensitivity Groups",
      subtitle = "Comparing Taxonomic Profiles of Psychrophilic, Opportunistic, and Thermophilic strategies",
      x = paste0("GMM Cluster (", tools::toTitleCase(summary_method), " IS)"),
      y = "Relative Proportion of ASVs (%)",
      caption = caption_base
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
      plot.caption = element_text(hjust = 0, size = 11, color = "gray30", margin = margin(t = 15)),
      axis.text.x = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # -------------------------------------------------------------------
  # Section 5. Save Results & Export (원터치 제어)
  # -------------------------------------------------------------------
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_tax, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    write.csv(tax_composition_clean, file_csv, row.names = FALSE)
    log_msg(paste("[SUCCESS] Plot and CSV saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot was generated in RStudio Viewer only (No files saved).")
  }
  
  print(p_tax)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

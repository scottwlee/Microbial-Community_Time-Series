################################################################################
# [Phase 5 - Part 5 (Version_2.1): Macro-ecological Traits (Seasonal Phenology)]
# 목적: 10년 치 전체 데이터를 월별(Month) 축으로 투영하여, 
#       GMM으로 분류된 세 그룹(Negative, Neutral, Positive)의 계절적 출현 패턴을 분석함.
# 특징: 
#   1) [자동 패키지 설치] 환경에 누락된 필수 라이브러리를 감지하고 자동 설치.
#   2) Cyclic Cubic Splines(bs="cc")가 적용된 비선형 GAM 트렌드를 통해
#      여름철 호열성 우점, 겨울철 호냉성 우점, 기회주의적 스파이크 패턴을 입증함.
#   3) [enable_save_outputs] 원터치 마스터 스위치로 파일/로그 저장 여부 완벽 제어.
#   4) 이전 분석과 일관성을 맞추기 위해 Y축 명칭에 'Absolute' 명시.
#   5) [색상 통일] 지정된 그룹별 고유 색상 및 글로벌 추세선 색상 일괄 적용.
################################################################################

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 1) 필수 패키지 목록 정의 및 자동 설치
required_packages <- c("dplyr", "ggplot2", "mgcv", "cowplot", "mclust", "lubridate")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

# 2) 패키지 로드
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
  library(lubridate)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택 (옵션: "median", "mean")
# - 의미: ASV의 10년 치 온도 상호작용 강도를 대변할 단일 값 결정.
summary_method  <- "median"  

# [2] GMM 클러스터 개수 (옵션: 정수, 기본값 3)
# - 의미: 데이터를 몇 개의 생태적 지위로 분할할 것인가.
g_clusters      <- 3 

# [3] 제로 풍부도 배제 여부 (옵션: TRUE, FALSE)
# - 의미: 휴면기(Abundance = 0) 데이터를 평균 계산에 포함할지 여부.
exclude_zeros   <- FALSE     

# [4] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
# - 의미: TRUE일 경우 디렉토리 생성, 로그 기록, TIFF 파일 저장을 실제로 수행함.
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 파일명 동적 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir    <- file.path(base_dir, "05_Phase5_Output/05_Seasonal_Phenology")

# [원터치 제어] 마스터 스위치가 TRUE일 때만 디렉토리 생성
if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

zero_str <- ifelse(exclude_zeros, "ZeroExc", "ZeroInc")
log_file  <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_", zero_str, "_Phenology_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part5_", toupper(summary_method), "_", zero_str, "_Phenology_GAM.tiff"))

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
  log_msg("Loading and Merging Data...")
  
  if(!file.exists(file_abundance)) stop("Error: Abundance Input file not found.")
  if(!file.exists(file_temp_is)) stop("Error: IS Input file not found.")
  
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  df_merged <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  if (exclude_zeros) {
    df_merged <- df_merged %>% filter(Absolute_Abundance > 0)
    zero_text <- "Zero-excluded"
  } else {
    zero_text <- "Zero-included"
  }
  
  # -------------------------------------------------------------------
  # Section 2. Apply GMM Clustering to Establish Groups
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
    dplyr::select(ASV, Cluster)
  
  # -------------------------------------------------------------------
  # Section 3. Monthly Aggregation & Integration
  # -------------------------------------------------------------------
  log_msg("Extracting seasonal variables and calculating monthly phenology...")
  
  df_phenology <- df_merged %>%
    inner_join(asv_clusters, by = "ASV") %>%
    mutate(
      Date = as.Date(Date),
      Month = as.numeric(format(Date, "%m"))
    )
  
  # 각 ASV별로 월간 평균(Mean) 출현량을 집계하여 오버플로팅(Overplotting) 방지
  df_monthly_mean <- df_phenology %>%
    group_by(ASV, Month, Cluster) %>%
    summarise(
      Monthly_Mean_Abundance = mean(Absolute_Abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Log10 치환을 적용하여 시각화 스케일 최적화
    mutate(Log_Abundance = log10(Monthly_Mean_Abundance + 1))
  
  # -------------------------------------------------------------------
  # Section 4. Visualization with Cyclic GAM Splines
  # -------------------------------------------------------------------
  log_msg("Generating Seasonal Phenology Plot with Cyclic GAM...")
  
  method_title <- tools::toTitleCase(summary_method)
  
  # [색상 설정 업데이트] 사용자가 지정한 그룹별 고유 색상 및 공통 추세선 색상 정의
  custom_colors <- c("1_Negative" = "#347433", "2_Neutral" = "#999999", "3_Positive" = "#DC2525")
  trend_color   <- "#0065F8"  
  ci_color      <- "#0065F8"  
  
  caption_base <- sprintf(
    "Note: Clusters defined via 1D GMM (G=%d) on Temp_IS.\nData points represent monthly mean absolute abundance per ASV (%s).\nGAM lines fitted with Cyclic Cubic Splines (bs='cc') to model continuous seasonal shifts.",
    g_clusters, zero_text
  )
  
  p_phenology <- ggplot(df_monthly_mean, aes(x = Month, y = Log_Abundance, color = Cluster, fill = Cluster)) +
    geom_jitter(alpha = 0.3, size = 1.5, width = 0.2) +
    # 해당 플롯은 군집별(Cluster)로 구분된 추세선을 그리므로 color 매핑에 의해 custom_colors가 자동 적용됨.
    # 만약 단일 글로벌 추세선이 추가될 경우 color = trend_color, fill = ci_color를 활용할 수 있음.
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), 
                alpha = 0.2, linewidth = 1.5) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(
      title = "Seasonal Phenology & Temporal Niche Partitioning",
      subtitle = paste0("Tracking 10-year microbial bloom dynamics across ", method_title, " IS Groups"),
      x = "Month of the Year",
      y = "Log10(Monthly Mean Absolute Abundance + 1)", 
      caption = caption_base
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray20", margin = margin(b = 15)),
      plot.caption = element_text(hjust = 0, size = 11, color = "gray30", margin = margin(t = 15)),
      panel.grid.minor.x = element_blank()
    )
  
  # -------------------------------------------------------------------
  # Section 5. Export (원터치 제어)
  # -------------------------------------------------------------------
  # [원터치 제어] 마스터 스위치가 TRUE일 때만 플롯 파일 저장 실행
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_phenology, device = "tiff", dpi = 600, width = 10, height = 8, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", file_plot))
  } else {
    log_msg("[SAFE MODE] Plot was generated in RStudio Viewer only (No file saved).")
  }
  
  print(p_phenology)
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

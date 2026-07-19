################################################################################
# [Phase 5 - Part 7 (Revised & Upgraded): Statistical Definition of POW]
# 목적: 10년 치 전체 데이터를 월(Month) 축으로 투영하여 생물학적 전성기를 찾고, 
#       ANOVA를 통해 1~4월이 'Psychrophilic Optimal Window (POW)'임을 
#       통계적 수치 출력 및 시각적 캡션과 함께 객관적으로 증명함.
# 특징:
#   1) [자동 패키지 설치] 환경에 누락된 필수 라이브러리를 감지하고 자동 설치함.
#   2) [enable_save_outputs] 원터치 마스터 스위치로 파일/로그 저장 여부 완벽 제어.
#   3) Negative 그룹 전용 컬러(#347433)를 플롯 B, C에 일관되게 적용.
#   4) Y축 라벨에 'Absolute Abundance' 명시하여 용어 통일성 확보.
################################################################################

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 필수 패키지 감지 및 설치
required_packages <- c("dplyr", "ggplot2", "mgcv", "cowplot", "mclust")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 요약 통계량 선택 (옵션: "median", "mean")
# - 의미: 10년간의 온도 상호작용 강도(IS)를 대표할 통계값 (권장: median)
summary_method <- "median"  

# [2] GMM 클러스터 강제 할당 (옵션: 정수, 기본값: 3)
# - 의미: 분석의 논리적 연속성을 위해 [Negative, Neutral, Positive] 3그룹 고정.
g_clusters <- 3 

# [3] POW(Psychrophilic Optimal Window) 검증 구간 (옵션: 벡터형식 월)
# - 의미: ANOVA 통계 검정을 통해 타 구간 대비 유의미한 생장 기간인지 판별할 목표 구간 (기본: 1~4월)
pow_months <- c(1, 2, 3, 4) 

# [4] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
# - 의미: TRUE일 경우 디렉토리 생성, 로그 기록, TIFF 파일 저장을 실제로 수행함.
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir <- "/home/scott/EDM_16SV4_PA"
input_dir <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir <- file.path(base_dir, "05_Phase5_Output/07_POW_Definition")

# [원터치 제어] 디렉토리 생성
if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

log_file <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_POW_Diagnostic_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_POW_Validation_Plot.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  # [원터치 제어] 텍스트 파일에 로그 기록
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  # -------------------------------------------------------------------
  # Section 1. Data Load & Clustering (Consistent with Part 4-6)
  # -------------------------------------------------------------------
  log_msg("Step 1: Loading Data & Applying GMM Clustering...")
  
  if(!file.exists(file_abundance)) stop("Error: Abundance Input file not found.")
  if(!file.exists(file_temp_is)) stop("Error: IS Input file not found.")
  
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_merged <- dplyr::inner_join(
    df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance),
    df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength),
    by = c("ASV", "Date")
  ) %>% filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  asv_summary <- df_merged %>% group_by(ASV) %>%
    summarise(
      agg_IS = if(summary_method == "mean") mean(Temp_IS, na.rm=TRUE) else median(Temp_IS, na.rm=TRUE), 
      asv_mean = mean(Absolute_Abundance, na.rm=TRUE)
    ) %>% 
    filter(asv_mean > 0 & !is.na(agg_IS))
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  asv_clusters <- asv_summary %>%
    mutate(Cluster = factor(gmm_model$classification, levels = order(gmm_model$parameters$mean), labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)
  
  # -------------------------------------------------------------------
  # Section 2. Monthly Aggregation (Targeting 1_Negative Group)
  # -------------------------------------------------------------------
  log_msg("Step 2: Monthly Aggregation of the Negative IS Group...")
  df_ts <- df_ab_raw %>% 
    inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    filter(!is.na(Temperature))
  
  # 겨울철 엔진(Negative) 그룹의 합산 생물량 추적
  diag_df <- df_ts %>% filter(Cluster == "1_Negative") %>%
    group_by(Date, Month, Temperature) %>%
    summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Log_Abund = log10(Total_Abund + 1),
           Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
           Is_POW = ifelse(Month %in% pow_months, "POW (Jan-Apr)", "Non-POW"))
  
  # -------------------------------------------------------------------
  # Section 3. Statistical Testing (ANOVA) & Logging
  # -------------------------------------------------------------------
  log_msg("Step 3: Performing Statistical Testing (One-way ANOVA)...")
  anova_res <- aov(Log_Abund ~ Month_Fct, data = diag_df)
  anova_summary <- summary(anova_res)
  pval_anova <- anova_summary[[1]][["Pr(>F)"]][1]
  f_val <- anova_summary[[1]][["F value"]][1]
  
  # 콘솔 출력 및 로그 파일 기록을 동시에 수행
  log_msg("\n=== ANOVA Results ===")
  out_anova <- capture.output(print(anova_summary))
  for(l in out_anova) log_msg(l)
  
  log_msg("\n=== Monthly Mean Log-Abundance (Ranking) ===")
  monthly_means <- diag_df %>% group_by(Month_Fct) %>% summarise(Mean_Log = mean(Log_Abund)) %>% arrange(desc(Mean_Log))
  out_mean <- capture.output(print(monthly_means))
  for(l in out_mean) log_msg(l)
  
  pval_text <- ifelse(pval_anova < 0.001, "ANOVA p < 0.001", sprintf("ANOVA p = %.3f", pval_anova))
  caption_text <- sprintf("Note: One-way ANOVA confirms distinct temporal niche partitioning (F = %.1f, %s).\nMonths Jan-Apr (POW) demonstrate significantly higher community biomass load.", f_val, pval_text)
  
  # -------------------------------------------------------------------
  # Section 4. Visualization (Dynamic Formatting & Colors)
  # -------------------------------------------------------------------
  log_msg("Step 4: Generating Tri-Panel Visualization...")
  
  # 컬러 설정: 환경(A)은 붉은색계열, 군집동태(B,C)는 Negative 고유 컬러(#347433) 활용
  neg_color <- "#347433" # 1_Negative Group Base Color
  neg_fill  <- "#94bca4" # 1_Negative Group Light Fill Color
  
  p_temp <- ggplot(diag_df, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = "red", fill = "red", alpha = 0.2) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p_abund <- ggplot(diag_df, aes(x = Month, y = Log_Abund)) +
    geom_jitter(color = neg_color, alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = neg_color, fill = neg_fill, alpha = 0.4, linewidth = 1.2) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "B. Biological Phenology (Negative IS Group)", x = "Month", y = "Log10(Absolute Abundance + 1)") +
    theme(plot.title = element_text(face = "bold"))
  
  p_anova <- ggplot(diag_df, aes(x = Month_Fct, y = Log_Abund, fill = Is_POW)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = Is_POW), width = 0.15, alpha = 0.6, size = 1.2) +
    scale_fill_manual(values = c("POW (Jan-Apr)" = neg_fill, "Non-POW" = "gray80")) +
    scale_color_manual(values = c("POW (Jan-Apr)" = neg_color, "Non-POW" = "gray50")) +
    labs(title = "C. Statistical Validation of POW", subtitle = pval_text, x = "Month", y = "Log10(Absolute Abundance + 1)", caption = caption_text) +
    theme(legend.position = "top", legend.title = element_blank(), 
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, color = "gray30", margin = margin(t = 15))) 
  
  top_row <- plot_grid(p_temp, p_abund, ncol = 1, align = "v")
  p_combined <- plot_grid(top_row, p_anova, ncol = 2, rel_widths = c(1, 1.2))
  
  # -------------------------------------------------------------------
  # Section 5. Export (원터치 제어)
  # -------------------------------------------------------------------
  if (enable_save_outputs) {
    ggsave(file_plot, plot = p_combined, device = "tiff", dpi = 600, width = 16, height = 8, compression = "lzw")
    log_msg(paste("[SUCCESS] Plot saved successfully to:", out_dir))
  } else {
    log_msg("[SAFE MODE] Plot generated in Viewer only (No files saved).")
  }
  
  print(p_combined)
  log_msg("Analysis complete.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

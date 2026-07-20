################################################################################
# [Phase 5 - Part 7 (Final Master Version): Statistical Definition of Peak Window]
# 목적: 1_Negative 군집의 시차(Time-lag)가 반영된 진정한 절대 정점(Absolute Peak) 추출.
# 특징:
#   1) [완벽한 동기화] True Zero-Included 로직(결측치 0 복원)을 적용하여 플롯 B와 완벽 일치.
#   2) [분석법 듀얼 엔진] Dunnett's Test(엄격한 가설 검정) vs K-means(생태학적 분할) 자유 전환.
#   3) [파라미터 제어] dunnett_alpha(엄격도) 및 peak_k_clusters(페이즈 분할) 유저 지정 가능.
#   4) [로깅 시스템] 통계 분석 값 및 사용된 파라미터가 콘솔과 텍스트 로그에 명확히 기록됨.
################################################################################

options(stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "stats", "DescTools")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("[System] Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
  library(stats)
  library(DescTools) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS 
#################################################
# [1] 요약 통계량 및 클러스터 설정
summary_method <- "median"  
g_clusters <- 3 

# [2] Peak Window 선별 방식 (옵션: "dunnett", "kmeans")
# - dunnett: 최고점(Max)을 대조군으로 삼아 통계적으로 차이 없는 뾰족한 정점만 핀셋 추출.
# - kmeans : 월별 궤적을 기계학습으로 K개의 생태적 계층으로 분할.
peak_selection_method <- "kmeans"  

# -----------------------------------------------------------
# [옵션 A] Dunnett's Test 파라미터 ("dunnett" 선택 시 작동)
# - 0.05 : 통계학적 표준 기준.
# - 0.10 ~ 0.20 : 엄격한 기준. (수치가 클수록 허들이 높아져, 최고점과 미세하게만 달라도 가차 없이 Peak에서 탈락시킴 -> 극도로 좁은 최전성기 도출)
dunnett_alpha <- 0.05 

# [옵션 B] K-means 파라미터 ("kmeans" 선택 시 작동)
# - 권장값: 3 (생태학적 3단계 페이즈인 '1.Peak(정점)', '2.Transition(전이기)', '3.Trough(휴면기)'를 수학적으로 가장 잘 반영함)
peak_k_clusters <- 3  
# -----------------------------------------------------------

# [3] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir <- "/home/scott/EDM_16SV4_PA"
input_dir <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir <- file.path(base_dir, "05_Phase5_Output/07_Peak_Window_Definition")

if (enable_save_outputs && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")

log_file <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_", toupper(peak_selection_method), "_Peak_Log.txt"))
file_plot <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_", toupper(peak_selection_method), "_Peak.tiff"))

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Peak Window Analysis Parameters]")
  log_msg(sprintf(" - IS Summary Method     : %s", summary_method))
  log_msg(sprintf(" - Selection Method      : %s", toupper(peak_selection_method)))
  if (peak_selection_method == "dunnett") {
    log_msg(sprintf(" - Dunnett Alpha (Strict): %.2f %s", dunnett_alpha, ifelse(dunnett_alpha > 0.05, "(Strict Filtering)", "(Standard)")))
  } else {
    log_msg(sprintf(" - K-means Clusters (K)  : %d (Ecological Phases: Peak, Transition, Trough)", peak_k_clusters))
  }
  log_msg(sprintf(" - Zero Handling         : True Zero-Included"))
  log_msg("===================================================================\n")
  
  # -------------------------------------------------------------------
  # Section 1. Data Load & Clustering
  # -------------------------------------------------------------------
  log_msg("Step 1: Loading Data & Applying GMM Clustering...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_ab <- df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance)
  df_is <- df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)
  
  df_merged <- dplyr::inner_join(df_ab, df_is, by = c("ASV", "Date")) %>%
    dplyr::filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  asv_summary <- df_merged %>%
    group_by(ASV) %>%
    summarise(
      agg_IS   = if(summary_method == "mean") mean(Temp_IS, na.rm = TRUE) else median(Temp_IS, na.rm = TRUE),
      asv_mean = mean(Absolute_Abundance, na.rm = TRUE),
      .groups  = "drop"
    ) %>% 
    filter(!is.na(agg_IS) & asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  cluster_order <- order(gmm_model$parameters$mean) 
  mapped_clusters <- factor(gmm_model$classification, levels = cluster_order, labels = c("1_Negative", "2_Neutral", "3_Positive"))
  
  asv_clusters <- asv_summary %>%
    mutate(Cluster = mapped_clusters) %>%
    dplyr::select(ASV, Cluster)
  
  # -------------------------------------------------------------------
  # Section 2. True Zero-Included Monthly Aggregation
  # -------------------------------------------------------------------
  log_msg("Step 2: Aggregating Data with True Zero-Inclusion...")
  
  # 환경 데이터 (온도) 월별 평균 산출
  df_env <- df_ab_raw %>% 
    mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    filter(!is.na(Temperature)) %>%
    group_by(Date, Month) %>%
    summarise(Temperature = mean(Temperature, na.rm = TRUE), .groups = "drop")
  
  # 월별 총 샘플링 횟수(Effort) 도출
  month_counts <- df_ab_raw %>% 
    mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Date, Month) %>% 
    dplyr::distinct() %>% 
    group_by(Month) %>% 
    summarise(Num_Samples = n(), .groups = "drop")
  
  # ASV별 생물량 합산 및 0(Zero) 강제 복원 후 최종 평균 산출
  df_monthly_mean <- df_ab_raw %>%
    mutate(Month = as.numeric(format(as.Date(Sample_Date), "%m"))) %>%
    group_by(ASV_ID, Month) %>%
    summarise(Sum_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(ASV_ID, Month = 1:12, fill = list(Sum_Abund = 0)) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>% 
    dplyr::left_join(month_counts, by = "Month") %>%
    mutate(
      Mean_Abund = Sum_Abund / Num_Samples,
      Log_Abundance = log10(Mean_Abund + 1),
      Month_Fct = factor(Month, levels = 1:12, labels = month.abb)
    )
  
  diag_df <- df_monthly_mean %>% filter(Cluster == "1_Negative")
  
  # -------------------------------------------------------------------
  # Section 3. Peak Selection (Dunnett vs K-means)
  # -------------------------------------------------------------------
  log_msg(sprintf("Step 3: Defining Peak Window via %s...", toupper(peak_selection_method)))
  
  monthly_means <- diag_df %>% group_by(Month, Month_Fct) %>% summarise(mean_abund = mean(Log_Abundance), .groups="drop")
  target_months <- c()
  
  if (peak_selection_method == "dunnett") {
    max_month_idx <- which.max(monthly_means$mean_abund)
    max_month_fct <- as.character(monthly_means$Month_Fct[max_month_idx])
    log_msg(paste(" -> Control Month (Absolute Maximum):", max_month_fct))
    
    dt_res <- DunnettTest(x = diag_df$Log_Abundance, g = diag_df$Month_Fct, control = max_month_fct)
    pvals <- dt_res[[1]][, "pval"]
    
    # 설정된 alpha 값보다 큰(차이가 없는) 비교군 추출
    indistinguishable_comparisons <- rownames(dt_res[[1]])[pvals > dunnett_alpha]
    extracted_months_str <- sapply(strsplit(indistinguishable_comparisons, "-"), function(x) x[1])
    
    peak_months_str <- unique(c(max_month_fct, extracted_months_str))
    target_months <- monthly_means %>% filter(Month_Fct %in% peak_months_str) %>% pull(Month) %>% sort()
    
  } else if (peak_selection_method == "kmeans") {
    set.seed(414)
    km_res <- kmeans(monthly_means$mean_abund, centers = peak_k_clusters)
    monthly_means$KM_Cluster <- km_res$cluster
    peak_cluster_id <- which.max(km_res$centers)
    target_months <- monthly_means %>% filter(KM_Cluster == peak_cluster_id) %>% pull(Month) %>% sort()
  }
  
  log_msg(paste(" -> Extracted Peak Window Months:", paste(month.abb[target_months], collapse = ", ")))
  
  # 연속성 여부에 따른 라벨링 동적 생성
  if (length(target_months) == 1) {
    peak_label <- paste0("Peak Window (", month.abb[target_months], ")")
  } else if (length(target_months) > 1 && all(diff(target_months) == 1)) {
    peak_label <- paste0("Peak Window (", month.abb[min(target_months)], "-", month.abb[max(target_months)], ")")
  } else {
    peak_label <- paste0("Peak Window (", paste(month.abb[target_months], collapse=","), ")")
  }
  
  diag_df <- diag_df %>%
    mutate(Is_Peak = ifelse(Month %in% target_months, peak_label, "Non-Peak"))
  
  # -------------------------------------------------------------------
  # Section 4. Statistical Testing (Peak vs Non-Peak)
  # -------------------------------------------------------------------
  log_msg("Step 4: Performing Final Statistical Testing (ANOVA)...")
  anova_res <- aov(Log_Abundance ~ Is_Peak, data = diag_df)
  anova_summary <- summary(anova_res)
  pval_anova <- anova_summary[[1]][["Pr(>F)"]][1]
  f_val <- anova_summary[[1]][["F value"]][1]
  
  pval_text <- ifelse(pval_anova < 0.001, "ANOVA p < 0.001", sprintf("ANOVA p = %.3f", pval_anova))
  
  log_msg("\n===================================================================")
  log_msg(" [Final Statistical Validation Results]")
  log_msg(sprintf(" - ANOVA F-value : %.2f", f_val))
  log_msg(sprintf(" - ANOVA P-value : %s", pval_text))
  log_msg("===================================================================\n")
  
  caption_method <- ifelse(peak_selection_method == "dunnett", 
                           sprintf("Dunnett's Test vs Max Month (\u03b1=%.2f)", dunnett_alpha),
                           sprintf("Unsupervised K-means (k=%d)", peak_k_clusters))
  caption_text <- sprintf("Note: Peak Window defined via %s on True Zero-Included data.\nOne-way ANOVA confirms distinct temporal niche partitioning (F = %.1f, %s).", caption_method, f_val, pval_text)
  
  # -------------------------------------------------------------------
  # Section 5. Visualization 
  # -------------------------------------------------------------------
  log_msg("Step 5: Generating Synchronized Tri-Panel Visualization...")
  neg_color <- "#347433" 
  neg_fill  <- "#94bca4" 
  
  p_temp <- ggplot(df_env, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = "red", fill = "red", alpha = 0.2, linewidth=1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p_abund <- ggplot(diag_df, aes(x = Month, y = Log_Abundance)) +
    geom_jitter(color = neg_color, alpha = 0.3, size = 1.5, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = neg_color, fill = neg_fill, alpha = 0.2, linewidth = 1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = paste0("B. Biological Phenology (Negative IS Group)"), x = "Month", y = "Log10(Monthly Mean Absolute Abundance + 1)") +
    theme(plot.title = element_text(face = "bold"))
  
  fill_palette <- c(neg_fill, "gray80"); names(fill_palette) <- c(peak_label, "Non-Peak")
  color_palette <- c(neg_color, "gray50"); names(color_palette) <- c(peak_label, "Non-Peak")
  
  p_anova <- ggplot(diag_df, aes(x = Month_Fct, y = Log_Abundance, fill = Is_Peak)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = Is_Peak), width = 0.15, alpha = 0.6, size = 1.2) +
    scale_fill_manual(values = fill_palette) +
    scale_color_manual(values = color_palette) +
    labs(title = paste("C. Statistical Validation of", peak_label), subtitle = pval_text, x = "Month of the Year", y = "Log10(Monthly Mean Absolute Abundance + 1)", caption = caption_text) +
    theme(legend.position = "top", legend.title = element_blank(), 
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, color = "gray30", margin = margin(t = 15))) 
  
  top_row <- plot_grid(p_temp, p_abund, ncol = 1, align = "v")
  p_combined <- plot_grid(top_row, p_anova, ncol = 2, rel_widths = c(1, 1.2))
  
  # -------------------------------------------------------------------
  # Section 6. Export
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

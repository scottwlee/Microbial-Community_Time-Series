################################################################################
# [Phase 5 - Part 7 & 9 (Ultimate Master Version): Auto-Linked POW & Climate Impact]
# 목적: Tukey's HSD 그룹핑(Letter Grouping)을 통해 통계적으로 완벽하게 정의된 
#       최전성기(POW)를 추출하고, 이를 Part 9의 시차 보정 분석으로 자동 연동함.
# 특징:
#   1) [Tukey HSD 그룹핑] 'a' 그룹(최고 생물량과 통계적 차이가 없는 달)을 자동 추출.
#   2) [Two-way ANOVA] 10년 시계열 특성을 반영하여 연도(Year)를 블록 효과로 통제.
#   3) [정밀 True Zero] 샘플링 미수행(NA)과 미발견(0)을 철저히 구분하는 시계열 복원.
#   4) [패키지 도입] 그룹핑을 위해 R 농생물통계 표준 패키지인 `agricolae` 사용.
################################################################################

options(stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# Section 0. Environment Setup & Package Auto-Installation
# -------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "cowplot", "mclust", "stats", "agricolae")
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
  library(agricolae) 
})
theme_set(theme_cowplot())

#################################################
# USER SETTINGS (스위치 및 파라미터 제어)
#################################################
# [1] 공통 GMM 클러스터 설정
summary_method <- "median"  
g_clusters     <- 3 
target_cluster <- "1_Negative" 

# [2] [Part 7] Peak Window 선별 방식 (옵션: "tukey", "kmeans")
# - tukey: Two-way ANOVA 후 Tukey's HSD 사후검정 실시 -> 최상위 'a' 그룹 추출
# - kmeans: 월별 궤적을 순수 기계학습(거리 기반)으로 K개의 생태적 계층으로 분할
peak_selection_method <- "tukey"  

# [2-A] K-means 파라미터 ("kmeans" 선택 시 작동)
peak_k_clusters <- 3  

# [3] 결과 저장 마스터 스위치 (옵션: TRUE, FALSE)
enable_save_outputs <- FALSE

# -------------------------------------------------------------------
# 경로 및 동적 파일명 설정
# -------------------------------------------------------------------
base_dir   <- "/home/scott/EDM_16SV4_PA"
input_dir  <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
smap_dir   <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap")

out_dir_p7 <- file.path(base_dir, "05_Phase5_Output/07_Peak_Window_Definition")
out_dir_p9 <- file.path(base_dir, "05_Phase5_Output/10_Lag_Adjusted_Impact")

if (enable_save_outputs) {
  if (!dir.exists(out_dir_p7)) dir.create(out_dir_p7, recursive = TRUE)
  if (!dir.exists(out_dir_p9)) dir.create(out_dir_p9, recursive = TRUE)
}

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is   <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_smap_sum  <- file.path(smap_dir, "Phase3_Part3_MDR_Smap_Summary.csv")

log_file_p7  <- file.path(out_dir_p7, paste0("Part7_", toupper(summary_method), "_", toupper(peak_selection_method), "_Peak_Log.txt"))
plot_file_p7 <- file.path(out_dir_p7, paste0("Part7_", toupper(summary_method), "_", toupper(peak_selection_method), "_Peak.tiff"))

log_file_p9  <- file.path(out_dir_p9, paste0("Part9_", toupper(summary_method), "_POW_Impact_Log.txt"))
plot_file_p9 <- file.path(out_dir_p9, paste0("Part9_", toupper(summary_method), "_POW_Impact.tiff"))
csv_file_p9  <- file.path(out_dir_p9, paste0("Part9_", toupper(summary_method), "_POW_Impact_Data.csv"))

log_msg <- function(msg, target_log = NULL) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  if (enable_save_outputs && !is.null(target_log)) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = target_log, append = TRUE)
  }
}

tryCatch({
  log_msg("\n===================================================================")
  log_msg(" [Integrated Pipeline: POW Definition -> Lag-Adjusted Impact]")
  log_msg(sprintf(" - Target Cluster        : %s", target_cluster))
  log_msg(sprintf(" - Peak Selection Method : %s", toupper(peak_selection_method)))
  log_msg("===================================================================\n")
  
  # -------------------------------------------------------------------
  # Section 1. Data Load & Common GMM Clustering 
  # -------------------------------------------------------------------
  log_msg("[COMMON] Step 1: Loading Data & Applying True Zero-Included GMM...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  df_smap_raw <- read.csv(file_smap_sum, stringsAsFactors = FALSE)
  
  colnames(df_smap_raw)[1] <- "ASV_Hash"
  asv_dict <- df_ab_raw %>% dplyr::select(ASV_Hash, ASV_ID) %>% dplyr::distinct()
  df_smap <- df_smap_raw %>% dplyr::inner_join(asv_dict, by = "ASV_Hash") %>% dplyr::select(ASV = ASV_ID, Best_TP)
  
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
  
  neg_color <- "#347433" 
  neg_fill  <- "#94bca4" 
  
  #####################################################################
  # [PART 7 EXECUTION] - Peak Window Definition (Cluster-level Block ANOVA)
  #####################################################################
  log_msg("\n>>> [PART 7] STARTING PEAK WINDOW DEFINITION <<<", log_file_p7)
  
  df_env <- df_ab_raw %>% 
    dplyr::mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(!is.na(Temperature)) %>%
    dplyr::group_by(Date, Month) %>%
    dplyr::summarise(Temperature = mean(Temperature, na.rm = TRUE), .groups = "drop")
  
  # 정밀 True Zero 로직: '실제 샘플링이 수행된 연도-월'만 기준 골격으로 생성
  sampled_ym <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::select(Year, Month) %>% dplyr::distinct()
  
  df_cluster_abund <- df_ab_raw %>%
    dplyr::mutate(Date = as.Date(Sample_Date), Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    dplyr::filter(Cluster == target_cluster) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop")
  
  # 샘플링은 되었으나 타겟 종이 안 나온 달(NA)을 0으로 맵핑하여 군집 단위 시계열 생성
  df_cluster_monthly <- sampled_ym %>%
    dplyr::left_join(df_cluster_abund, by = c("Year", "Month")) %>%
    dplyr::mutate(
      Total_Abund = tidyr::replace_na(Total_Abund, 0),
      Log_Abundance = log10(Total_Abund + 1),
      Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
      Year_Fct = factor(Year)
    )
  
  # Peak Selection Logic (Tukey's HSD or K-means)
  target_months <- c()
  
  if (peak_selection_method == "tukey") {
    # Two-way ANOVA (Randomized Complete Block Design: Year = Block)
    aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
    
    # Tukey's HSD with Compact Letter Display (CLD)
    tukey_res <- HSD.test(aov_model, "Month_Fct", group = TRUE)
    top_group <- tukey_res$groups
    
    log_msg("\n[Tukey's HSD Letter Grouping Results]", log_file_p7)
    for(i in 1:nrow(top_group)) {
      log_msg(sprintf(" -> %s: Group '%s' (Mean Log_Abund: %.2f)", rownames(top_group)[i], top_group$groups[i], top_group$Log_Abundance[i]), log_file_p7)
    }
    
    # 'a' 문자를 포함하는 모든 그룹(예: "a", "ab", "abc")을 Peak Window로 판정
    peak_months_str <- rownames(top_group)[grepl("a", top_group$groups)]
    target_months <- which(month.abb %in% peak_months_str)
    
    # 보고용 통계량
    f_val <- summary(aov_model)[[1]][["F value"]][1]
    pval_anova <- summary(aov_model)[[1]][["Pr(>F)"]][1]
    
  } else if (peak_selection_method == "kmeans") {
    monthly_means <- df_cluster_monthly %>% dplyr::group_by(Month, Month_Fct) %>% dplyr::summarise(mean_abund = mean(Log_Abundance), .groups="drop")
    set.seed(414)
    km_res <- kmeans(monthly_means$mean_abund, centers = peak_k_clusters)
    monthly_means$KM_Cluster <- km_res$cluster
    target_months <- monthly_means %>% dplyr::filter(KM_Cluster == which.max(km_res$centers)) %>% dplyr::pull(Month) %>% sort()
    
    aov_model <- aov(Log_Abundance ~ Month_Fct + Year_Fct, data = df_cluster_monthly)
    f_val <- summary(aov_model)[[1]][["F value"]][1]
    pval_anova <- summary(aov_model)[[1]][["Pr(>F)"]][1]
  }
  
  target_months <- sort(target_months)
  peak_months_str <- paste(month.abb[target_months], collapse = ", ")
  log_msg(sprintf("\n[Part 7] Extracted POW Months: %s", peak_months_str), log_file_p7)
  
  peak_label <- ifelse(length(target_months) == 1, paste0("Peak Window (", month.abb[target_months], ")"), 
                       ifelse(all(diff(target_months) == 1), paste0("Peak Window (", month.abb[min(target_months)], "-", month.abb[max(target_months)], ")"), 
                              paste0("Peak Window (", peak_months_str, ")")))
  
  df_cluster_monthly <- df_cluster_monthly %>% dplyr::mutate(Is_Peak = ifelse(Month %in% target_months, peak_label, "Non-Peak"))
  
  pval_text_p7 <- ifelse(pval_anova < 0.001, "p < 0.001", sprintf("p = %.3f", pval_anova))
  log_msg(sprintf("[Part 7 STATS] Block ANOVA (Month Effect) F-value: %.2f, %s", f_val, pval_text_p7), log_file_p7)
  
  # Part 7 Plotting (3-Panel)
  p7_temp <- ggplot(df_env, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = "red", fill = "red", alpha = 0.2, linewidth=1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p7_abund <- ggplot(df_cluster_monthly, aes(x = Month, y = Log_Abundance)) +
    geom_jitter(color = neg_color, alpha = 0.5, size = 1.5, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = neg_color, fill = neg_fill, alpha = 0.2, linewidth = 1.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = paste0("B. Biological Phenology (", target_cluster, ")"), x = "Month", y = "Log10(Cluster Abundance + 1)") + theme(plot.title = element_text(face = "bold"))
  
  fill_palette <- c(neg_fill, "gray80"); names(fill_palette) <- c(peak_label, "Non-Peak")
  color_palette <- c(neg_color, "gray50"); names(color_palette) <- c(peak_label, "Non-Peak")
  
  caption_p7 <- sprintf("Note: Peak Window defined via %s.\nTwo-way ANOVA (Block: Year) confirms temporal niche partitioning (F = %.1f, %s).", toupper(peak_selection_method), f_val, pval_text_p7)
  p7_anova <- ggplot(df_cluster_monthly, aes(x = Month_Fct, y = Log_Abundance, fill = Is_Peak)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + geom_jitter(aes(color = Is_Peak), width = 0.15, alpha = 0.6, size = 1.2) +
    scale_fill_manual(values = fill_palette) + scale_color_manual(values = color_palette) +
    labs(title = paste("C. Statistical Validation of", peak_label), subtitle = paste("Two-way ANOVA (Month Effect):", pval_text_p7), x = "Month of the Year", y = "Log10(Cluster Abundance + 1)", caption = caption_p7) +
    theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, color = "gray30", margin = margin(t = 15))) 
  
  p7_combined <- plot_grid(plot_grid(p7_temp, p7_abund, ncol = 1, align = "v"), p7_anova, ncol = 2, rel_widths = c(1, 1.2))
  if(enable_save_outputs) { ggsave(plot_file_p7, plot = p7_combined, device = "tiff", dpi = 600, width = 16, height = 8, compression = "lzw") }
  print(p7_combined)
  
  #####################################################################
  # [LINKAGE BRIDGE] 
  #####################################################################
  eco_window_months <- target_months 
  log_msg(sprintf("\n>>> [LINKAGE] POW (%s) automatically passed to Part 9 <<<", paste(eco_window_months, collapse=",")), log_file_p9)
  
  #####################################################################
  # [PART 9 EXECUTION] - Lag-Adjusted Climate Impact
  #####################################################################
  log_msg("\n>>> [PART 9] STARTING LAG-ADJUSTED CLIMATE IMPACT <<<", log_file_p9)
  
  df_env_daily <- df_ab_raw %>% dplyr::select(Date = Sample_Date, Temperature) %>% dplyr::distinct() %>% dplyr::filter(!is.na(Temperature)) %>% dplyr::arrange(Date)
  
  df_aligned <- df_ab %>%
    tidyr::complete(ASV, Date = all_sample_dates, fill = list(Absolute_Abundance = 0)) %>%
    dplyr::arrange(ASV, as.Date(Date)) %>% 
    dplyr::inner_join(asv_clusters, by = "ASV") %>% 
    dplyr::inner_join(df_smap, by = "ASV") %>%
    dplyr::mutate(Response_Lag_Weeks = abs(Best_TP)) %>%
    dplyr::group_by(ASV) %>%
    dplyr::mutate(shift_n = as.integer(Response_Lag_Weeks[1]), Shifted_Abundance = dplyr::lead(Absolute_Abundance, n = shift_n[1])) %>%
    dplyr::ungroup() %>% dplyr::filter(!is.na(Shifted_Abundance)) %>%
    dplyr::inner_join(df_env_daily, by = "Date") %>% dplyr::mutate(Date = as.Date(Date))
  
  analysis_df <- df_aligned %>%
    dplyr::mutate(Year = as.numeric(format(Date, "%Y")), Month = as.numeric(format(Date, "%m"))) %>%
    dplyr::filter(Month %in% eco_window_months & Cluster == target_cluster) %>% 
    dplyr::group_by(Date, Year, Temperature) %>%
    dplyr::summarise(Total_Aligned_Abund = sum(Shifted_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Log_Abund = log10(Total_Aligned_Abund + 1), Year_Fct = as.factor(Year))
  
  pval_abund <- summary(lm(Log_Abund ~ Year, data = analysis_df))$coefficients[2, 4]
  pval_direct <- summary(mgcv::gam(Log_Abund ~ s(Temperature, k=6), data = analysis_df, method="REML"))$s.table[1, 4]
  format_pval <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  
  log_msg(sprintf("[Part 9 STATS] Abundance vs Year (Decline)   : %s", format_pval(pval_abund)), log_file_p9)
  log_msg(sprintf("[Part 9 STATS] Temp vs Abundance (GAM Limit) : %s", format_pval(pval_direct)), log_file_p9)
  
  p9_abund <- ggplot(analysis_df, aes(x = Year, y = Log_Abund)) +
    geom_boxplot(aes(group=Year_Fct), fill=neg_fill, alpha=0.3, color=neg_color, outlier.shape=NA) +
    geom_jitter(color=neg_color, alpha=0.6, width=0.15) + geom_smooth(method="lm", color=neg_color, fill=neg_fill, alpha=0.3) +
    labs(title="A. Suppression of Psychrophilic Engine", subtitle=paste("Decline of Phase-Aligned Abundance (", format_pval(pval_abund), ")"), x="Year", y="Log10(Lag-Adjusted POW Abundance + 1)") +
    scale_x_continuous(breaks = min(analysis_df$Year):max(analysis_df$Year)) + theme(plot.title = element_text(face="bold"))
  
  caption_p9 <- sprintf("Note: Analyzed strictly within the Dynamically Linked POW (Months %s).\nAbundances are phase-aligned using True Zero-Included Best_TP.", paste(eco_window_months, collapse="-"))
  p9_direct <- ggplot(analysis_df, aes(x = Temperature, y = Log_Abund)) +
    geom_point(color=neg_color, size=2.5, alpha=0.6) + geom_smooth(method="gam", formula=y~s(x, bs="cs", k=6), color="black", linetype="dashed", linewidth=1.2) +
    labs(title="B. High-Resolution POW Thermal Limitation", subtitle=paste("GAM correlation between Trigger Temp & Abundance (", format_pval(pval_direct), ")"), x="POW Trigger Temp (°C)", y="Log10(Lag-Adjusted POW Abundance + 1)", caption=caption_p9) +
    theme(plot.title = element_text(face="bold"), plot.caption = element_text(hjust=0, color="gray30", margin=margin(t=15)))
  
  p9_combined <- plot_grid(p9_abund, p9_direct, ncol=2, align="h")
  
  if(enable_save_outputs) { 
    ggsave(plot_file_p9, plot = p9_combined, device = "tiff", dpi = 600, width = 13, height = 6, compression = "lzw") 
    write.csv(analysis_df, csv_file_p9, row.names = FALSE)
  }
  print(p9_combined)
  
  log_msg("\n[SUCCESS] Integrated Pipeline Completed Successfully.")
  
}, error = function(e) { log_msg(paste("ERROR:", e$message)); stop(e) })

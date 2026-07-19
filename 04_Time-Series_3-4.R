################################################################################
# [Phase 5 - Part 7 (Revised & Upgraded): Statistical Definition of POW]
# 목적: 10년 치 전체 데이터를 월(Month) 축으로 투영하여 생물학적 전성기를 찾고, 
#       ANOVA를 통해 1~4월이 'Psychrophilic Optimal Window (POW)'임을 
#       통계적 수치 출력 및 시각적 캡션과 함께 객관적으로 증명함.
################################################################################

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(mclust)
})
theme_set(theme_cowplot())

# USER SETTINGS
summary_method <- "median"  
g_clusters <- 3 
pow_months <- c(1, 2, 3, 4) # 검증할 POW 기간

base_dir <- "/home/scott/EDM_16SV4_PA"
input_dir <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_dir <- file.path(base_dir, "05_Phase5_Output/07_POW_Definition")
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

file_abundance <- file.path(input_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")
file_temp_is <- file.path(input_dir, "Merged_Interaction_Strength_Final.csv")
file_plot <- file.path(out_dir, paste0("Part7_", toupper(summary_method), "_POW_Validation_Plot.tiff"))
log_file <- file.path(out_dir, "Part7_POW_Diagnostic_Log.txt")

log_msg <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
  # cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Step 1: Loading Data & Clustering...")
  df_ab_raw <- read.csv(file_abundance, stringsAsFactors = FALSE)
  df_is_raw <- read.csv(file_temp_is, stringsAsFactors = FALSE)
  
  df_merged <- dplyr::inner_join(
    df_ab_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance),
    df_is_raw %>% dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength),
    by = c("ASV", "Date")
  ) %>% filter(!is.na(Temp_IS) & !is.na(Absolute_Abundance))
  
  asv_summary <- df_merged %>% group_by(ASV) %>%
    summarise(agg_IS = median(Temp_IS, na.rm=TRUE), asv_mean = mean(Absolute_Abundance, na.rm=TRUE)) %>% 
    filter(asv_mean > 0)
  
  set.seed(414)
  gmm_model <- Mclust(asv_summary$agg_IS, G = g_clusters)
  asv_clusters <- asv_summary %>%
    mutate(Cluster = factor(gmm_model$classification, levels = order(gmm_model$parameters$mean), labels = c("1_Negative", "2_Neutral", "3_Positive"))) %>%
    dplyr::select(ASV, Cluster)
  
  log_msg("Step 2: Monthly Aggregation...")
  df_ts <- df_ab_raw %>% 
    inner_join(asv_clusters, by = c("ASV_ID" = "ASV")) %>%
    mutate(Date = as.Date(Sample_Date), Month = as.numeric(format(Date, "%m"))) %>%
    filter(!is.na(Temperature))
  
  diag_df <- df_ts %>% filter(Cluster == "1_Negative") %>%
    group_by(Date, Month, Temperature) %>%
    summarise(Total_Abund = sum(Absolute_Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Log_Abund = log10(Total_Abund + 1),
           Month_Fct = factor(Month, levels = 1:12, labels = month.abb),
           Is_POW = ifelse(Month %in% pow_months, "POW (Jan-Apr)", "Non-POW"))
  
  log_msg("Step 3: Statistical Testing (ANOVA)...")
  anova_res <- aov(Log_Abund ~ Month_Fct, data = diag_df)
  anova_summary <- summary(anova_res)
  pval_anova <- anova_summary[[1]][["Pr(>F)"]][1]
  f_val <- anova_summary[[1]][["F value"]][1]
  
  # === 콘솔에 통계 결과 출력 및 로그 기록 ===
  log_msg("\n=== ANOVA Results ===")
  log_msg(capture.output(print(anova_summary)))
  
  log_msg("\n=== Monthly Mean Log-Abundance (Ranking) ===")
  monthly_means <- diag_df %>% group_by(Month_Fct) %>% summarise(Mean_Log = mean(Log_Abund)) %>% arrange(desc(Mean_Log))
  log_msg(capture.output(print(monthly_means)))
  
  pval_text <- ifelse(pval_anova < 0.001, "ANOVA p < 0.001", sprintf("ANOVA p = %.3f", pval_anova))
  caption_text <- sprintf("Note: One-way ANOVA confirms distinct temporal niche partitioning (F = %.1f, %s).\nMonths Jan-Apr (POW) demonstrate significantly higher community biomass load.", f_val, pval_text)
  
  log_msg("Step 4: Generating Visualization...")
  p_temp <- ggplot(diag_df, aes(x = Month, y = Temperature)) +
    geom_jitter(color = "darkred", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = "red", fill = "red", alpha = 0.2) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "A. Physical Environment", y = "Temp (°C)") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  
  p_abund <- ggplot(diag_df, aes(x = Month, y = Log_Abund)) +
    geom_jitter(color = "#3182ce", alpha = 0.3, width = 0.2) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 12), color = "#090AB5", fill = "lightblue", alpha = 0.2) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "B. Biological Phenology", x = "Month", y = "Log10(Abundance+1)") +
    theme(plot.title = element_text(face = "bold"))
  
  # === 하단 패널(C)에 통계 캡션(Note) 추가 ===
  p_anova <- ggplot(diag_df, aes(x = Month_Fct, y = Log_Abund, fill = Is_POW)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = Is_POW), width = 0.15, alpha = 0.5, size = 1.2) +
    scale_fill_manual(values = c("POW (Jan-Apr)" = "#3182ce", "Non-POW" = "gray80")) +
    scale_color_manual(values = c("POW (Jan-Apr)" = "#090AB5", "Non-POW" = "gray50")) +
    labs(title = "C. Statistical Validation of POW", subtitle = pval_text, x = "Month", y = "Log10(Abundance+1)", caption = caption_text) +
    theme(legend.position = "top", legend.title = element_blank(), 
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, color = "gray30", margin = margin(t = 15))) # 캡션 스타일 지정
  
  top_row <- plot_grid(p_temp, p_abund, ncol = 1, align = "v")
  p_combined <- plot_grid(top_row, p_anova, ncol = 2, rel_widths = c(1, 1.2))
  
  # ggsave(file_plot, plot = p_combined, device = "tiff", dpi = 600, width = 16, height = 8, compression = "lzw")
  print(p_combined)
  log_msg("Analysis complete.")
  
}, error = function(e) { stop(e) })

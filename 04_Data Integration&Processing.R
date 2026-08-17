################################################################################
# [Phase 4 Master Script] 
# 목적: Phase 3의 S-map 계수들을 자동 인식 및 통합하고, 
#       이에 맞춰 절대 풍부도를 산출하여 Phase 5(거시 생태학 분석)용 데이터를 완성합니다.
# 
# [Master Pipeline 목차]
#   Section 1. Environment & Path Settings (환경 및 경로 설정)
#   Section 2. Interaction Strength Data Integration (상호작용 강도 병합)
#   Section 3. Time-Series Visualization of Interaction Strength (IS 시계열 플롯)
#   Section 4. Absolute Abundance Calculation (절대 풍부도 산출)
#   Section 5. Time-Series Visualization of Absolute Abundance (풍부도 시계열 플롯 - True Zero 적용)
#   Section 6. Moving Window CV Calculation & Data Merging (CV 연산 및 데이터 병합)
#   Section 7. Scatter Plot (CV vs Interaction Strength) (최종 스캐터 플롯)
################################################################################

options(stringsAsFactors = FALSE)

################################################################################
# Section 1. Environment & Path Settings
################################################################################
cat("[INFO] Checking and installing required packages for Phase 4...\n")
required_packages <- c("dplyr", "tidyr", "ggplot2", "scales", "lubridate", "stringr", 
                       "zoo", "mgcv", "cowplot", "parallel", "pbapply")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
  library(phyloseq)
})
theme_set(theme_cowplot())

base_dir       <- "/home/scott/EDM_16SV4_PA"

input_dir_coef <- file.path(base_dir, "03_Phase3_Output/Phase3_Part3_MDR_Smap/Coefficients")
input_ps_raw   <- file.path(base_dir, "02_Phase2_Output/Bac_Phyloseq_Raw_Updated.rds")
input_ps_filt  <- file.path(base_dir, "02_Phase2_Output/Filtered_Phyloseq_Updated.rds")

out_data_dir   <- file.path(base_dir, "04_Phase4_Output/01_Data_Integration")
out_plot_dir   <- file.path(base_dir, "04_Phase4_Output/02_Visualization")
out_cv_dir     <- file.path(base_dir, "04_Phase4_Output/03_Part4_Moving_Window_CV")
out_scat_dir   <- file.path(base_dir, "04_Phase4_Output/04_Part5_Scatter_Plot")

for(d in c(out_data_dir, out_plot_dir, out_cv_dir, out_scat_dir)) {
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- file.path(base_dir, "04_Phase4_Output", paste0("Phase4_Master_Log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
write_log <- function(msg) {
  timestamp_msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(timestamp_msg, "\n")
  cat(timestamp_msg, "\n", file = log_file, append = TRUE)
}

write_log("=== Phase 4 Master Pipeline Started ===")

n_threads  <- min(33, parallel::detectCores())
start_date <- "2012-03-07"
x_date_range <- c("2012-01-01", "2021-12-31")

################################################################################
# Section 2. Interaction Strength Data Integration
################################################################################
write_log("[Section 2] Interaction Strength Data Integration")

file_temp_is <- file.path(out_data_dir, "Merged_Interaction_Strength_Final.csv")

ps_filt <- readRDS(input_ps_filt)
tax_df <- as.data.frame(tax_table(ps_filt))
tax_df$ASV_Hash <- rownames(tax_df)

mapping_df <- tax_df %>% 
  dplyr::select(ASV_Hash, ASV) %>% 
  dplyr::rename(ASV_ID = ASV)

extract_interaction_data <- function(file) {
  data_raw <- read.csv(file)
  cols <- colnames(data_raw)
  if(!("time" %in% cols)) return(NULL)
  c_cols <- grep("^c_[1-9][0-9]*$", cols, value = TRUE)
  if(length(c_cols) == 0) return(NULL)
  c_nums <- as.numeric(gsub("c_", "", c_cols))
  max_c_col <- paste0("c_", max(c_nums))
  asv_hash <- gsub("Coef_|.csv|_2", "", basename(file)) 
  
  data_raw %>%
    dplyr::select(time, dplyr::all_of(max_c_col)) %>%
    dplyr::rename(Interaction_Strength = dplyr::all_of(max_c_col)) %>%
    dplyr::mutate(ASV_Hash = asv_hash) %>%
    dplyr::filter(!is.na(Interaction_Strength))
}

file_list <- sort(list.files(input_dir_coef, pattern = "\\.csv$", full.names = TRUE))
write_log(sprintf("Extracting coefficients from %d dynamically detected files using %d threads...", length(file_list), n_threads))

cl <- makeCluster(n_threads)
clusterEvalQ(cl, library(dplyr))
merged_list <- pblapply(file_list, extract_interaction_data, cl = cl)
stopCluster(cl)

mapped_data <- dplyr::bind_rows(Filter(Negate(is.null), merged_list)) %>%
  dplyr::mutate(Sample_Date = as.Date(start_date) + weeks(time - 1)) %>%
  dplyr::left_join(mapping_df, by = "ASV_Hash") %>%
  dplyr::mutate(
    ASV_num = as.numeric(str_extract(ASV_ID, "\\d+")),
    ASV_ID = reorder(factor(ASV_ID), ASV_num)
  ) %>%
  dplyr::select(Sample_Date, time, ASV_ID, ASV_Hash, Interaction_Strength)

write.csv(mapped_data, file_temp_is, row.names = FALSE)
write_log("-> Saved: Merged_Interaction_Strength_Final.csv")

################################################################################
# Section 3. Time-Series Visualization of Interaction Strength
################################################################################
write_log("[Section 3] Time-Series Visualization (Temp_IS)")

years <- year(as.Date(x_date_range[1])):year(as.Date(x_date_range[2]))
shading_ranges <- data.frame(
  xmin = as.Date(paste0(years, "-01-01")), 
  xmax = as.Date(paste0(years + 1, "-01-01")), 
  ymin = -Inf, ymax = Inf, year = years
) %>% dplyr::filter(year %% 2 == 0)

p_is <- ggplot(mapped_data, aes(x = Sample_Date, y = Interaction_Strength, color = ASV_ID)) +
  geom_rect(data = shading_ranges, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="grey85", alpha=0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.0) +
  geom_vline(xintercept = as.Date(x_date_range[1]), color = "black", linewidth = 1.0) +
  geom_line(linewidth = 0.3, alpha = 0.6) +
  scale_x_date(limits = as.Date(x_date_range), date_breaks = "2 years", date_labels = "%Y", expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 1)) +
  labs(x = "Sample Date", y = "Interaction Strength", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_line(color="black"), legend.position="bottom",
        legend.key.width = unit(1,"cm"), legend.text=element_text(size=9)) +
  guides(color = guide_legend(ncol=10, override.aes=list(linewidth=1, alpha=1)))

output_is_plot_path <- file.path(out_plot_dir, "Phase4_Section3_Interaction_Strength_Plot.tiff")
ggsave(output_is_plot_path, plot = p_is, device = "tiff", dpi = 600, width = 14, height = 8, compression = "lzw")
write_log(paste("-> Saved:", output_is_plot_path))

################################################################################
# Section 4. Absolute Abundance Calculation
################################################################################
write_log("[Section 4] Absolute Abundance Calculation")

file_abundance <- file.path(out_data_dir, "Target_ASVs_Absolute_Abundance_Calculated.csv")

asv_mapping_df <- read.csv(file_temp_is) %>%
  dplyr::select(ASV_Hash, ASV_ID) %>%
  dplyr::distinct()

ps_raw <- readRDS(input_ps_raw)
total_reads_raw <- sample_sums(ps_raw)
meta_raw <- data.frame(sample_data(ps_raw))
meta_raw$Sample_ID <- rownames(meta_raw)
meta_raw$Total_Raw_Reads <- total_reads_raw[rownames(meta_raw)]
meta_raw$FC_Total <- meta_raw$Syn_abundance + meta_raw$HB_abundance
sample_calc_info <- meta_raw %>% dplyr::select(Sample_ID, Total_Raw_Reads, FC_Total)

ps_target <- prune_taxa(asv_mapping_df$ASV_Hash, ps_filt)
df_target_melt <- psmelt(ps_target) %>% 
  dplyr::rename(ASV_Hash = OTU, Raw_ASV_Count = Abundance) %>% 
  dplyr::select(-Sample)

abund_data <- df_target_melt %>%
  dplyr::left_join(sample_calc_info, by = "Sample_ID") %>%
  dplyr::left_join(asv_mapping_df, by = "ASV_Hash") %>%
  dplyr::mutate(
    Absolute_Abundance = (Raw_ASV_Count / Total_Raw_Reads) * FC_Total,
    Sample_Date = as.Date(as.character(Date)),
    ASV_Num = as.numeric(str_extract(ASV_ID, "\\d+")),
    ASV_ID = reorder(factor(ASV_ID), ASV_Num)
  ) %>%
  dplyr::filter(!is.na(Absolute_Abundance), !is.na(Sample_Date))

write.csv(abund_data, file_abundance, row.names = FALSE)
write_log("-> Saved: Target_ASVs_Absolute_Abundance_Calculated.csv")

################################################################################
# Section 5. Time-Series Visualization of Absolute Abundance
################################################################################
write_log("[Section 5] Time-Series Visualization of Absolute Abundance")

# Phase 5의 "True Zero" 로직을 시각화에도 동일하게 적용하여 빈 날짜를 0으로 강제 렌더링
all_sample_dates <- unique(abund_data$Sample_Date)

abund_data_plot <- abund_data %>%
  tidyr::complete(ASV_ID, Sample_Date = all_sample_dates, fill = list(Absolute_Abundance = 0))

p_ab <- ggplot(abund_data_plot, aes(x = Sample_Date, y = Absolute_Abundance, color = ASV_ID)) +
  geom_rect(data = shading_ranges, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, fill="grey85", alpha=0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.0) +
  geom_vline(xintercept = as.Date(x_date_range[1]), color = "black", linewidth = 1.0) +
  geom_line(linewidth = 0.3, alpha = 0.6) +
  scale_x_date(limits = as.Date(x_date_range), date_breaks = "2 years", date_labels = "%Y", expand = c(0,0)) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Sample Date", y = "Absolute Abundance (cells/mL)", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_line(color="black"), legend.position="bottom",
        legend.key.width=unit(1,"cm"), legend.text=element_text(size=9)) +
  guides(color = guide_legend(ncol=10, override.aes=list(linewidth=1, alpha=1)))

output_ab_plot_path <- file.path(out_plot_dir, "Phase4_Section5_Absolute_Abundance_Plot.tiff")
ggsave(output_ab_plot_path, plot = p_ab, device = "tiff", dpi = 600, width = 14, height = 8, compression = "lzw")
write_log(paste("-> Saved:", output_ab_plot_path))

################################################################################
# Section 6. Moving Window CV Calculation & Data Merging
################################################################################
write_log("[Section 6] Moving Window CV Calculation & Data Merging")

window_size <- 26; min_valid_points <- 13; min_nonzero <- 2
cv_fun_strict <- function(x) {
  x_valid <- x[!is.na(x)]
  if (length(x_valid) < min_valid_points) return(NA_real_)
  if (sum(x_valid > 0) < min_nonzero) return(NA_real_)
  m <- mean(x_valid)
  if (is.na(m) || m == 0) return(NA_real_)
  return(sd(x_valid) / m)
}

abund_for_cv <- abund_data %>% 
  dplyr::select(ASV = ASV_ID, Date = Sample_Date, Absolute_Abundance) %>% 
  dplyr::arrange(ASV, Date)

asv_list <- split(abund_for_cv, abund_for_cv$ASV)

cl <- makeCluster(n_threads)
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, varlist = c("window_size", "cv_fun_strict", "min_valid_points", "min_nonzero"))

calc_cv_list <- pblapply(asv_list, function(sub_df) {
  sub_df$CV <- rollapply(sub_df$Absolute_Abundance, width=window_size, FUN=cv_fun_strict, align="left", fill=NA_real_)
  return(sub_df %>% dplyr::select(ASV, Date, CV))
}, cl = cl)
stopCluster(cl)

df_cv <- dplyr::bind_rows(calc_cv_list)

df_temp_is <- mapped_data %>% 
  dplyr::select(ASV = ASV_ID, Date = Sample_Date, Temp_IS = Interaction_Strength)

df_merged <- dplyr::inner_join(df_cv, df_temp_is, by = c("ASV", "Date")) %>%
  dplyr::filter(!is.na(CV) & !is.na(Temp_IS)) %>% 
  dplyr::filter(CV > 0) 

q_temp <- quantile(df_merged$Temp_IS, probs = c(0.01, 0.99), na.rm = TRUE)
q_cv   <- quantile(df_merged$CV, probs = c(0.01, 0.99), na.rm = TRUE)
df_final <- df_merged %>% 
  dplyr::filter(Temp_IS >= q_temp[1] & Temp_IS <= q_temp[2]) %>% 
  dplyr::filter(CV >= q_cv[1] & CV <= q_cv[2])

write.csv(df_final, file.path(out_cv_dir, "Phase4_Section6_Merged_Filtered_CV_IS.csv"), row.names = FALSE)
write_log(sprintf("-> Saved Merged CV Data (%d rows).", nrow(df_final)))

################################################################################
# Section 7. Scatter Plot (CV vs Interaction Strength)
################################################################################
write_log("[Section 7] Scatter Plot (CV vs Interaction Strength)")

df_final_plot <- df_final %>%
  dplyr::mutate(ASV_Num = as.numeric(str_extract(ASV, "\\d+")), ASV = reorder(factor(ASV), ASV_Num))

p_scatter <- ggplot(data = df_final_plot, aes(x = Temp_IS, y = CV, color = ASV)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Distribution of Interaction Strength and Abundance CV",
       x = "Interaction Strength", y = "Abundance CV", color = NULL) +
  theme_cowplot() +
  theme(plot.title=element_text(face="bold", size=16), axis.title=element_text(face="bold"),
        panel.grid.major=element_line(color="grey90", linetype="dashed"),
        legend.position="bottom", legend.text=element_text(size=9)) +
  guides(color = guide_legend(ncol=10, override.aes=list(size=3, alpha=1)))

output_scatter_path <- file.path(out_scat_dir, "Phase4_Section7_Raw_Scatter_CV_vs_TempIS.tiff")
ggsave(output_scatter_path, plot = p_scatter, device = "tiff", dpi = 600, width = 14, height = 10, compression = "lzw")
write_log(paste("-> Saved:", output_scatter_path))

write_log("=== Phase 4 Master Pipeline Completed Successfully! ===")

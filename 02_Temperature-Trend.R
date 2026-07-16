################################################################################
# 마이크로바이옴 파이프라인 Phase 2: 시계열 분석 및 상관관계 (Publication Ready)
################################################################################

#################################################
# USER SETTINGS (사용자 설정 영역)
#################################################
# 1. 경로 설정 
# Input: Phase 1에서 생성된 단일 Output 폴더
# Output: 본 Phase에서 생성할 단일 Output 폴더
base_dir       <- "/home/scott/EDM_16SV4_PA"
input_dir      <- file.path(base_dir, "01_Phase1_Output")
out_dir        <- file.path(base_dir, "02_Phase2_Output")

# 2. 분석 파라미터
outlier_temp_threshold <- 35.0  # 온도 이상치 기준 (°C)
abs_abund_threshold    <- 1e6   # 절대 풍부도 상한선 임계값

# 3. 시스템 설정
n_threads <- min(33, parallel::detectCores())
save_plot <- TRUE

#################################################
# Section 1. Load Packages, Directories & Setup Logging
#################################################
required_packages <- c("phyloseq", "tidyverse", "dplyr", "ggplot2", "scales")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Output 디렉토리 단일 생성 (모든 결과물이 한 곳에 저장됨)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 로그 파일 초기화 및 파라미터 기록 (단일 Output 폴더 내)
log_file <- file.path(out_dir, paste0("phase2_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
write_log <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"), file = log_file, append = TRUE)
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

write_log("=== Phase 2: 시계열 데이터 교정 및 상관관계 분석 시작 ===")
write_log(paste("Input Directory:", input_dir))
write_log(paste("Output Directory:", out_dir))
write_log(paste("Outlier Temperature Threshold:", outlier_temp_threshold, "°C"))
write_log(paste("Absolute Abundance Threshold:", abs_abund_threshold))
write_log(paste("Threads assigned:", n_threads))

#################################################
# Section 2. Load Phyloseq Objects (Read-Only Input)
#################################################
# [데이터 무결성] 절대 오리지널 데이터 수정 금지: Phase 1의 결과물만 읽어옵니다.
file_raw <- file.path(input_dir, "Bac_Phyloseq_Raw.rds")
file_filtered <- file.path(input_dir, "Filtered_Phyloseq.rds")

if(!file.exists(file_raw)) stop("Error: Raw Phyloseq 파일을 찾을 수 없습니다. 경로를 확인하세요.")
if(!file.exists(file_filtered)) stop("Error: Filtered Phyloseq 파일을 찾을 수 없습니다. 경로를 확인하세요.")

ps_raw <- readRDS(file_raw)
ps_filtered <- readRDS(file_filtered)
write_log("필터링 전(Raw)과 후(Filtered) 원본 Phyloseq 객체 로드 완료.")

#################################################
# Section 3. Metadata Preprocessing & Outlier Correction
#################################################
# 베이스라인이 될 Raw 데이터를 함께 받아, 정확한 절대 풍부도와 보정 온도를 계산
update_phyloseq_metadata <- function(physeq_obj, baseline_obj) {
  meta <- data.frame(sample_data(physeq_obj))
  
  if(!"Sample_ID" %in% colnames(meta)) {
    meta$Sample_ID <- rownames(meta)
  }
  
  meta$Date <- as.Date(substr(gsub("Pro", "", meta$Sample_ID), 1, 8), format = "%Y%m%d")
  meta <- meta %>% arrange(Date)
  
  # 온도 이상치(outlier) 보정 로직
  temp_vec <- meta$Temperature
  for (i in 2:(length(temp_vec) - 1)) {
    if (!is.na(temp_vec[i]) && temp_vec[i] >= outlier_temp_threshold) {
      prev_val <- temp_vec[i - 1]
      next_val <- temp_vec[i + 1]
      if (!is.na(prev_val) && !is.na(next_val)) {
        temp_vec[i] <- mean(c(prev_val, next_val), na.rm = TRUE)
      }
    }
  }
  # 보정된 온도를 메타데이터에 덮어쓰기
  meta$Temperature <- temp_vec
  
  # 군집(ASV)의 완벽한 절대 풍부도 역산출
  total_reads_raw <- sample_sums(baseline_obj)
  current_reads <- sample_sums(physeq_obj)
  aligned_current <- current_reads[meta$Sample_ID]
  aligned_baseline <- total_reads_raw[meta$Sample_ID]
  
  proportion <- aligned_current / aligned_baseline
  fc_total <- meta$Syn_abundance + meta$HB_abundance
  meta$ASV_Total_Reads <- proportion * fc_total
  
  # 보정된 온도와 풍부도가 Phyloseq 객체 내부에 영구 반영됨
  rownames(meta) <- meta$Sample_ID
  sample_data(physeq_obj) <- sample_data(meta)
  
  return(physeq_obj)
}

write_log("메타데이터(날짜, 온도 이상치 평활화, 실제 ASV 절대 풍부도) 교정 중...")

ps_raw_updated <- update_phyloseq_metadata(ps_raw, ps_raw)
ps_filtered_updated <- update_phyloseq_metadata(ps_filtered, ps_raw)

# 업데이트된 객체를 단일 Output 폴더에 저장 (향후 분석용)
saveRDS(ps_raw_updated, file = file.path(out_dir, "Bac_Phyloseq_Raw_Updated.rds"))
saveRDS(ps_filtered_updated, file = file.path(out_dir, "Filtered_Phyloseq_Updated.rds"))
write_log("교정된 메타데이터(보정 온도 포함)가 반영된 새로운 Phyloseq 파일 저장 완료.")

#################################################
# Section 4. Combined Dual-Axis Time Series Visualization
#################################################
write_log("이중 축 시계열 패턴 시각화 (Combined Plot) 생성 중...")

# 패키지 충돌 방지를 위한 Base R 방식 데이터 추출 및 병합
df_raw <- as(sample_data(ps_raw_updated), "data.frame")
df_raw <- df_raw[, c("Date", "Temperature", "ASV_Total_Reads")]
colnames(df_raw)[3] <- "Reads_Raw"

df_filt <- as(sample_data(ps_filtered_updated), "data.frame")
df_filt <- df_filt[, c("Date", "ASV_Total_Reads")]
colnames(df_filt)[2] <- "Reads_Filt"

df_combined <- merge(df_raw, df_filt, by = "Date", all = FALSE)
df_combined <- na.omit(df_combined)

# 유사도 지표 계산
cor_pearson <- cor(df_combined$Reads_Raw, df_combined$Reads_Filt, method = "pearson")
cor_spearman <- cor(df_combined$Reads_Raw, df_combined$Reads_Filt, method = "spearman")
mean_prop <- mean(df_combined$Reads_Filt / df_combined$Reads_Raw, na.rm = TRUE)

# 스케일링 팩터 계산
max_temp <- max(df_combined$Temperature, na.rm = TRUE)
max_abund <- max(df_combined$Reads_Raw, na.rm = TRUE)
scale_factor <- max_abund / max_temp

# Combined Plot 생성
p_time_combined <- ggplot(df_combined, aes(x = Date)) +
  geom_line(aes(y = Temperature, color = "Temperature"), linewidth = 0.6) +
  geom_line(aes(y = Reads_Raw / scale_factor, color = "Pre-Filtering Abundance"), linewidth = 0.6, alpha = 0.8) +
  geom_line(aes(y = Reads_Filt / scale_factor, color = "Post-Filtering (Core) Abundance"), linewidth = 0.7, alpha = 0.9) +
  
  scale_y_continuous(
    name = "Temperature (°C)",
    sec.axis = sec_axis(~ . * scale_factor, name = "Total ASV Abundance", labels = scientific_format(digits = 2))
  ) +
  
  # Legend 표기 순서 고정 및 색상 매핑
  scale_color_manual(
    name = "Legend",
    breaks = c("Temperature", "Pre-Filtering Abundance", "Post-Filtering (Core) Abundance"),
    values = c(
      "Temperature" = "red",
      "Pre-Filtering Abundance" = "blue",
      "Post-Filtering (Core) Abundance" = "#D55E00"
    )
  ) +
  
  labs(
    title = "Combined Time Series of Temperature and Total Abundance",
    x = "Date",
    caption = sprintf(
      "Note:\n1. Pearson correlation (Pre vs Post): %.4f (%.1f%%)\n2. Spearman correlation (Pre vs Post): %.4f (%.1f%%)\n3. Mean retained biomass proportion: %.1f%%",
      cor_pearson, cor_pearson * 100, cor_spearman, cor_spearman * 100, mean_prop * 100
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left = element_text(color = "red", face = "bold"),
    axis.title.y.right = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0, face = "italic", color = "grey30", size = 11, margin = margin(t = 15))
  )

# 즉시 출력 보장 및 단일 Output 폴더에 저장
print(p_time_combined)

if(save_plot) {
  ggsave(file.path(out_dir, "01_TimeSeries_Combined.tiff"), plot = p_time_combined, width = 11, height = 6, dpi = 300)
}
write_log("Combined 이중 축 시계열 시각화 완료.")

#################################################
# Section 5. Correlation Analysis (Temp vs Abundance)
#################################################
fmt_rp <- function(cor_test_obj) {
  r  <- unname(cor_test_obj$estimate)
  p  <- cor_test_obj$p.value
  paste0("r = ", sprintf("%.3f", r), ", p = ", 
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 2), sprintf("%.4f", p)))
}

analyze_and_plot_correlation <- function(physeq_obj, condition_name) {
  df <- data.frame(sample_data(physeq_obj)) %>% select(Temperature, ASV_Total_Reads) %>% na.omit()
  
  plot_cor <- function(data, title, rp_text, filename) {
    p <- ggplot(data, aes(x = Temperature, y = ASV_Total_Reads)) +
      geom_point(alpha = 0.6, color = "darkblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = paste(title, "-", condition_name), subtitle = rp_text, x = "Temperature (°C)", y = "Total ASV Abundance") +
      scale_y_continuous(labels = scientific_format(digits = 2)) +
      theme_classic(base_size = 14) +
      theme(
        axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13)
      ) +
      annotate("text", x = -Inf, y = Inf, label = rp_text, hjust = -0.1, vjust = 1.5, size = 5, color = "black", fontface = "bold")
    
    print(p)
    if(save_plot) ggsave(file.path(out_dir, filename), plot = p, width = 7, height = 5, dpi = 300)
  }
  
  # A) Original Correlation
  cor_A <- cor.test(df$Temperature, df$ASV_Total_Reads, method = "pearson")
  plot_cor(df, "Correlation (Original)", fmt_rp(cor_A), paste0("Corr_A_Original_", condition_name, ".tiff"))
  
  # B) IQR Method
  Q1 <- quantile(df$ASV_Total_Reads, 0.25, na.rm = TRUE)
  Q3 <- quantile(df$ASV_Total_Reads, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  df_IQR <- df %>% filter(ASV_Total_Reads > (Q1 - 1.5 * IQR_value), ASV_Total_Reads < (Q3 + 1.5 * IQR_value))
  
  if(nrow(df_IQR) > 2) {
    cor_B <- cor.test(df_IQR$Temperature, df_IQR$ASV_Total_Reads, method = "pearson")
    plot_cor(df_IQR, "Correlation (IQR Filtered)", fmt_rp(cor_B), paste0("Corr_B_IQR_", condition_name, ".tiff"))
  }
  
  # C) Absolute Threshold Method
  df_abs <- df %>% filter(ASV_Total_Reads < abs_abund_threshold)
  if(nrow(df_abs) > 2) {
    cor_C <- cor.test(df_abs$Temperature, df_abs$ASV_Total_Reads, method = "pearson")
    plot_cor(df_abs, paste0("Correlation (Total < ", format(abs_abund_threshold, scientific=TRUE), ")"), 
             fmt_rp(cor_C), paste0("Corr_C_Threshold_", condition_name, ".tiff"))
  }
}

write_log("상관관계 분석 및 산점도 출력 중...")
analyze_and_plot_correlation(ps_raw_updated, "Pre-Filtering")
analyze_and_plot_correlation(ps_filtered_updated, "Post-Filtering")

#################################################
# Section 6. Final Results & Session Info Logging
#################################################
write_log("분석 환경(패키지 및 세션 정보)을 저장합니다.")
sink(file.path(out_dir, "session_info_phase2.txt"))
print(sessionInfo())
sink()

write_log("=== Phase 2 분석 정상 종료 ===")

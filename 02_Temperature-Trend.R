################################################################################
# 마이크로바이옴 파이프라인 Phase 2: 시계열 분석 및 환경 변수 트렌드
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

# 3. 시스템 및 출력 설정
n_threads   <- min(33, parallel::detectCores())
save_output <- TRUE       # TRUE: 데이터, 플롯, 로그 파일 모두 저장 / FALSE: 저장 없이 출력만 실행
plot_style  <- "single"   # "dual": 이중 축(Dual-axis) 통합 플롯 / "single": 온도와 풍부도를 각각 분리된 플롯으로 생성

#################################################
# Section 1. Load Packages, Directories & Setup Logging
#################################################
# [신규 추가] 2단 구획 박스 생성을 위한 ggh4x 패키지 추가
required_packages <- c("phyloseq", "tidyverse", "dplyr", "ggplot2", "scales", "lubridate", "ggh4x")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Output 디렉토리 생성
if(save_output) {
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  log_file <- file.path(out_dir, paste0("phase2_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
}

# 로그 함수
write_log <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
  if(save_output) {
    cat(paste0("[", Sys.time(), "] ", message, "\n"), file = log_file, append = TRUE)
  }
}

write_log("=== Phase 2: 시계열 데이터 교정 및 분석 시작 ===")
write_log(paste("Input Directory:", input_dir))
write_log(paste("Output Directory:", out_dir))
write_log(paste("Output Saving Enabled:", save_output))

#################################################
# Section 2. Load Phyloseq Objects (Read-Only Input)
#################################################
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
update_phyloseq_metadata <- function(physeq_obj, baseline_obj) {
  meta <- data.frame(sample_data(physeq_obj))
  
  if(!"Sample_ID" %in% colnames(meta)) {
    meta$Sample_ID <- rownames(meta)
  }
  
  meta$Date <- as.Date(substr(gsub("Pro", "", meta$Sample_ID), 1, 8), format = "%Y%m%d")
  meta <- meta %>% arrange(Date)
  
  # 온도 이상치 보정 로직
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
  meta$Temperature <- temp_vec
  
  # 군집(ASV)의 완벽한 절대 풍부도 역산출
  total_reads_raw <- sample_sums(baseline_obj)
  current_reads <- sample_sums(physeq_obj)
  aligned_current <- current_reads[meta$Sample_ID]
  aligned_baseline <- total_reads_raw[meta$Sample_ID]
  
  proportion <- aligned_current / aligned_baseline
  fc_total <- meta$Syn_abundance + meta$HB_abundance
  meta$ASV_Total_Reads <- proportion * fc_total
  
  rownames(meta) <- meta$Sample_ID
  sample_data(physeq_obj) <- sample_data(meta)
  
  return(physeq_obj)
}

write_log("메타데이터(날짜, 온도 이상치 평활화, 실제 ASV 절대 풍부도) 교정 중...")

ps_raw_updated <- update_phyloseq_metadata(ps_raw, ps_raw)
ps_filtered_updated <- update_phyloseq_metadata(ps_filtered, ps_raw)

if(save_output) {
  saveRDS(ps_raw_updated, file = file.path(out_dir, "Bac_Phyloseq_Raw_Updated.rds"))
  saveRDS(ps_filtered_updated, file = file.path(out_dir, "Filtered_Phyloseq_Updated.rds"))
  write_log("교정된 메타데이터가 반영된 새로운 Phyloseq 파일 저장 완료.")
}

#################################################
# Section 4. Time Series Visualization (Dual or Single Axis)
#################################################
write_log("연도별 교차 음영 및 트렌드가 적용된 시계열 시각화 생성 중...")

df_raw <- as(sample_data(ps_raw_updated), "data.frame")
df_raw <- df_raw[, c("Date", "Temperature", "ASV_Total_Reads")]
colnames(df_raw)[3] <- "Reads_Raw"

df_filt <- as(sample_data(ps_filtered_updated), "data.frame")
df_filt <- df_filt[, c("Date", "ASV_Total_Reads")]
colnames(df_filt)[2] <- "Reads_Filt"

df_combined <- merge(df_raw, df_filt, by = "Date", all = FALSE)
df_combined <- na.omit(df_combined)

# -------------------------------------------------------------
# [통계 검정] 전체 기간 수온 트렌드 분석 (Anomaly-based Linear Regression)
# 계절성(Seasonality)에 의한 오차를 제거하여 정확한 장기 트렌드 p-value를 산출함.
# -------------------------------------------------------------
df_combined$DecYear <- decimal_date(df_combined$Date)

# 1. 월별 기준점(Baseline) 및 수온 편차(Anomaly) 계산
df_combined <- df_combined %>%
  mutate(Month_Fct = format(Date, "%m")) %>%
  group_by(Month_Fct) %>%
  mutate(
    Monthly_Mean_Temp = mean(Temperature, na.rm = TRUE),
    Temp_Anomaly = Temperature - Monthly_Mean_Temp
  ) %>%
  ungroup()

# 2. 계절성이 제거된 '수온 편차(Anomaly)'를 이용한 엄격한 선형 회귀
fit_temp_anomaly <- lm(Temp_Anomaly ~ DecYear, data = df_combined)
slope_temp <- coef(fit_temp_anomaly)[2]
pval_temp <- summary(fit_temp_anomaly)$coefficients[2, 4]

# 국제 표준 유의성 별표 
get_stars_global <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("") 
}
sig_stars <- get_stars_global(pval_temp)

cat("\n==================================================\n")
cat("[Overall Temperature Trend Statistics (De-seasonalized Anomaly)]\n")
cat(sprintf("Annual Slope : %+.3f °C / year %s\n", slope_temp, sig_stars))
cat(sprintf("P-value      : %.3e\n", pval_temp))
cat("==================================================\n\n")
write_log(sprintf("Overall Temp Trend (Anomaly): %+.3f °C/year %s, p=%.3e", slope_temp, sig_stars, pval_temp))

write_log(sprintf("Overall Temp Trend: %+.3f °C/year, p=%.3e", slope_temp, pval_temp))

cor_pearson <- cor(df_combined$Reads_Raw, df_combined$Reads_Filt, method = "pearson")
mean_prop <- mean(df_combined$Reads_Filt / df_combined$Reads_Raw, na.rm = TRUE)

note_text_combined <- sprintf(
  "Statistics Note:\n- Overall Temp Trend: %+.3f °C/yr (p = %.3e)\n- Pearson correlation (Pre vs Post): %.4f (%.1f%%)\n- Mean retained biomass proportion: %.1f%%\n\nPlot Elements:\n- Grey shaded backgrounds indicate even years.\n- Solid GREEN line (#347433) represents the overall linear regression trend.\n- GREEN shaded area around the solid line represents the 95%% Confidence Interval (CI).",
  slope_temp, pval_temp, cor_pearson, cor_pearson * 100, mean_prop * 100
)

note_text_temp <- sprintf(
  "Statistics Note:\n- Overall Temp Trend: %+.3f °C/yr (p = %.3e)\n\nPlot Elements:\n- Grey shaded backgrounds indicate even years.\n- Solid GREEN line (#347433) represents the overall linear regression trend.\n- GREEN shaded area represents the 95%% Confidence Interval (CI).",
  slope_temp, pval_temp
)

note_text_abund <- sprintf(
  "Statistics Note:\n- Pearson correlation (Pre vs Post): %.4f (%.1f%%)\n- Mean retained biomass proportion: %.1f%%\n\nPlot Elements:\n- Grey shaded backgrounds indicate even years for temporal reference.",
  cor_pearson, cor_pearson * 100, mean_prop * 100
)

write_log("\n==================================================")
write_log("[LOGGED PLOT NOTES: TIME SERIES]")
write_log("--- Plot: Combined Time Series ---")
for (line in strsplit(note_text_combined, "\n")[[1]]) if(trimws(line) != "") write_log(line)
write_log("==================================================\n")

min_year <- as.numeric(format(min(df_combined$Date), "%Y"))
max_year <- as.numeric(format(max(df_combined$Date), "%Y"))
years_seq <- min_year:max_year

shading_ranges <- data.frame(
  xmin = as.Date(paste0(years_seq, "-01-01")),
  xmax = as.Date(paste0(years_seq + 1, "-01-01")),
  ymin = -Inf,
  ymax = Inf,
  year = years_seq
) %>% filter(year %% 2 == 0)

my_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

if(plot_style == "dual") {
  max_temp <- max(df_combined$Temperature, na.rm = TRUE)
  max_abund <- max(df_combined$Reads_Raw, na.rm = TRUE)
  scale_factor <- max_abund / max_temp
  
  p_time <- ggplot(df_combined, aes(x = Date)) +
    geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_smooth(aes(y = Temperature), method = "lm", color = "#347433", fill = "#347433", alpha = 0.15, linetype = "solid", se = TRUE, linewidth = 0.8) +
    geom_line(aes(y = Temperature, color = "Temperature"), linewidth = 0.6) +
    geom_line(aes(y = Reads_Raw / scale_factor, color = "Pre-Filtering Abundance"), linewidth = 0.6, alpha = 0.8) +
    geom_line(aes(y = Reads_Filt / scale_factor, color = "Post-Filtering (Core) Abundance"), linewidth = 0.7, alpha = 0.9) +
    scale_x_date(expand = c(0, 0)) +
    scale_y_continuous(name = "Temperature (°C)", sec.axis = sec_axis(~ . * scale_factor, name = "Total ASV Abundance", labels = scientific_format(digits = 2))) +
    scale_color_manual(breaks = c("Temperature", "Pre-Filtering Abundance", "Post-Filtering (Core) Abundance"), values = c("Temperature" = "red", "Pre-Filtering Abundance" = "black", "Post-Filtering (Core) Abundance" = "#D55E00")) +
    labs(title = "Combined Time Series of Temperature and Total Abundance", x = "Sample Date") + 
    my_theme +
    theme(axis.title.y.left = element_text(color = "red", face = "bold"), axis.title.y.right = element_text(color = "black", face = "bold"))
  
  print(p_time)
  if(save_output) ggsave(file.path(out_dir, "01_TimeSeries_DualAxis.tiff"), plot = p_time, width = 11, height = 7, dpi = 300)
  
} else if (plot_style == "single") {
  p_temp <- ggplot(df_combined, aes(x = Date)) +
    geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_smooth(aes(y = Temperature), method = "lm", color = "#347433", fill = "#347433", alpha = 0.15, linetype = "solid", se = TRUE, linewidth = 0.8) +
    geom_line(aes(y = Temperature), color = "red", linewidth = 0.6) +
    scale_x_date(expand = c(0, 0)) +
    scale_y_continuous(name = "Temperature (°C)") +
    labs(title = "Time Series of Temperature", x = "Sample Date") + 
    my_theme +
    theme(axis.title.y = element_text(color = "red", face = "bold"))
  
  p_abund <- ggplot(df_combined, aes(x = Date)) +
    geom_rect(data = shading_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
    geom_line(aes(y = Reads_Raw, color = "Pre-Filtering Abundance"), linewidth = 0.6, alpha = 0.8) +
    geom_line(aes(y = Reads_Filt, color = "Post-Filtering (Core) Abundance"), linewidth = 0.7, alpha = 0.9) +
    scale_x_date(expand = c(0, 0)) +
    scale_y_continuous(name = "Total ASV Abundance", labels = scientific_format(digits = 2)) +
    scale_color_manual(breaks = c("Pre-Filtering Abundance", "Post-Filtering (Core) Abundance"), values = c("Pre-Filtering Abundance" = "black", "Post-Filtering (Core) Abundance" = "#D55E00")) +
    labs(title = "Time Series of Total Abundance", x = "Sample Date") + 
    my_theme +
    theme(axis.title.y = element_text(color = "black", face = "bold"))
  
  print(p_temp)
  print(p_abund)
  if(save_output) {
    ggsave(file.path(out_dir, "01_TimeSeries_Temperature.tiff"), plot = p_temp, width = 11, height = 5.5, dpi = 300)
    ggsave(file.path(out_dir, "02_TimeSeries_Abundance.tiff"), plot = p_abund, width = 11, height = 6, dpi = 300)
  }
}

#################################################
# Section 5. Monthly Interannual Temperature Trends (Nested Strip 적용)
#################################################
write_log("월별 수온 변화 트렌드 분석 및 시각화 진행 중 (2단 구획 분할 모드)...")

df_temp <- data.frame(sample_data(ps_raw_updated)) %>% 
  select(Date, Temperature) %>% 
  na.omit()

df_temp$Year <- as.numeric(format(df_temp$Date, "%Y"))
df_temp$Month <- factor(format(df_temp$Date, "%m"), 
                        levels = sprintf("%02d", 1:12), 
                        labels = month.name) 

# 국제 표준 별표 생성 함수 (ns는 생략하여 깔끔함 유지)
get_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("") 
}

# APA 스타일 p-value 변환 함수
format_apa_pval <- function(p) {
  if (p < 0.001) return("< .001")
  if (p < 0.01)  return("< .01")
  if (p < 0.05)  return("< .05")
  p_str <- sprintf("%.3f", p)
  return(paste0("= ", sub("^0", "", p_str)))
}

# 월별 선형 회귀 통계량 계산 및 통계 텍스트 생성
trend_stats <- df_temp %>%
  group_by(Month) %>%
  summarise(
    slope = coef(lm(Temperature ~ Year))[2],
    p_val = summary(lm(Temperature ~ Year))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    stars = get_stars(p_val),
    p_text = format_apa_pval(p_val),
    # R plotmath 문법: 통계 텍스트 렌더링용 (이탤릭 p 적용, Slope에 별표 결합)
    stat_text = paste0("plain('Slope: ", sprintf("%+.3f", slope), stars, ", ') * italic(p) ~ '", p_text, "'")
  ) %>%
  ungroup()

# stat_text를 Factor로 변환하여 Month와 동일한 순서(1~12월) 보장
trend_stats$stat_text <- factor(trend_stats$stat_text, levels = trend_stats$stat_text)

# 원본 데이터에 통계 텍스트 컬럼 매핑
df_temp <- df_temp %>%
  left_join(trend_stats %>% select(Month, stat_text), by = "Month")

shading_ranges_num <- data.frame(
  xmin = years_seq - 0.5, xmax = years_seq + 0.5, ymin = -Inf, ymax = Inf, year = years_seq
) %>% filter(year %% 2 == 0)

monthly_caption <- "Plot Elements & Statistics:\n- Grey shaded backgrounds indicate even years.\n- Solid GREEN line (#347433) represents the linear regression trend for each month.\n- Significance thresholds: * p < .05, ** p < .01, *** p < .001 (ns is completely omitted for a clean look).\n- Statistics are cleanly separated into a nested 2-layer strip above each data plotting area."

write_log("\n==================================================")
write_log("[LOGGED PLOT NOTES: MONTHLY TRENDS]")
write_log("--- Plot: Monthly Interannual Temperature Trends (3x4 Grid) ---")
for (line in strsplit(monthly_caption, "\n")[[1]]) if(trimws(line) != "") write_log(line)
write_log("==================================================\n")

# 월별 수온 변화 플롯 생성 (ggh4x의 Nested Strip으로 2단 분리 구획 적용)
p_monthly_temp <- ggplot(df_temp, aes(x = Year, y = Temperature)) +
  geom_rect(data = shading_ranges_num,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.4) +
  geom_point(alpha = 0.5, color = "red", size = 2) +
  geom_smooth(method = "lm", color = "#347433", fill = "#347433", alpha = 0.15, linetype = "solid", se = TRUE) +
  
  # [핵심] ggh4x::facet_wrap2 적용: 2개의 변수를 묶어 위아래 2개의 박스로 완벽히 분리 생성
  facet_wrap2(vars(Month, stat_text), nrow = 3, ncol = 4, 
              labeller = labeller(stat_text = label_parsed), # 통계 텍스트만 수식 변환
              strip = strip_nested(
                # 상/하단 박스의 배경과 테두리 지정 (서로 맞닿는 부분에 가로선이 형성됨)
                background_x = list(
                  element_rect(fill = "grey90", color = "black", linewidth = 0.5), # 1단 (월)
                  element_rect(fill = "grey95", color = "black", linewidth = 0.5)  # 2단 (통계) - 살짝 더 밝은 회색
                ),
                # 상/하단 텍스트 스타일 각각 지정
                text_x = list(
                  element_text(face = "bold", size = 12, margin = margin(t=4, b=4)),
                  element_text(face = "plain", size = 11, margin = margin(t=4, b=4))
                ),
                by_layer_x = TRUE
              )) +
  
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Monthly Interannual Temperature Trends",
    x = "Sample Date",
    y = "Temperature (°C)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank()
  )

print(p_monthly_temp)

if(save_output) {
  ggsave(file.path(out_dir, "03_Monthly_Temperature_Trends.tiff"), plot = p_monthly_temp, width = 12, height = 9, dpi = 300)
}
write_log("월별 수온 변화 플롯(2단 분할 스트립 적용) 3x4 그리드 출력 완료.")

#################################################
# Section 6. Final Results & Session Info Logging
#################################################
if(save_output) {
  write_log("분석 환경(패키지 및 세션 정보)을 저장합니다.")
  sink(file.path(out_dir, "session_info_phase2.txt"))
  print(sessionInfo())
  sink()
}

write_log("=== Phase 2 분석 정상 종료 ===")

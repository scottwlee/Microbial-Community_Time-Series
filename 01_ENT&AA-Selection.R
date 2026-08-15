################################################################################
# 마이크로바이옴 데이터 전처리 및 엔트로피 분석 파이프라인 [Phase_1]
################################################################################

#################################################
# USER SETTINGS (사용자 설정 영역)
#################################################
# 1. 입출력 경로 설정 (★ 절대 경로 확인 필수)
input_dir      <- "/home/scott/Data/EDM_16SV4_PA"
base_out_dir   <- "/home/scott/EDM_16SV4_PA"

# 2. 파일명 설정
file_otu       <- "phy_asv-table.csv"
file_tax       <- "phy_taxonomy.csv"
file_meta      <- "phy_metadata.csv"
file_fasta     <- "sequences.fasta" # FASTA 파일 추가

# 3. 분석 파라미터 설정 (Cut-off 스위치 및 임계값)
target_taxa      <- "Cyanobacteria"

# [핵심 옵션] 자동 vs 수동 임계값 스위치
# TRUE (기본값): 데이터 기반 수학적 변곡점(Elbow point) 알고리즘이 알아서 최적의 임계값을 찾아냅니다. (아래 수동 값 무시됨)
# FALSE: 알고리즘을 사용하지 않고, 연구자가 아래에 입력한 'quant_ent'와 'quant_abund' 값을 기준으로 데이터를 정확히 잘라냅니다.
use_auto_cutoff  <- TRUE  

# 수동 설정 시 적용될 분위수 (use_auto_cutoff가 FALSE일 때만 작동합니다)
# 예: 0.93은 상위 7%, 0.90은 상위 10%를 의미합니다.
quant_ent        <- 0.93  # 엔트로피 필터링 분위수
quant_abund      <- 0.90  # 절대 풍부도 필터링 분위수

# 4. 시스템 설정
n_threads      <- min(33, parallel::detectCores())
save_plot      <- TRUE

#################################################
# Section 1. Load Packages
#################################################
required_packages <- c("infotheo", "tibble", "tidyverse", "dplyr", 
                       "readxl", "ggplot2", "phyloseq", "parallel", "pbapply", "mgcv")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# FASTA 처리를 위한 Biostrings 패키지 로드 (Bioconductor 의존성)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if(!require("Biostrings", character.only = TRUE)) {
  BiocManager::install("Biostrings", ask = FALSE)
  library("Biostrings", character.only = TRUE)
}

#################################################
# Section 2. Initialization & Logging Setup
#################################################
# Output 디렉토리 단일 생성 (모든 결과물이 한 곳에 저장됨)
out_dir <- file.path(base_out_dir, "01_Phase1_Output")

if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 로그 파일 생성 (생성된 단일 폴더 내 저장)
log_file <- file.path(out_dir, paste0("phase1_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
write_log <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"), file = log_file, append = TRUE)
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

write_log("=== 마이크로바이옴 파이프라인 분석 시작 ===")
write_log(paste("Input 경로:", input_dir))
write_log(paste("Output 경로:", out_dir))
write_log(sprintf("임계값(Cut-off) 적용 모드: %s", ifelse(use_auto_cutoff, "AUTO (수학적 변곡점)", "MANUAL (사용자 지정 분위수)")))

#################################################
# Section 3. Load Data & Create Raw Phyloseq
#################################################
path_otu   <- file.path(input_dir, file_otu)
path_tax   <- file.path(input_dir, file_tax)
path_meta  <- file.path(input_dir, file_meta)
path_fasta <- file.path(input_dir, file_fasta)

if(!file.exists(path_otu)) stop("Error: OTU file not found. 경로를 확인하세요.")
if(!file.exists(path_fasta)) stop("Error: FASTA file not found. 경로를 확인하세요.")

# [데이터 무결성] 원본 데이터(raw)는 읽기 전용으로 유지
otu_raw  <- read.csv(path_otu)
tax_raw  <- read.csv(path_tax)
meta_raw <- read.csv(path_meta)

# 행렬 변환
otu_mat <- otu_raw %>% tibble::column_to_rownames("ID") %>% as.matrix()
tax_mat <- tax_raw %>% tibble::column_to_rownames("ID") %>% as.matrix()
samples_df <- meta_raw %>% tibble::column_to_rownames("Sample_ID")

# FASTA 파일 읽기
fasta_seqs <- Biostrings::readDNAStringSet(path_fasta)

# Taxonomy ID와 FASTA ID 일치 여부 확인 및 콘솔 출력
total_tax_ids <- nrow(tax_mat)
matched_ids   <- sum(rownames(tax_mat) %in% names(fasta_seqs))

cat("\n==================================================\n")
cat("[FASTA Sequence Matching Result]\n")
cat(sprintf("Total Taxonomy IDs: %d\n", total_tax_ids))
cat(sprintf("Matched FASTA IDs: %d\n", matched_ids))
cat("==================================================\n\n")
write_log(sprintf("FASTA 서열 매칭 완료: 총 %d개 Taxonomy 중 %d개 ID 일치", total_tax_ids, matched_ids))

# Phyloseq 객체 구성요소 생성
OTU    <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX    <- tax_table(tax_mat)
samples<- sample_data(samples_df)
REFSEQ <- refseq(fasta_seqs) # FASTA 서열 정보 추가

# 첫 번째 Phyloseq 객체 생성 (REFSEQ 포함)
Bac_Phyloseq_Raw <- phyloseq(OTU, TAX, samples, REFSEQ)

# [출력 1] 오리지널 Phyloseq 객체 별도 저장
saveRDS(Bac_Phyloseq_Raw, file = file.path(out_dir, "Bac_Phyloseq_Raw.rds"))
write_log("Original Phyloseq 객체 (FASTA 포함) 저장 완료.")

#################################################
# Section 4. Data Processing (Normalization & Scaling)
#################################################
normalize_and_scale <- function(physeq_obj, abundance_vector) {
  otu_df <- as.data.frame(otu_table(physeq_obj))
  tax_df <- as.data.frame(tax_table(physeq_obj))
  common_names <- intersect(rownames(otu_df), rownames(tax_df))
  otu_matched <- otu_df[common_names, ]
  rownames(otu_matched) <- tax_df[common_names, "ASV"]
  otu_norm <- apply(otu_matched, 2, function(x) x / sum(x, na.rm = TRUE))
  otu_scaled <- sweep(as.data.frame(otu_norm), 2, abundance_vector, `*`)
  return(otu_scaled)
}

taxa_syn <- taxa_names(Bac_Phyloseq_Raw)[tax_table(Bac_Phyloseq_Raw)[, "L2"] == target_taxa]
Syn_phyloseq <- prune_taxa(taxa_syn, Bac_Phyloseq_Raw)
otu_table_Syn_scaled <- normalize_and_scale(Syn_phyloseq, samples_df$Syn_abundance)

taxa_hb <- taxa_names(Bac_Phyloseq_Raw)[tax_table(Bac_Phyloseq_Raw)[, "L2"] != target_taxa]
HB_phyloseq <- prune_taxa(taxa_hb, Bac_Phyloseq_Raw)
otu_table_HB_scaled <- normalize_and_scale(HB_phyloseq, samples_df$HB_abundance)

#################################################
# Section 5. Entropy Calculation
#################################################
ASV_table_Total_AA <- rbind(otu_table_Syn_scaled, otu_table_HB_scaled)
ASV_table_Total_AA <- tibble::rownames_to_column(ASV_table_Total_AA, var = "ASV")
otu_mat_ASV <- ASV_table_Total_AA %>% tibble::column_to_rownames("ASV") %>% as.matrix()
tax_mat_ASV <- tax_raw %>% tibble::column_to_rownames("ASV") %>% as.matrix()

OTU_AA <- otu_table(otu_mat_ASV, taxa_are_rows = TRUE)
TAX_AA <- tax_table(tax_mat_ASV)
ASV_Total_AA_phyloseq <- phyloseq(OTU_AA, TAX_AA, samples)

Bac_Phyloseq1 <- prune_taxa(taxa_sums(ASV_Total_AA_phyloseq) > 0, ASV_Total_AA_phyloseq)
otu_table_df_for_ent <- as.data.frame(otu_table(Bac_Phyloseq1))

entropy_ps <- function(ts_object) {
  n_bin_for_ts <- ceiling(sqrt(length(ts_object)))
  binned_ts <- infotheo::discretize(c(ts_object), disc = "equalwidth", nbins = n_bin_for_ts)
  infotheo::entropy(binned_ts)
}

write_log("멀티코어를 활용한 엔트로피 계산 중 (Progress Bar 확인)...")
cl <- makeCluster(n_threads)
bac_ent <- pbapply(otu_table_df_for_ent, 1, entropy_ps, cl = cl)
stopCluster(cl)

otu_absolute_abund <- data.frame(otu_table(Bac_Phyloseq1))
taxa_Mean_AA_abund <- rowMeans(otu_absolute_abund)

#################################################
# Section 5-1. Scientific Threshold Calculation (Auto vs Manual)
#################################################
# 변곡점(Elbow point)을 찾는 함수 정의 (자동 모드에서 사용됨)
find_elbow_point <- function(x, y) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  p1 <- c(x_norm[1], y_norm[1])
  p2 <- c(x_norm[length(x_norm)], y_norm[length(y_norm)])
  a <- p1[2] - p2[2]
  b <- p2[1] - p1[1]
  c <- p1[1] * p2[2] - p2[1] * p1[2]
  distances <- abs(a * x_norm + b * y_norm + c) / sqrt(a^2 + b^2)
  return(which.max(distances))
}

if (use_auto_cutoff) {
  write_log("데이터 기반의 과학적 변곡점(Elbow point) 자동 탐색을 시작합니다.")
  
  # 1) Threshold 계산 (Auto 모드)
  sorted_ent <- sort(bac_ent, decreasing = FALSE)
  x_ent <- seq_along(sorted_ent)
  elbow_idx_ent <- find_elbow_point(x_ent, sorted_ent)
  optimal_ent_val <- sorted_ent[elbow_idx_ent]
  optimal_quant_ent <- ecdf(bac_ent)(optimal_ent_val)
  
  sorted_abund <- sort(taxa_Mean_AA_abund, decreasing = TRUE)
  x_abund <- seq_along(sorted_abund)
  elbow_idx_abund <- find_elbow_point(x_abund, log10(sorted_abund + 1e-9)) 
  optimal_abund_val <- sorted_abund[elbow_idx_abund]
  optimal_quant_abund <- ecdf(taxa_Mean_AA_abund)(optimal_abund_val)
  
  cat("\n==================================================\n")
  cat("[Data-Driven Auto Thresholds Identified]\n")
  cat(sprintf("Optimal Entropy Quantile: %.4f (Value: %.4f)\n", optimal_quant_ent, optimal_ent_val))
  cat(sprintf("Optimal Abundance Quantile: %.4f (Value: %.4f)\n", optimal_quant_abund, optimal_abund_val))
  cat("==================================================\n\n")
  
  # 동적 Threshold 적용
  q_ent_cutoff <- optimal_ent_val
  q_abund_cutoff <- optimal_abund_val
  quant_ent <- optimal_quant_ent
  quant_abund <- optimal_quant_abund
  
} else {
  write_log("사용자 지정(Manual) 분위수 기반 수동 임계값 적용을 시작합니다.")
  
  # 1) Threshold 계산 (Manual 모드: 사용자가 입력한 quant_ent, quant_abund을 기반으로 실제 Cut-off 값 역산출)
  q_ent_cutoff <- as.numeric(quantile(bac_ent, probs = quant_ent, names = FALSE))
  q_abund_cutoff <- as.numeric(quantile(taxa_Mean_AA_abund, probs = quant_abund, names = FALSE))
  
  # 시각화를 위한 elbow_idx 위치 강제 매칭 (플롯의 수직/수평선 표시용)
  sorted_abund <- sort(taxa_Mean_AA_abund, decreasing = TRUE)
  elbow_idx_abund <- which.min(abs(sorted_abund - q_abund_cutoff))
  
  cat("\n==================================================\n")
  cat("[User-Defined Manual Thresholds Applied]\n")
  cat(sprintf("Manual Entropy Quantile: %.4f (Value: %.4f)\n", quant_ent, q_ent_cutoff))
  cat(sprintf("Manual Abundance Quantile: %.4f (Value: %.4f)\n", quant_abund, q_abund_cutoff))
  cat("==================================================\n\n")
}

# [로그 개선] Plot Notes를 로그 파일에 기록 (모드에 따라 동적 변환)
write_log("-------------------------------------------------------------------")
write_log(sprintf("[LOGGED PLOT NOTES: THRESHOLD DETECTION - %s MODE]", ifelse(use_auto_cutoff, "AUTO", "MANUAL")))
if (use_auto_cutoff) {
  write_log(" -> Plot 1 (01_Entropy_Density_Cutoff.tiff): The cut-off threshold was mathematically determined by finding the 'elbow point' (maximum perpendicular distance to the secant line) of the sorted entropy curve.")
  write_log(" -> Plot 2 (02_Rank_Abundance_Cutoff.tiff): The cut-off threshold was mathematically determined by finding the 'elbow point' (maximum perpendicular distance to the secant line) on the log-transformed abundance curve.")
  write_log(" -> Plot 2 (02_Rank_Abundance_Cutoff.tiff): The blue dotted line represents the mathematical elbow point threshold on the log10 scale.")
} else {
  write_log(sprintf(" -> Plot 1 (01_Entropy_Density_Cutoff.tiff): The cut-off threshold was manually set to the %.1fth percentile.", quant_ent * 100))
  write_log(sprintf(" -> Plot 2 (02_Rank_Abundance_Cutoff.tiff): The cut-off threshold was manually set to the %.1fth percentile.", quant_abund * 100))
  write_log(" -> Plot 2 (02_Rank_Abundance_Cutoff.tiff): The blue dotted line represents the manually set threshold on the log10 scale.")
}
write_log("-------------------------------------------------------------------")

# 2) 예비 분석(Preliminary) 플롯 생성 및 출력
prelim_df <- data.frame(ASV = names(bac_ent), Entropy = bac_ent, Abundance = taxa_Mean_AA_abund)

p_density_ent <- ggplot(prelim_df, aes(x = Entropy)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = q_ent_cutoff, color = "red", linetype = "dashed", linewidth = 1) +
  # [수정됨] 한 줄 표기(줄바꿈 제거), 볼드 해제, 소수점 3자리(%.3f)
  annotate("text", x = q_ent_cutoff, y = Inf, 
           label = sprintf(" Cut-off: %.3f (%.1f%%)", q_ent_cutoff, quant_ent * 100), 
           color = "red", hjust = -0.1, vjust = 2) +
  theme_bw(base_size = 14) + 
  # [수정됨] 긴 caption 제거, Subtitle 동적 변환
  labs(
    title = "Density Plot of ASV Shannon Entropy", 
    subtitle = ifelse(use_auto_cutoff, "Data-driven elbow point cut-off applied", "User-defined manual cut-off applied")
  )

p_rank_abund <- prelim_df %>%
  arrange(desc(Abundance)) %>%
  mutate(Rank = row_number()) %>%
  ggplot(aes(x = Rank, y = log10(Abundance + 1e-9))) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = elbow_idx_abund, color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = log10(q_abund_cutoff + 1e-9), color = "blue", linetype = "dotted", linewidth = 0.8) +
  # [수정됨] Cut-off Value -> Cut-off로 변경, 한 줄 표기, 볼드 해제, 소수점 3자리(%.3f)
  annotate("text", x = elbow_idx_abund, y = Inf, 
           label = sprintf(" Cut-off: %.3f (%.1f%%)", q_abund_cutoff, quant_abund * 100), 
           color = "red", hjust = -0.1, vjust = 2) +
  theme_bw(base_size = 14) + 
  # [수정됨] 긴 caption 제거, Subtitle 동적 변환
  labs(
    title = "Rank-Abundance Curve", 
    subtitle = ifelse(use_auto_cutoff, "Mathematical elbow point detection for cut-off", "User-defined manual cut-off applied")
  )

# 플롯 즉시 출력 보장
print(p_density_ent)
print(p_rank_abund)

if(save_plot) {
  ggsave(file.path(out_dir, "01_Entropy_Density_Cutoff.tiff"), plot = p_density_ent, width = 7, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "02_Rank_Abundance_Cutoff.tiff"), plot = p_rank_abund, width = 7, height = 5, dpi = 300)
}

#################################################
# Section 6. Final Filtering & Publication Plotting
#################################################
plot_df <- data.frame(
  ASV = rownames(otu_absolute_abund),
  Abundance = taxa_Mean_AA_abund,
  Sqrt_Abundance = sqrt(taxa_Mean_AA_abund),
  Entropy = bac_ent
) %>%
  mutate(Category = ifelse(Entropy > q_ent_cutoff & Abundance > q_abund_cutoff, 
                           "Core/Active ASVs (Filtered)", "Low Ent/Abund (Discarded)"))

p_pub <- ggplot(plot_df, aes(x = Sqrt_Abundance, y = Entropy)) +
  geom_point(aes(fill = Category, color = Category), shape = 21, size = 2.5, alpha = 0.7, stroke = 0.3) +
  scale_fill_manual(values = c("Core/Active ASVs (Filtered)" = "#D55E00", "Low Ent/Abund (Discarded)" = "#E6E6E6")) +
  scale_color_manual(values = c("Core/Active ASVs (Filtered)" = "#A64B00", "Low Ent/Abund (Discarded)" = "#B3B3B3")) +
  geom_vline(xintercept = sqrt(q_abund_cutoff), linetype = "dashed", color = "#0072B2", linewidth = 0.7) +
  geom_hline(yintercept = q_ent_cutoff, linetype = "dashed", color = "#0072B2", linewidth = 0.7) +
  labs(
    title = "Identification of Core ASVs based on Entropy and Abundance",
    x = "Square Root of Absolute Abundance", 
    y = "Shannon Entropy",                   
    # [수정됨] 캡션 내 동적 텍스트 적용 (Auto/Manual 문구 분기)
    caption = sprintf(
      "Note:\n1. Thresholds: %.1fth percentile for Entropy (%.2f), %.1fth percentile for Abundance (%.2f).\n2. Retained ASVs: %d out of %d total ASVs.\n3. Cut-off values were determined by %s.",
      quant_ent * 100, q_ent_cutoff, quant_abund * 100, q_abund_cutoff,
      sum(plot_df$Category == "Core/Active ASVs (Filtered)"), nrow(plot_df),
      ifelse(use_auto_cutoff, "detecting the mathematical elbow points of the data distributions", "user-defined manual percentiles")
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom", 
    legend.title = element_blank(),
    plot.caption = element_text(hjust = 0, face = "italic", color = "grey30", size = 10, margin = margin(t = 15)),
    aspect.ratio = 1/1.2 
  )

# 최종 플롯 즉시 출력 보장
print(p_pub)

if(save_plot) {
  ggsave(file.path(out_dir, "AA_Shannon_Ent_ASV.tiff"), 
         plot = p_pub, width = 8, height = 7, dpi = 300)
}

#################################################
# Section 7. Create Filtered Phyloseq & Save Results
#################################################
# 필터링 통과 대상 ASV 이름 분리
target_asvs <- plot_df$ASV[plot_df$Category == "Core/Active ASVs (Filtered)"]
high_entropy_abundance_ASVs <- otu_absolute_abund[target_asvs, ]

# CSV 추출
write.csv(high_entropy_abundance_ASVs, file.path(out_dir, "ASV_table_high_ent_AA.csv"), row.names = TRUE)

# 원본 데이터에 맞게 ASV 이름을 다시 기존 ID로 매핑(Mapping)
tax_df_raw <- as.data.frame(tax_table(Bac_Phyloseq_Raw))
original_ids <- rownames(tax_df_raw)[tax_df_raw$ASV %in% target_asvs]

# 방어적 프로그래밍: 필터링 결과가 비어있는지 확인
if(length(original_ids) == 0) {
  stop("Error: No ASVs passed the filtering thresholds. Please check the data or cut-off values.")
}

# 매핑된 원본 ID를 사용하여 두 번째 객체 가지치기(Pruning)
Filtered_Phyloseq <- prune_taxa(original_ids, Bac_Phyloseq_Raw)

# [출력 2] 필터링이 완료된 Phyloseq 객체 별도 저장
saveRDS(Filtered_Phyloseq, file = file.path(out_dir, "Filtered_Phyloseq.rds"))
write_log(paste("Filtered Phyloseq 객체 생성 완료 (통과된 핵심 ASV 수:", ntaxa(Filtered_Phyloseq), ")"))

# [재현성 확보] 분석 환경(세션 정보) 로깅 추가
write_log("분석 환경 세션 정보를 저장합니다.")
sink(file.path(out_dir, "session_info_phase1.txt"))
print(sessionInfo())
sink()

write_log("=== Phase 1 분석 정상 종료 ===")

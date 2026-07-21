################################################################################
# Helper Function for Phase 3 - Part 3: MDR S-map (Multiview Distance Regularized S-map)
# File Name: F2_HelperFunction_MDR_Block.R
# Purpose: Phase 3의 Part 1(Simplex)과 Part 2(UIC)로써 얻어진 인과성 결과
#          (Optimal E, Max TE Time-lag)를 바탕으로 MDR S-map 입력용 
#          '다변량 최적 시차 블록(Multiview Block)'을 생성하는 핵심 로직.
# Note: rEDM 0.7.3 레거시 버전에 맞춰 소문자 make_block을 사용합니다.
################################################################################

library(dplyr)
library(rEDM)

make_block_mvd_TOed <- function (block,
                                 uic_res,
                                 effect_var,
                                 E_effect_var,
                                 cause_var_colname = "cause_var",
                                 include_var = "strongest_only",
                                 p_threshold = 0.050,
                                 sort_tp = TRUE,
                                 silent = FALSE) {
  
  # 1. 입력 데이터 및 컬럼 검증
  x_names <- colnames(block)
  
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]
  if (!(effect_var %in% x_names)) stop("Error: 'effect_var'가 데이터 블록에 존재하지 않습니다.")
  if (!(cause_var_colname %in% colnames(uic_res))) stop("Error: 'cause_var_colname'이 uic_res 데이터에 없습니다.")
  if (!("tp" %in% colnames(uic_res))) stop("Error: 'tp' 컬럼이 uic_res 데이터에 필요합니다.")
  if (!("te" %in% colnames(uic_res))) stop("Error: 'te' 컬럼이 uic_res 데이터에 필요합니다.")
  if (!("pval" %in% colnames(uic_res))) stop("Error: 'pval' 컬럼이 uic_res 데이터에 필요합니다.")
  if (!is.data.frame(block)) stop("Error: 입력 'block'은 반드시 data.frame 형식이어야 합니다.")
  
  # 2. 타겟 변수(미생물)의 내재적 동역학(Optimal E)을 바탕으로 기본 상태 공간(State space) 생성
  if (include_var == "tp0_only") {
    block_mvd <- data.frame(block[,effect_var])
    colnames(block_mvd) <- sprintf("%s_tp0", effect_var)
  } else {
    # [핵심 수정] rEDM 0.7.3 버전에 맞는 소문자 make_block 사용
    block_mvd <- data.frame(rEDM::make_block(block[,effect_var], max_lag = E_effect_var)[,-1])
    colnames(block_mvd) <- sprintf("%s_tp%s", effect_var, 0:(-(E_effect_var-1)))
  }
  
  # 3. 유의미한 인과성 데이터(Phase 3 - Part 2 결과) 필터링
  if (!silent) message(sprintf("Pre-screening: pval <= %.3f 및 tp <= 0 (과거 시점) 데이터만 유지합니다.", p_threshold))
  uic_res <- uic_res[uic_res$pval <= p_threshold & uic_res$tp <= 0, ]
  
  if (nrow(uic_res) < 1) stop("Error: 통계적으로 유의미한 원인 변수(온도)가 발견되지 않았습니다.")
  
  # 4. 시차(tp) 기준 정렬 로직
  if (sort_tp) {
    sort_id <- 0
    for (col_i in unique(uic_res[,cause_var_colname])) {
      sort_id_i <- order(-uic_res[uic_res[,cause_var_colname] == col_i, "tp", drop = TRUE]) 
      sort_id <- c(sort_id, max(sort_id) + sort_id_i)
    }
    sort_id <- sort_id[-1]
    uic_res <- uic_res[sort_id,]
  }
  
  # 5. 환경 변수(온도)의 인과성 시차 데이터를 원래 블록에 덧붙이기 (Add-on)
  if (include_var == "all_significant") {
    for (i in 1:nrow(uic_res)) {
      block_new <- dplyr::lag(block[, uic_res[i, cause_var_colname]], n = abs(uic_res[i,"tp"]))
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i, cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "strongest_only") {
    for (cause_i in unique(uic_res[,cause_var_colname])) {
      uic_res_tmp <- uic_res[uic_res[,cause_var_colname] == cause_i,]
      if (!exists("uic_res_new")) {
        uic_res_new <- uic_res_tmp[which.max(uic_res_tmp$te),] 
      } else {
        uic_res_new <- rbind(uic_res_new, uic_res_tmp[which.max(uic_res_tmp$te),]) 
      }
    }
    uic_res <- uic_res_new
    
    for (i in 1:nrow(uic_res)) {
      block_new <- dplyr::lag(block[, unlist(uic_res[i,cause_var_colname])], n = abs(unlist(uic_res[i,"tp"]))) 
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i,cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "tp0_only") {
    block_new <- block[, unique(uic_res[,cause_var_colname])]
    block_new <- data.frame(block_new)
    colnames(block_new) <- sprintf("%s_tp0", unique(uic_res[,cause_var_colname]))
    block_mvd <- cbind(block_mvd, block_new)
  }
  
  return(block_mvd)
}

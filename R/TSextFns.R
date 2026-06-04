#####################################################################
# Helper functions for the TimeSignature gene-selection extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# These wrap and re-use the original TimeStampFns.R pipeline so that
# the extension experiments (i) train an elastic-net TimeSignature on
# an arbitrary *restricted* gene pool, (ii) apply it with the paper's
# two-sample antipodal calibration, and (iii) score it with the three
# clinically-relevant metrics (MAE, nAUC, % within 2 h).
#
# Requires TimeStampFns.R to be sourced first (uses recalibrateExprs,
# trainTimeStamp, predTimeStamp, time2XY, XY2dectime, timeErr).
######################################################################


#---------------------------------------------------------------------
# Two-sample ("antipodal") calibration, factored out of TSexample.R.
# For every sample we find, within the same subject, the timepoint
# closest to 12 h away and centre the sample on the mean of the pair.
# This mimics the clinical setting of just two blood draws ~12 h apart.
#
#   exprMat : genes x samples expression matrix
#   meta    : data.frame with rownames == sample IDs and columns
#             `samp`, `sID` (subject) and a time column
#   timeCol : name of the decimal-24h time column (default "LocalTime")
# Returns a genes x samples matrix of within-pair-centred values, with
# any NA (subject with a single draw) set to 0 (gene assumed flat).
#---------------------------------------------------------------------
twoPointCalibrate <- function(exprMat, meta, timeCol = "LocalTime") {

	timesBySubj <- split(meta[, c("samp", timeCol)], meta$sID)

	findAntipode <- function(tbs) {
		do.call(rbind, lapply(rownames(tbs), function(s) {
			thisTime  <- tbs[s, timeCol]
			avail     <- tbs[[timeCol]]; names(avail) <- rownames(tbs)
			avail[s]  <- NA                      # never pair with self
			antip     <- rownames(tbs)[which.min(timeErr(avail, thisTime + 12))]
			data.frame(samp = s, asamp = antip, stringsAsFactors = FALSE)
		}))
	}

	anti <- do.call(rbind, lapply(timesBySubj, findAntipode))
	rownames(anti) <- anti$samp
	anti <- anti[rownames(meta), ]
	partner <- anti$asamp; names(partner) <- anti$samp

	calib <- sapply(colnames(exprMat), function(s) {
		ctr <- rowMeans(exprMat[, c(s, partner[s]), drop = FALSE], na.rm = TRUE)
		exprMat[, s] - ctr
	})
	calib[is.na(calib)] <- 0
	calib
}


#---------------------------------------------------------------------
# Train an elastic-net TimeSignature on a restricted gene pool.
#
#   trainWSN : genes x samples, within-subject-normalised training data
#   subjIDs  : subject ID per training sample
#   times    : decimal-24h time per training sample
#   genes    : character vector of genes to restrict to (NULL = all)
#   alpha    : elastic-net mixing parameter (0.5, as in the paper)
#   s        : glmnet penalty; NULL -> CV-selected lambda.min
#   seed     : RNG seed for reproducible CV folds
# Returns the timestamp object from trainTimeStamp(), with $genePool
# and $nSelected (number of non-zero predictors) attached.
#---------------------------------------------------------------------
trainTSpool <- function(trainWSN, subjIDs, times, genes = NULL,
                        alpha = 0.5, s = NULL, seed = 194) {
	if (!is.null(genes)) {
		genes    <- intersect(genes, rownames(trainWSN))
		trainWSN <- trainWSN[genes, , drop = FALSE]
	}
	set.seed(seed)
	ts <- trainTimeStamp(
		expr = trainWSN, subjIDs = subjIDs, times = times,
		trainFrac = 1, recalib = FALSE, a = alpha, s = s, plot = FALSE
	)
	ts$genePool  <- rownames(trainWSN)
	ts$nSelected <- if (is.matrix(ts$coef)) nrow(ts$coef) else length(ts$coef)
	ts
}


#---------------------------------------------------------------------
# Three headline accuracy metrics for circadian-time prediction.
# nAUC reproduces the normalised area-under-the-CDF used in Fig 1.
#---------------------------------------------------------------------
evalMetrics <- function(trueTime, predTime) {
	err <- timeErr(trueTime, predTime)
	err <- err[!is.na(err)]
	hrsoff   <- seq(0, 12, length = 49)
	pctOK    <- sapply(hrsoff, function(h) mean(err <= h))   # CDF of |error|
	nAUC     <- sum(pctOK[-1] * diff(hrsoff / 12))
	list(n      = length(err),
	     MAE    = mean(err),
	     medErr = median(err),
	     nAUC   = nAUC,
	     pct2h  = 100 * mean(err <= 2))
}


#---------------------------------------------------------------------
# Apply a trained pool model across studies and tabulate the metrics.
#
#   ts       : timestamp object from trainTSpool()
#   calibMat : genes x samples two-point-calibrated matrix (all samples)
#   meta     : metadata with `study` and the truth time column
#   studies  : which studies to score (default the three validation sets)
#   truthCol : ground-truth time column (default "LocalTime")
#   s        : glmnet penalty used for prediction
# Returns a data.frame, one row per study (+ a pooled "Valid.all" row).
#---------------------------------------------------------------------
evalAcrossStudies <- function(ts, calibMat, meta,
                              studies = c("V1", "V2", "V3"),
                              truthCol = "LocalTime",
                              s = ts$cv.fit$lambda.min) {

	x   <- calibMat[ts$genePool, , drop = FALSE]
	prd <- predTimeStamp(ts, newx = x, s = s)
	names(prd) <- colnames(calibMat)

	rows <- lapply(studies, function(st) {
		idx <- rownames(meta)[meta$study == st]
		m   <- evalMetrics(meta[idx, truthCol], prd[idx])
		data.frame(study = st, n = m$n, nGene = length(ts$genePool),
		           nSel = ts$nSelected, MAE = m$MAE, medErr = m$medErr,
		           nAUC = m$nAUC, pct2h = m$pct2h, stringsAsFactors = FALSE)
	})
	# pooled across the validation studies
	idxAll <- rownames(meta)[meta$study %in% studies]
	mAll   <- evalMetrics(meta[idxAll, truthCol], prd[idxAll])
	rows[[length(rows) + 1]] <- data.frame(
		study = "Valid.all", n = mAll$n, nGene = length(ts$genePool),
		nSel = ts$nSelected, MAE = mAll$MAE, medErr = mAll$medErr,
		nAUC = mAll$nAUC, pct2h = mAll$pct2h, stringsAsFactors = FALSE)

	do.call(rbind, rows)
}


library(ggplot2)

#---------------------------------------------------------------------
# Gene-count vs accuracy trade-off plot (ggplot2 Publication-Ready ver)
#   df      : data.frame with columns nGene, nAUC, MAE (pooled rows)
#   target  : clinical nAUC floor to draw (default 0.80)
#---------------------------------------------------------------------
tradeoffPlot <- function(df, target = 0.80, maeTarget = 2.5,
                         main = "Gene count vs. accuracy trade-off") {
  
  df <- df[order(df$nGene), ]
  
  # 1. 이중 축(Dual-axis)을 위한 동적 스케일링 설정
  # 왼쪽 y축(nAUC) 범위 설정 (Combined plot과 동일하게 맞춤)
  p_min <- min(0.75, min(df$nAUC, na.rm = TRUE))
  p_max <- max(0.88, max(df$nAUC, na.rm = TRUE))
  
  # 오른쪽 y축(MAE) 범위 설정
  s_min <- min(1.5, min(df$MAE, na.rm = TRUE))
  s_max <- max(3.0, max(df$MAE, na.rm = TRUE))
  
  # 변환 공식: MAE_scaled = MAE * slope + intercept
  slope <- (p_max - p_min) / (s_max - s_min)
  intercept <- p_min - slope * s_min
  
  # MAE 값을 nAUC 그리는 공간에 맞춰 스케일링
  df$MAE_scaled <- df$MAE * slope + intercept
  
  # 2. ggplot 객체 생성
  p <- ggplot(df, aes(x = nGene)) +
    
    # 가이드라인 (Clinical floor nAUC)
    geom_hline(yintercept = target, linetype = "dashed", color = "firebrick", alpha = 0.7) +
    annotate("text", x = min(df$nGene), y = target - 0.005, 
             label = sprintf("clinical floor nAUC = %.2f", target), 
             color = "firebrick", size = 3.5, hjust = 0, fontface = "italic") +
    
    # MAE (오른쪽 축 대응) - 주황색
    geom_line(aes(y = MAE_scaled, color = "MAE (h)"), alpha = 0.6) +
    geom_point(aes(y = MAE_scaled, color = "MAE (h)", shape = "MAE (h)"), size = 2.5) +
    
    # nAUC (왼쪽 축 대응) - 네이비
    geom_line(aes(y = nAUC, color = "nAUC"), alpha = 0.6) +
    geom_point(aes(y = nAUC, color = "nAUC", shape = "nAUC"), size = 2.5) +
    
    # x축(로그 스케일) 및 이중 y축 설정
    scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000)) +
    scale_y_continuous(
      name = "nAUC",
      limits = c(p_min, p_max),
      # sec_axis를 통해 오른쪽 축의 숫자를 역변환하여 원래 MAE 값으로 표시
      sec.axis = sec_axis(~ (. - intercept) / slope, name = "MAE (h)")
    ) +
    
    # 색상 및 범례 매핑
    scale_color_manual(name = NULL, values = c("nAUC" = "navy", "MAE (h)" = "darkorange")) +
    scale_shape_manual(name = NULL, values = c("nAUC" = 16, "MAE (h)" = 17)) +
    
    # 3. 테마 설정 (격자 제거, 이중 축 색상 지정, 여백 확보)
    labs(
      title = main, 
      x = "Number of genes in candidate pool (log scale)"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.98, 0.10),      # 두 그래프 선을 피해 우측 중앙 쯤에 배치
      legend.justification = c("right", "center"),
      legend.background = element_blank(),  # 범례 배경 투명화
      legend.key = element_blank(),         # 아이콘 배경 투명화
      
      # 텍스트 여백 및 색상 설정 (각 축의 색상을 데이터 선 색상과 일치시킴!)
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 10), color = "navy"),
      axis.text.y.left = element_text(color = "navy"),
      axis.title.y.right = element_text(margin = margin(l = 10), color = "darkorange"),
      axis.text.y.right = element_text(color = "darkorange")
    )
  
  print(p)
  invisible(df)
}
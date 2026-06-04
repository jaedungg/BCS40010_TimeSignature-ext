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


#---------------------------------------------------------------------
# Gene-count vs accuracy trade-off plot (the project deliverable).
#   df      : data.frame with columns nGene, nAUC, MAE (pooled rows)
#   target  : clinical nAUC floor to draw (default 0.80)
#---------------------------------------------------------------------
tradeoffPlot <- function(df, target = 0.80, maeTarget = 2.5,
                         main = "Gene count vs. accuracy trade-off") {
	df <- df[order(df$nGene), ]
	opar <- par(mar = c(4.2, 4.2, 3, 4.2)); on.exit(par(opar))
	plot(df$nGene, df$nAUC, log = "x", type = "b", pch = 19, col = "navy",
	     ylim = c(min(0.6, df$nAUC), max(0.9, df$nAUC)),
	     xlab = "Number of genes in pool (log scale)",
	     ylab = "nAUC", main = main)
	abline(h = target, lty = 2, col = "red")
	text(min(df$nGene), target, sprintf("nAUC = %.2f", target),
	     pos = 3, col = "red", cex = 0.8)
	# overlay MAE on a second axis
	usr <- par("usr"); yr <- range(df$MAE)
	sc  <- function(v) usr[3] + (v - yr[1]) / diff(yr) * diff(usr[3:4])
	lines(df$nGene, sc(df$MAE), type = "b", pch = 17, col = "darkorange")
	axis(4, at = sc(pretty(df$MAE)), labels = pretty(df$MAE), col.axis = "darkorange")
	mtext("MAE (h)", side = 4, line = 2.5, col = "darkorange")
	legend("bottomright", bty = "n", lty = 1, pch = c(19, 17),
	       col = c("navy", "darkorange"), legend = c("nAUC", "MAE (h)"))
	invisible(df)
}

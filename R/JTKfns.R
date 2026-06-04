#####################################################################
# JTK_Cycle-style rhythmicity detection for TimeSignature extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
# Extension: "Gene Selection Optimization via JTK_Cycle Pre-Filtering"
#
# Compact, fully vectorised re-implementation of the core idea of
# JTK_Cycle (Hughes, Hogenesch & Kornacker, J Biol Rhythms 2010):
#
#   * the expression profile of each gene is compared, by Kendall's
#     rank correlation (tau), against a family of cosine reference
#     templates of fixed period (24 h) spanning a grid of phases;
#   * the optimal phase is the one maximising |tau|; the test
#     statistic is Kendall's S at that phase;
#   * a Bonferroni correction is applied over the phases searched
#     (this is JTK_Cycle's ADJ.P), and a Benjamini-Hochberg FDR is
#     applied across genes.
#
# This is the same non-parametric, rank-based rhythmicity filter that
# TimeMachine (Huang & Braun 2024) and tauFisher (Duan & Ngo 2024)
# place upstream of their regression step.  We run it on the
# within-subject-normalised training data so that population-level
# circadian rhythms are detected rather than between-subject baseline
# differences.
######################################################################


#---------------------------------------------------------------------
# Bin samples into circadian-time bins and average each gene.
#   exprMat : genes x samples matrix
#   times   : decimal 24h sampling time for each column
#   binHours: width of each time bin (h); 1 -> 24 bins over the day
# Returns a list with the binned expression (bins x genes) and the
# representative (centre) time of each occupied bin.
#---------------------------------------------------------------------
binByTime <- function(exprMat, times, period = 24, binHours = 1) {
	t <- times %% period
	breaks <- seq(0, period, by = binHours)
	bin <- cut(t, breaks = breaks, include.lowest = TRUE, labels = FALSE)
	occupied <- sort(unique(bin[!is.na(bin)]))
	binTime <- (breaks[-length(breaks)] + breaks[-1]) / 2          # bin centres
	exprBin <- t(sapply(occupied, function(b) {
		rowMeans(exprMat[, which(bin == b), drop = FALSE], na.rm = TRUE)
	}))
	rownames(exprBin) <- paste0("t", round(binTime[occupied], 1))
	colnames(exprBin) <- rownames(exprMat)
	list(expr = exprBin, time = binTime[occupied])
}


#---------------------------------------------------------------------
# JTK_Cycle-style rhythmicity test.
#   exprMat : genes x samples matrix (within-subject normalised data
#             is recommended; see header)
#   times   : decimal 24h sampling time for each column
#   period  : assumed period of the rhythm (h)
#   nPhase  : number of cosine phase templates searched over [0,period)
#   binHours: circadian-time bin width used to build the profiles
#
# Returns a data.frame, one row per gene, sorted by ascending FDR:
#   gene   - gene symbol
#   tau    - Kendall's tau at the optimal phase (signed)
#   phase  - estimated peak phase (h)
#   amp    - peak-to-trough amplitude of the binned profile
#   ADJ.P  - phase-Bonferroni p-value  (== JTK_Cycle's ADJ.P)
#   BH.Q   - Benjamini-Hochberg FDR across genes
#---------------------------------------------------------------------
jtkCycle <- function(exprMat, times, period = 24, nPhase = 24, binHours = 1) {

	# 1. build one circadian profile per gene -------------------------
	b       <- binByTime(exprMat, times, period = period, binHours = binHours)
	exprBin <- b$expr                       # bins x genes
	binTime <- b$time
	nb      <- nrow(exprBin)
	nGene   <- ncol(exprBin)
	stopifnot(nb >= 4)

	# 2. enumerate the C(nb,2) ordered bin pairs ----------------------
	pr     <- combn(nb, 2)
	nPair  <- ncol(pr)
	# sign of (gene value in earlier-listed bin  -  later-listed bin),
	# for every pair x every gene  ->  nPair x nGene
	geneSign <- sign(exprBin[pr[1, ], , drop = FALSE] -
	                 exprBin[pr[2, ], , drop = FALSE])

	# 3. Kendall S of every gene against every cosine phase template ---
	phases <- seq(0, period, length.out = nPhase + 1)[-(nPhase + 1)]
	S <- matrix(0, nGene, nPhase)
	for (k in seq_along(phases)) {
		ref     <- cos(2 * pi * (binTime - phases[k]) / period)
		refSign <- sign(ref[pr[1, ]] - ref[pr[2, ]])
		S[, k]  <- as.vector(crossprod(geneSign, refSign))
	}

	# 4. pick optimal phase, form tau and the JTK p-value -------------
	bestK  <- max.col(abs(S), ties.method = "first")
	Sbest  <- S[cbind(seq_len(nGene), bestK)]
	varS   <- nPair * (2 * nb + 5) / 9          # = nb(nb-1)(2nb+5)/18, with nPair=nb(nb-1)/2
	pPhase <- 2 * pnorm(-abs(Sbest) / sqrt(varS))
	ADJ.P  <- pmin(1, pPhase * nPhase)          # Bonferroni over phases
	tau    <- Sbest / nPair

	# peak phase reported on the conventional "time of peak" convention
	peak   <- phases[bestK]
	peak[Sbest < 0] <- (peak[Sbest < 0] + period / 2) %% period
	amp    <- apply(exprBin, 2, function(z) max(z) - min(z))

	out <- data.frame(
		gene  = colnames(exprBin),
		tau   = tau,
		phase = peak,
		amp   = amp,
		ADJ.P = ADJ.P,
		BH.Q  = p.adjust(ADJ.P, method = "BH"),
		row.names = colnames(exprBin),
		stringsAsFactors = FALSE
	)
	out[order(out$ADJ.P, -abs(out$tau)), ]
}


#---------------------------------------------------------------------
# Convenience: names of the genes passing a JTK threshold.
#   jtk     - data.frame returned by jtkCycle()
#   thresh  - significance cut-off
#   use     - which column to threshold on ("ADJ.P" or "BH.Q")
#---------------------------------------------------------------------
jtkGenes <- function(jtk, thresh = 0.1, use = c("ADJ.P", "BH.Q")) {
	use <- match.arg(use)
	jtk$gene[jtk[[use]] < thresh]
}

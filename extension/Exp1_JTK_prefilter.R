#####################################################################
# Experiment 1 -- JTK_Cycle pre-filter, then elastic net
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Question: does restricting the elastic-net candidate pool to the
#           most rhythmic genes (JTK_Cycle pre-filter) improve gene-
#           selection efficiency while preserving accuracy?
#
# Design:   rank all 7,615 genes by JTK_Cycle rhythmicity on the
#           GSE39445 training data, then retrain the elastic-net
#           TimeSignature on the top-K rhythmic genes for a sweep of
#           pool sizes K.  Each model is validated, with two blood
#           draws per subject, on GSE48113 / GSE56931 / GSE113883.
#           The pool size is the design variable; the deliverable is
#           the gene-count vs accuracy trade-off curve.
#
# Note on thresholds: on this densely-sampled (~213-sample) training
# set the absolute JTK p-values are liberal, so the nominal p<0.1/
# 0.05/0.01 cut-offs select thousands of genes.  Because only the
# *ranking* is needed to define a pre-filtered pool, we drive the
# sweep by pool size K (equivalently, progressively stricter
# rhythmicity ranks) and annotate where the nominal cut-offs land.
#
# Assumes 00_setup.R has been sourced.
#####################################################################

if (!exists("jtkRanked")) source("00_setup.R")

#--- reference: how many genes pass the nominal JTK cut-offs ----------
cat("\n--- JTK_Cycle gene counts at the synopsis thresholds ---\n")
for (p in c(0.1, 0.05, 0.01))
	cat(sprintf("  ADJ.P < %.2f : %5d genes   |   BH.Q < %.2f : %5d genes\n",
	            p, sum(jtk$ADJ.P < p), p, sum(jtk$BH.Q < p)))

#--- pool-size sweep -------------------------------------------------
Kgrid <- c(10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200,
           300, 500, 1000, 2000, length(jtkRanked))

exp1 <- do.call(rbind, lapply(Kgrid, function(K) {
	pool <- head(jtkRanked, K)
	ts   <- trainTSpool(trainWSN, trainSubjs, trainTimes,
	                    genes = pool, alpha = 0.5, s = NULL, seed = 194)
	res  <- evalAcrossStudies(ts, calibAll, all.meta,
	                          studies = c("V1", "V2", "V3"))
	res$K <- K
	cat(sprintf("  K=%5d  selected=%2d  Valid.all: nAUC=%.3f  MAE=%.2f  %%2h=%.0f\n",
	            K, ts$nSelected,
	            res$nAUC[res$study == "Valid.all"],
	            res$MAE[res$study == "Valid.all"],
	            res$pct2h[res$study == "Valid.all"]))
	res
}))

#--- save & summarise ------------------------------------------------
write.csv(exp1, "output/Exp1_JTK_prefilter_results.csv", row.names = FALSE)

exp1.pooled <- exp1[exp1$study == "Valid.all", ]
cat("\n--- Experiment 1: pooled validation metrics by pool size ---\n")
print(exp1.pooled[, c("K", "nGene", "nSel", "MAE", "nAUC", "pct2h")],
      row.names = FALSE, digits = 3)

#--- trade-off curve -------------------------------------------------
pdf("output/Exp1_tradeoff.pdf", width = 7, height = 5)
tradeoffPlot(exp1.pooled, target = 0.80,
             main = "Exp 1: JTK pre-filtered pool size vs accuracy")
dev.off()

# minimum pool that still clears the clinical nAUC >= 0.80 floor
ok <- exp1.pooled[exp1.pooled$nAUC >= 0.80, ]
if (nrow(ok)) {
	minK <- ok[which.min(ok$nGene), ]
	cat(sprintf("\nMinimum JTK-filtered pool with nAUC >= 0.80: K=%d genes (%d selected), nAUC=%.3f, MAE=%.2f h\n",
	            minK$nGene, minK$nSel, minK$nAUC, minK$MAE))
}
cat("Exp 1 done -> output/Exp1_JTK_prefilter_results.csv, output/Exp1_tradeoff.pdf\n")

#####################################################################
# Master driver for the TimeSignature gene-selection extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Runs both experiments and assembles the headline deliverable:
# a single gene-count vs accuracy trade-off curve (Exp 1 JTK pool
# sweep) with the published gene sets (Exp 2) overlaid as points.
#
# Usage (working directory = TimeSignatR/extension/):
#   Rscript run_extension.R
#####################################################################

source("00_setup.R")
source("Exp1_JTK_prefilter.R")
source("Exp2_geneset_substitution.R")

#--- combined trade-off figure ---------------------------------------
e1 <- exp1.pooled[order(exp1.pooled$nGene), ]
e2 <- exp2.pooled[exp2.pooled$set != "fullPool", ]

pdf("output/Fig_tradeoff_combined.pdf", width = 8, height = 5.5)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(e1$nGene, e1$nAUC, log = "x", type = "b", pch = 19, col = "navy",
     ylim = c(0.75, 0.88), xlim = c(3, 8000),
     xlab = "Number of genes in candidate pool (log scale)",
     ylab = "nAUC (pooled validation: GSE48113/56931/113883)",
     main = "TimeSignature gene-count vs accuracy trade-off")
abline(h = 0.80, lty = 2, col = "red")
text(3, 0.805, "clinical floor nAUC = 0.80", pos = 4, col = "red", cex = 0.8)

# overlay the published gene sets (Exp 2)
cols <- c(TScore18 = "darkgreen", TimeMachine37 = "darkorange",
          tauFisher = "purple", clockOnly = "brown", intersection = "magenta")
points(e2$nGene, e2$nAUC, pch = 17, cex = 1.4, col = cols[e2$set])
text(e2$nGene, e2$nAUC, e2$set, pos = c(1, 3, 1, 1, 3), cex = 0.75,
     col = cols[e2$set])

legend("bottomright", bty = "n", cex = 0.85,
       pch = c(19, 17), col = c("navy", "black"),
       legend = c("Exp 1: JTK-filtered top-K pool",
                  "Exp 2: published gene sets"))
par(op); dev.off()

#--- master summary table --------------------------------------------
master <- rbind(
	data.frame(experiment = "Exp1-JTKpool",
	           label = paste0("top", e1$K),
	           nGene = e1$nGene, nSel = e1$nSel,
	           MAE = e1$MAE, nAUC = e1$nAUC, pct2h = e1$pct2h),
	data.frame(experiment = "Exp2-geneset",
	           label = exp2.pooled$set,
	           nGene = exp2.pooled$nGene, nSel = exp2.pooled$nSel,
	           MAE = exp2.pooled$MAE, nAUC = exp2.pooled$nAUC,
	           pct2h = exp2.pooled$pct2h)
)
write.csv(master, "output/master_summary.csv", row.names = FALSE)

cat("\n=====================================================\n")
cat("Extension complete. Outputs in extension/output/:\n")
cat("  Exp1_JTK_prefilter_results.csv\n")
cat("  Exp1_tradeoff.pdf\n")
cat("  Exp2_geneset_substitution_results.csv\n")
cat("  Exp2_geneset_nAUC.pdf\n")
cat("  Fig_tradeoff_combined.pdf   <- headline deliverable\n")
cat("  master_summary.csv\n")
cat("=====================================================\n")

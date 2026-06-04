#####################################################################
# Experiment 2 -- direct substitution of published gene sets
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Question: which published circadian gene set transfers best when
#           dropped into the TimeSignature elastic-net framework?
#
# Sets compared (all trained on GSE39445, validated on GSE48113 /
# GSE56931 / GSE113883 with two blood draws per subject):
#   * TScore18      -- Braun 2018 core-18 predictors (Table 1)
#   * TimeMachine37 -- Huang & Braun 2024 JTK panel (37 genes)
#   * tauFisher     -- Duan & Ngo 2024 core-clock + top-10 rhythmic
#   * clockOnly     -- canonical core-clock genes only
#   * intersection  -- genes shared by the three published sets
# A full-pool model (all 7,615 genes) is included as the reference.
#
# Assumes 00_setup.R has been sourced.
#####################################################################

if (!exists("jtkRanked")) source("00_setup.R")

#--- assemble the gene sets (literature + JTK-derived panels) --------
sets <- buildGeneSets(jtk, available = rownames(all.expr),
                      nTimeMachine = 37, nTauTop = 10)
sets$fullPool <- rownames(all.expr)         # original TimeSignature reference

cat("\n--- Experiment 2 gene sets (size, members) ---\n")
for (nm in names(sets)) {
	g <- sets[[nm]]
	cat(sprintf("  %-13s n=%4d", nm, length(g)))
	if (length(g) <= 20) cat("  : ", paste(g, collapse = ", "))
	cat("\n")
}

#--- train + validate each set ---------------------------------------
exp2 <- do.call(rbind, lapply(names(sets), function(nm) {
	ts  <- trainTSpool(trainWSN, trainSubjs, trainTimes,
	                   genes = sets[[nm]], alpha = 0.5, s = NULL, seed = 194)
	res <- evalAcrossStudies(ts, calibAll, all.meta,
	                         studies = c("V1", "V2", "V3"))
	res$set <- nm
	cat(sprintf("  %-13s nGene=%4d selected=%2d  Valid.all: nAUC=%.3f MAE=%.2f %%2h=%.0f\n",
	            nm, length(sets[[nm]]), ts$nSelected,
	            res$nAUC[res$study == "Valid.all"],
	            res$MAE[res$study == "Valid.all"],
	            res$pct2h[res$study == "Valid.all"]))
	res
}))

write.csv(exp2, "output/Exp2_geneset_substitution_results.csv", row.names = FALSE)

#--- summary table (pooled validation) -------------------------------
exp2.pooled <- exp2[exp2$study == "Valid.all", ]
exp2.pooled <- exp2.pooled[order(exp2.pooled$nGene), ]
cat("\n--- Experiment 2: pooled validation metrics by gene set ---\n")
print(exp2.pooled[, c("set", "nGene", "nSel", "MAE", "nAUC", "pct2h")],
      row.names = FALSE, digits = 3)

#--- bar chart of nAUC per set ---------------------------------------
pdf("output/Exp2_geneset_nAUC.pdf", width = 7.5, height = 5)
op <- par(mar = c(7, 4.2, 3, 1))
bp <- barplot(exp2.pooled$nAUC, names.arg = sprintf("%s\n(n=%d)",
              exp2.pooled$set, exp2.pooled$nGene), las = 2,
              ylim = c(0, 1), col = "steelblue", ylab = "nAUC",
              main = "Exp 2: published gene sets in the TS framework")
abline(h = 0.80, lty = 2, col = "red")
text(bp, exp2.pooled$nAUC, sprintf("%.2f", exp2.pooled$nAUC), pos = 3, cex = 0.8)
par(op); dev.off()

best <- exp2.pooled[which.max(exp2.pooled$nAUC), ]
cat(sprintf("\nBest-transferring set: %s (n=%d genes, nAUC=%.3f, MAE=%.2f h)\n",
            best$set, best$nGene, best$nAUC, best$MAE))
cat("Exp 2 done -> output/Exp2_geneset_substitution_results.csv, output/Exp2_geneset_nAUC.pdf\n")

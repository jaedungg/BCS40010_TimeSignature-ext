#####################################################################
# Shared setup for the TimeSignature gene-selection extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
#
# Sourcing this file:
#   * loads the bundled TSexampleData (7,615 common genes, 4 studies),
#   * builds the within-subject-normalised training data (GSE39445),
#   * runs the JTK_Cycle rhythmicity ranking on the training data,
#   * builds the two-sample antipodal-calibration matrix for ALL
#     samples once (it is gene-independent, so each experiment just
#     subsets its rows).
#
# Run with the working directory set to TimeSignatR/extension/.
#
# Datasets:  TrTe = GSE39445 (Moller-Levet, train/test) -> TRAINING
#            V1   = GSE48113 (Archer)        -> validation
#            V2   = GSE56931 (Arnardottir)   -> validation
#            V3   = GSE113883 (RNA-seq)      -> validation
######################################################################

suppressMessages({
	library(limma)
	library(glmnet)
})

source("../R/TimeStampFns.R")
source("../R/JTKfns.R")
source("../R/GeneSets.R")
source("../R/TSextFns.R")

# exact-reproducibility RNG, as in the paper / TSexample.R
RNGversion("3.5.1")
set.seed(194)

load("../example/DATA/TSexampleData.Rdata")   # -> all.expr, all.meta

#---------------------------------------------------------------------
# Training data: the GSE39445 training half (all.meta$train == 1),
# within-subject normalised across all of its timepoints.
#---------------------------------------------------------------------
trainDat   <- all.expr[, all.meta$train == 1]
trainSubjs <- all.meta[all.meta$train == 1, "ID"]
trainTimes <- all.meta[all.meta$train == 1, "LocalTime"]
trainWSN   <- recalibrateExprs(trainDat, trainSubjs)     # genes x samples

#---------------------------------------------------------------------
# JTK_Cycle rhythmicity ranking on the training data (computed once).
# Used by both experiments: Exp 1 sweeps the top-K rhythmic genes;
# Exp 2 builds the TimeMachine / tauFisher panels from this ranking.
#---------------------------------------------------------------------
cat("Running JTK_Cycle on training data ...\n")
jtk <- jtkCycle(trainWSN, trainTimes, period = 24, nPhase = 24, binHours = 1)
jtkRanked <- jtk$gene[order(jtk$ADJ.P, -abs(jtk$tau))]   # most-rhythmic first

#---------------------------------------------------------------------
# Two-sample antipodal calibration for every sample (gene-independent).
# Each experiment predicts on calibAll[poolGenes, studySamples].
#---------------------------------------------------------------------
cat("Building two-sample antipodal calibration ...\n")
calibAll <- twoPointCalibrate(all.expr, all.meta, timeCol = "LocalTime")

cat("Setup complete: ", nrow(all.expr), " genes, ",
    ncol(all.expr), " samples, studies = ",
    paste(names(table(all.meta$study)), collapse = "/"), "\n", sep = "")

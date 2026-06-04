# TimeSignatR
TimeSignature R package v1.0

This package is currently under development to make it user-friendly.  In the interim, the script R/TimeStampFns.R contains all the necessary functions for carrying out a TimeSignature analysis, and the example/ subdirectory contains example code to train and apply TimeSignature and reproduce Fig 1 from the paper.  The user is encouraged to read all scripts prior to use.

Going forward, the R/TimeStampFns.R script will be turned into a fully-documented package that can be installed using the usual R tools.

---

## BCS40010 course extension

The `extension/` directory adds a gene-selection study on top of the original package: **does a JTK_Cycle rhythmicity pre-filter let TimeSignature use a much smaller predictor set while keeping clinically acceptable accuracy (nAUC ≥ 0.80, MAE < 2.5 h)?** 
New code lives in `R/JTKfns.R`, `R/GeneSets.R`, `R/TSextFns.R` (re-using the unchanged `R/TimeStampFns.R` pipeline).
Run with `cd extension && Rscript run_extension.R`. See `extension/README.md` for the methods, results and the headline gene-count vs accuracy trade-off figure.


# TimeSignature extension — Gene-selection optimization via JTK_Cycle pre-filtering

**BCS40010 project**

This directory extends the original `TimeSignatR` package (Braun et al. 2018,
*PNAS*) following the One-Page Synopsis. The base method runs an elastic net
over all 7,768 (here 7,615 common) genes and recovers ~40 predictors *without*
any prior periodicity filter. Later methods — **TimeMachine** (Huang & Braun
2024) and **tauFisher** (Duan & Ngo 2024) — insert a **JTK_Cycle rhythmicity
pre-filter** before regression. We ask:

> Can JTK_Cycle pre-filtering give TimeSignature a **smaller** predictor set
> that still meets clinically acceptable accuracy (**nAUC ≥ 0.80, MAE < 2.5 h**)?

All models **train on GSE39445** (Möller-Levet training half) and **validate on
GSE48113 / GSE56931 / GSE113883**, using **two blood draws per subject**
(antipodal calibration), exactly as in Fig 1 of the paper.

---

## What was added (new / modified code)

| File | Role |
|------|------|
| `../R/JTKfns.R` | **New.** Compact, fully-vectorised JTK_Cycle-style rhythmicity test (Kendall-τ cosine-template search, phase-Bonferroni `ADJ.P`, BH-FDR). |
| `../R/GeneSets.R` | **New.** Published gene-set definitions: TimeSignature core-18 (Braun Table 1), canonical core-clock genes, and builders for the TimeMachine / tauFisher panels. |
| `../R/TSextFns.R` | **New.** Helpers that re-use the original pipeline: `twoPointCalibrate`, `trainTSpool` (elastic net on a restricted pool), `evalMetrics` (MAE / nAUC / %within-2h), `evalAcrossStudies`, `tradeoffPlot`. |
| `00_setup.R` | Loads data, builds within-subject-normalised training data, runs the JTK ranking, builds the antipodal-calibration matrix. |
| `Exp1_JTK_prefilter.R` | **Experiment 1** — JTK pre-filter then elastic net, swept over pool size. |
| `Exp2_geneset_substitution.R` | **Experiment 2** — drop published gene sets into the TS framework. |
| `run_extension.R` | Master driver; produces the combined trade-off figure + summary table. |

The original `R/TimeStampFns.R` and `example/` are **unchanged** — the extension
only re-uses their functions.

## How to run

```r
# from TimeSignatR/extension/
Rscript run_extension.R          # both experiments + combined figure
# or individually:
Rscript Exp1_JTK_prefilter.R
Rscript Exp2_geneset_substitution.R
```
Outputs are written to `extension/output/` (CSV tables + PDF figures).
Requires R (tested on 4.6.0) with `glmnet` and `limma`; uses the bundled
`example/DATA/TSexampleData.Rdata`. RNG is pinned (`RNGversion("3.5.1")`,
`set.seed(194)`) for reproducibility.

---

## Methods (brief)

* **JTK_Cycle filter** (`jtkCycle`): each gene's within-subject-normalised
  training profile is binned to circadian time and compared, by Kendall's τ,
  against 24 cosine phase templates (period 24 h). The optimal-phase Kendall S
  gives a phase-Bonferroni `ADJ.P`; genes are ranked most-rhythmic first.
* **Restricted elastic net** (`trainTSpool`): the unchanged `trainTimeStamp`
  (α = 0.5, mgaussian on sin/cos clock coordinates) is fit on the chosen gene
  pool; penalty λ is CV-selected (`lambda.min`).
* **Two-draw evaluation** (`twoPointCalibrate` + `evalAcrossStudies`): each test
  sample is centred on the mean of its ~12 h-antipodal partner, then scored by
  MAE, nAUC (normalised area under the |error| CDF, as in Fig 1) and the
  fraction predicted within ±2 h.

**Note on JTK thresholds.** On the densely-sampled (~213-sample) training set the
absolute JTK p-values are liberal (p < 0.1 selects thousands of genes), so only
the *ranking* is used to define the pre-filtered pool. The sweep is therefore
driven by pool size *K* (equivalently, progressively stricter rhythmicity rank),
which directly yields the synopsis deliverable — the gene-count vs accuracy
curve. The TimeMachine/tauFisher panels are reconstructed from each method's own
stated procedure (the exact TimeMachine SI Table S2 symbols are not bundled).

---

## Results (pooled validation, this run)

**Experiment 1 — JTK pre-filtered pool size vs accuracy**

| pool K | predictors selected | MAE (h) | nAUC | %≤2h |
|-------:|--------:|-----:|-----:|----:|
| 10 | 10 | 2.40 | 0.810 | 51 |
| 20 | 19 | 2.20 | 0.827 | 59 |
| 30 | 26 | 2.13 | **0.833** | 62 |
| 50 | 42 | 2.11 | 0.834 | 59 |
| 7615 (full) | 173 | 1.89 | 0.853 | 62 |

→ A JTK-filtered pool of only **~30 rhythmic genes** reaches **nAUC 0.833 /
MAE 2.13 h**, within ~0.02 nAUC of the full 7,615-gene model while screening
**< 0.5 %** of the transcriptome. Even **10 genes clear the nAUC ≥ 0.80 floor.**
Pre-filtering does not *beat* the full elastic net, but reaches near-equivalent
accuracy from a drastically smaller candidate pool.

**Experiment 2 — published gene sets in the TS framework**

| gene set | n | MAE (h) | nAUC | %≤2h |
|----------|--:|-----:|-----:|----:|
| intersection (DDIT4, PER1, NR1D1, NR1D2) | 4 | 2.12 | 0.834 | 58 |
| clock-only | 8 | 2.61 | 0.792 | 52 |
| tauFisher | 17 | 2.19 | 0.828 | 55 |
| TimeSignature core-18 | 18 | 2.13 | 0.833 | 57 |
| **TimeMachine-37** | 37 | **1.99** | **0.844** | 63 |
| full pool | 7615 | 1.89 | 0.853 | 62 |

→ Among published sets, **TimeMachine-37 transfers best** (nAUC 0.844). The
**4-gene intersection** of the three sets already clears the clinical floor
(nAUC 0.834), whereas **canonical clock genes alone fall short** (0.792) — the
data-driven rhythmic genes, not the textbook clock loop, carry the signal.

**Headline figure:** `output/Fig_tradeoff_combined.pdf` — the Exp 1 trade-off
curve with the Exp 2 gene sets overlaid.

---

## References
- Braun et al. (2018) *PNAS* **115**, E9247–E9256.
- Huang & Braun (2024) *PNAS* **121**, e2308114120.
- Duan & Ngo (2024) *Nat. Commun.* **15**, 3840.
- Hughes, Hogenesch & Kornacker (2010) *J. Biol. Rhythms* **25**, 372–380 (JTK_Cycle).

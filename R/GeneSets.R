#####################################################################
# Published circadian gene sets used in the TimeSignature extension
#
# BCS40010 project (Kim Na-hye / Lee Jae-hyun)
# Experiment 2: "Which published gene set transfers best to the
#                TimeSignature framework?"
#
# The fixed (literature) sets are coded here as gene-symbol vectors.
# The data-derived sets (TimeMachine's JTK panel, tauFisher's
# rhythmic genes) are built at run time from the JTK_Cycle ranking,
# because they are defined procedurally by their authors rather than
# as a fixed list (and the exact TimeMachine SI Table S2 symbols are
# not bundled with this repository).  See buildGeneSets() below.
######################################################################


#---------------------------------------------------------------------
# TimeSignature "core 18"  (Braun et al. 2018, PNAS, Table 1)
# The 18 genes selected as predictors a *majority* of the time
# (selection frequency >= 0.50) across the authors' 12 repeated runs.
#---------------------------------------------------------------------
TScore18 <- c(
	"DDIT4", "GHRL", "PER1", "EPHX2", "GNG2", "IL1B", "DHRS13",
	"NR1D1", "ZNF438", "NR1D2", "CD38", "TIAM2", "CD1C", "LLGL2",
	"GZMB", "CLEC10A", "PDK1", "GPCPD1"
)

#---------------------------------------------------------------------
# Canonical core-clock genes (the transcriptional/translational
# feedback loop).  tauFisher (Duan & Ngo 2024) anchors its predictor
# on exactly these genes; TimeMachine likewise tracks them.
# Aliases: BMAL1 == ARNTL.  (PER2/PER3/CLOCK/NPAS2/CIART are not in
# the 7,615-gene common set bundled with TimeSignatR.)
#---------------------------------------------------------------------
clockGenes <- c(
	"ARNTL", "BMAL1", "DBP", "NR1D1", "NR1D2",
	"PER1", "PER2", "PER3", "CRY1", "CRY2",
	"CLOCK", "NPAS2", "CIART", "TEF", "HLF"
)


#---------------------------------------------------------------------
# Build all Experiment-2 gene sets, restricted to the genes actually
# present on the expression matrix.
#
#   jtk        : data.frame from jtkCycle() (training-data ranking)
#   available  : rownames(all.expr) -- the assayable gene universe
#   nTimeMachine: size of the TimeMachine-style JTK panel (default 37,
#                 matching the "37 genes at JTK P<0.1" of Huang &
#                 Braun 2024)
#   nTauTop    : number of top rhythmic genes tauFisher adds to the
#                core-clock anchor (default 10, per Duan & Ngo 2024)
#
# Returns a named list of gene-symbol vectors:
#   TScore18      - the literature core-18 set
#   TimeMachine37 - core-clock genes + the most rhythmic genes, JTK-ranked
#   tauFisher     - core-clock genes + top-10 rhythmic genes
#   clockOnly     - core-clock genes only
#   intersection  - genes common to TScore18, TimeMachine37 and tauFisher
#---------------------------------------------------------------------
buildGeneSets <- function(jtk, available, nTimeMachine = 37, nTauTop = 10) {

	keep   <- function(g) intersect(g, available)
	ranked <- jtk$gene[order(jtk$ADJ.P, -abs(jtk$tau))]   # most-rhythmic first
	ranked <- intersect(ranked, available)

	core   <- keep(clockGenes)

	# TimeMachine: the authors take all genes passing JTK P<0.1; on
	# this densely-sampled data that is a large pool, so we reproduce
	# their *published panel size* (37) by taking the core-clock genes
	# plus the most rhythmic genes up to 37 total.
	tm <- unique(c(core, ranked))
	tm <- head(tm, nTimeMachine)

	# tauFisher: core-clock anchor + the top-N rhythmic genes.
	tau <- unique(c(core, head(ranked, nTauTop)))

	sets <- list(
		TScore18      = keep(TScore18),
		TimeMachine37 = tm,
		tauFisher     = tau,
		clockOnly     = core
	)
	sets$intersection <- Reduce(intersect,
		list(sets$TScore18, sets$TimeMachine37, sets$tauFisher))
	sets
}

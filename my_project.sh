#!/usr/bin/env bash
# =============================================================================
# config.sh -- User-editable pipeline settings
# Copy this file for each project and fill in your paths:
#   cp config.sh my_project.sh
#   bash run_pipeline.sh my_project.sh
# =============================================================================

# -- Project identity --------------------------------------------------------
# Unique name for this analysis run (letters, numbers, underscores recommended).
# The pipeline creates: ${PROJECT_OUT}/${PROJECT_NAME}/  for machine-facing
# per-project outputs (scanpy, downstream, logs). Human-facing outputs land in
# ${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}/ and ${SUMMARY_TABLES_OUT}/${PROJECT_NAME}/
# instead -- kept out of PROJECT_OUT so they're easy to browse across projects
# without descending into each project's own directory.
# Cell Ranger outputs live in CELLRANGER_OUT (shared across projects) so that
# multiple projects sharing the same samples never re-run cellranger count.
PROJECT_NAME="5k_human_PBMC"

# -- Input mode (1) FASTQ ----------------------------------------------------
# (1) "fastq"  (default): run Cell Ranger on FASTQ_DIR to produce feature matrices.
# FASTQ_DIR, GENOME_REF are only used when INPUT_MODE="fastq".
INPUT_MODE="fastq"
FASTQ_DIR="${RUNDIR}/0.data/sctranscriptome/fastq/5k_Human_Donor"                                      # Directory containing per-sample FASTQ subdirs
GENOME_REF="${RUNDIR}/0.refseq/refdata-gex-GRCh38-2024-A"         # Cell Ranger genome reference directory

# -- Input mode (2) MATRIX ---------------------------------------------------
# (2) "matrix": skip Cell Ranger entirely and load pre-computed feature-barcode
#           matrices from MATRIX_DIR instead -- use this when matrices were
#           already generated elsewhere (a different pipeline, a GEO/public
#           dataset download, etc.) and you only want Stages 3-7.
#           Layout expected under MATRIX_DIR (flat, one folder per sample):
#             MATRIX_DIR/<sample>/matrix.mtx.gz
#             MATRIX_DIR/<sample>/features.tsv.gz
#             MATRIX_DIR/<sample>/barcodes.tsv.gz
#           Each <sample> must match a "sample" value in METADATA_FILE below.
# INPUT_MODE="matrix"
# MATRIX_DIR="${RUNDIR}/1.scrnaseq/data/matrix/plastics_test"

# -- Paths -------------------------------------------------------------------
METADATA_FILE="${RUNDIR}/1.scrnaseq/metadata/5k_human_pbmc.tsv"                    # See format note below
OUT_DIR="${RUNDIR}/1.scrnaseq"                         # Base directory for all pipeline output

# -- Sample metadata format --------------------------------------------------
# TSV with a "sample" column -- run_pipeline.sh reads this column directly to
# get the sample list, in the order rows appear here. There is no separate
# SAMPLES=() array to keep in sync by hand: each "sample" value must match a
# subdirectory in FASTQ_DIR (INPUT_MODE="fastq") or MATRIX_DIR
# (INPUT_MODE="matrix"). In fastq mode, the FASTQ filename prefix (used for
# cellranger --sample) is auto-detected from the .fastq.gz filenames inside
# that directory -- no manual mapping needed there either.
# All other columns (condition, sex, batch, ...) are merged into per-cell obs.
#
# Example:
#   sample    condition  sex  batch
#   Sample1   control    F    1
#   Sample2   treated    M    2

# Pick a control condition from the metadata "condition" column.
CONTROL_CONDITION="control"

# -- Tool paths --------------------------------------------------------------
# Set these if cellranger / conda are not already on your PATH.
# The pipeline injects them into PATH at runtime (nothing is written to .bashrc).
# Leave empty ("") to rely on the system PATH instead.
CELLRANGER_DIR=""   # e.g. /opt/cellranger-10.0.0  (dir containing the binary)
CONDA_BIN=""        # e.g. /home/user/miniforge3/bin

# -- Resources ---------------------------------------------------------------
# THREADS/MEM_GB are the pipeline-wide TOTAL budget; Stage 1 divides it across
# concurrent Cell Ranger jobs via CELLRANGER_THREADS_PER_SAMPLE/
# CELLRANGER_MEM_GB_PER_SAMPLE below.
THREADS=16
MEM_GB=100

CELLRANGER_THREADS_PER_SAMPLE=16
CELLRANGER_MEM_GB_PER_SAMPLE=64

# -- QC thresholds -----------------------------------------------------------
# MIN_GENES/MAX_GENES are a coarse pre-filter, applied before the MAD-based
# outlier detection below. Real cell filtering is MAD-based, not these fixed
# thresholds.
MIN_GENES=200
MAX_GENES=6000
MAX_MITO=20    # percent mitochondrial reads; hard cap combined with N_MADS_MT below
N_MADS=5       # MAD threshold for total_counts/n_genes/pct_counts_in_top_20_genes outliers
N_MADS_MT=3    # MAD threshold for pct_counts_mt specifically (tighter than N_MADS)
REMOVE_DOUBLETS=false  # true to drop predicted_doublet cells during Stage 3b QC (see 03b_doublet_summary.pdf first)

# -- Analysis parameters -----------------------------------------------------
N_HVGS=3000
N_PCS=50
N_NEIGHBORS=15
LEIDEN_RES=0.5   # clustering resolution: 0.3 = coarse, 1.2 = fine

# -- Conda environment names -------------------------------------------------
PYTHON_ENV="scrna_py"
R_ENV="scrna_r"

# -- Run flags ---------------------------------------------------------------
# Override on the command line without editing this file:
#   SETUP_ENVS=true bash run_pipeline.sh config.sh
#   SKIP_CELLRANGER=true bash run_pipeline.sh config.sh
#   SKIP_AMBIENT_QC=true bash run_pipeline.sh config.sh
#   PREPROCESSING_ONLY=true bash run_pipeline.sh config.sh
SETUP_ENVS=false          # true on first run to create the two conda envs
SKIP_AMBIENT_QC=false     # true to skip Stage 2 (SoupX ambient RNA correction + scDblFinder doublet detection)
SKIP_CELLRANGER=false     # true if filtered_feature_bc_matrix dirs already exist
PREPROCESSING_ONLY=false  # true to stop after Stage 5 (through clustering/annotation) and skip Stage 6 (scCODA) + Stage 7 (R downstream)

#!/usr/bin/env Rscript
# =============================================================================
# downstream.R — Stage 7: Pseudobulk DEG (7a: DESeq2), trajectory (7b:
#                Monocle3), cell–cell communication (7c: LIANA), gene
#                program (7d: NMF), Gene Set Enrichment Analysis (7e:
#                clusterProfiler), kNN-neighborhood differential abundance
#                (7f: Milo). Entry point: CLI parsing, loads scanpy's
#                exported h5ad/CSV, validates CONTROL_CONDITION, then
#                source()s r/*.R and calls each Stage 7 submodule below in
#                order. Minimal packages: DESeq2, Monocle3, liana, NMF,
#                clusterProfiler, miloR — all Bioconductor/CRAN
# =============================================================================
suppressPackageStartupMessages({
    library(optparse)
    library(Matrix)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
})

# Self-locate this script's directory so r/*.R can be source()'d by an
# absolute path regardless of the caller's working directory -- run_pipeline.sh
# always invokes this via an absolute path, but Rscript (unlike Python) has no
# built-in "directory of the running script" concept, hence the commandArgs()
# idiom below.
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
SCRIPT_DIR <- if (length(.file_arg) > 0) dirname(normalizePath(.file_arg)) else getwd()

# ── CLI ───────────────────────────────────────────────────────────────────────
# Only paths are needed — condition and cell type columns are auto-detected
# from obs_metadata.csv, which inherits all columns from sample_metadata.tsv.
option_list <- list(
    make_option("--scanpy-dir",  type = "character"),
    make_option("--out-dir",     type = "character"),
    make_option("--viz-dir",     type = "character", default = NULL,
                help = "Directory for summary figure PDFs [default: <out-dir>/../summary_figures]"),
    make_option("--tables-dir",  type = "character", default = NULL,
                help = "Directory for summary data tables (TSV) [default: <out-dir>/../summary_tables]"),
    make_option("--threads",     type = "integer",   default = 1,
                help = paste("Worker count for parallel::mclapply() over cell types in Stage",
                              "7a/7e (DESeq2, GSEA) [default: %default]")),
    make_option("--control-condition", type = "character", default = "",
                help = paste("Reference/baseline level of the auto-detected condition column.",
                              "Required if that column has more than 2 levels."))
)
opt <- parse_args(OptionParser(option_list = option_list))

out_dir <- opt[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

viz_dir <- if (!is.null(opt[["viz-dir"]])) opt[["viz-dir"]] else
               file.path(dirname(out_dir), "summary_figures")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

tables_dir <- if (!is.null(opt[["tables-dir"]])) opt[["tables-dir"]] else
               file.path(dirname(out_dir), "summary_tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Logging (log_msg/log_warn/log_error), the uncaught-error handler, shared
# visual style (PALETTE_QUALITATIVE/COLOR_*/FIGURE_THEME), save_plot,
# verify_files, safe_name, and condition-column auto-detection
# (PIPELINE_COLS/detect_condition_col) -- all defined once here since every
# Stage 7 submodule below depends on them.
source(file.path(SCRIPT_DIR, "r", "common.R"))

# ── Load data exported by scanpy ──────────────────────────────────────────────
log_msg("Loading scanpy outputs...")

obs <- read_csv(file.path(opt[["scanpy-dir"]], "obs_metadata.csv"),
                show_col_types = FALSE)

# Load the raw count matrix from .h5ad via zellkonverter (or fall back to csv)
counts_path <- file.path(opt[["scanpy-dir"]], "counts_matrix.csv")
if (file.exists(counts_path)) {
    log_msg("Loading counts from CSV export...")
    counts_raw <- read_csv(counts_path, show_col_types = FALSE) |>
        tibble::column_to_rownames("cell_barcode") |>
        as.matrix() |>
        t()
    counts_lognorm <- log1p(counts_raw)  # approximate; prefer .h5ad path
} else {
    if (!requireNamespace("zellkonverter", quietly = TRUE))
        stop("Install zellkonverter: BiocManager::install('zellkonverter')")
    suppressPackageStartupMessages(library(zellkonverter))
    suppressPackageStartupMessages(library(SingleCellExperiment))
    sce <- readH5AD(file.path(opt[["scanpy-dir"]], "annotated.h5ad"))
    # adata.layers["counts"] → stored as assay "counts" by zellkonverter
    # adata.X (log-normalized) → stored as assay "X" by zellkonverter
    if (!"counts" %in% assayNames(sce))
        stop("'counts' assay not found in .h5ad. Check that adata.layers['counts'] was saved.")
    counts_raw     <- assay(sce, "counts")     # raw integers — for DESeq2
    counts_lognorm <- assay(sce, "X")          # log-normalized — for LIANA, NMF
    obs <- as.data.frame(colData(sce))
}

log_msg("Cells loaded: ", nrow(obs))

# -- Verify inputs -----------------------------------------------------------
if (nrow(obs) == 0)
    stop("[VERIFY FAIL] obs_metadata.csv is empty — was scanpy_analysis.py (Stage 3-5) completed?")
if (ncol(counts_raw) != nrow(obs))
    stop(sprintf("[VERIFY FAIL] Count matrix has %d cells but obs has %d rows — mismatch.",
                 ncol(counts_raw), nrow(obs)))
log_msg("[VERIFY SUCCESS] Inputs OK — ", nrow(obs), " cells, ", nrow(counts_raw), " genes")

# ── Auto-detect condition column ──────────────────────────────────────────────
# Cell type column is always "cell_type" — written by scanpy_analysis.py.
ctype_col <- "cell_type"
if (!ctype_col %in% names(obs))
    stop("'cell_type' column not found in obs_metadata.csv. Was annotation completed?")

cond_col <- detect_condition_col(obs, PIPELINE_COLS)

# as.character() is essential here, not cosmetic: obs[[cond_col]] arrives as
# an R factor via the annotated.h5ad/zellkonverter path (colData's categorical
# columns round-trip as factors, not character). A bare factor value passed as
# relevel()'s `ref=` skips its is.character() branch and is treated as an
# already-numeric level index instead -- silently wrong or an outright crash
# ("ref must be in 1L:nlev") depending on whether the index happens to be
# in-range. Worse, c(character_vector, factor_value) coerces the factor to its
# *integer code* as a string, which would silently corrupt DESeq2's
# `contrast=` argument. Every level extracted from obs/pb_meta/sample_meta
# must go through as.character() before being used as a relevel/contrast
# value -- this is the single point where that's enforced for condition_levels
# (r/stage7a_deseq2.R and r/stage7f_milo.R independently re-apply the same
# as.character() rule where they each extract their own level sets).
condition_levels <- as.character(unique(obs[[cond_col]]))
log_msg("Conditions: ", paste(condition_levels, collapse = ", "))
log_msg("Cell types: ", paste(unique(obs[[ctype_col]]), collapse = ", "))

# ── Control-condition validation ──────────────────────────────────────────────
# With >2 levels, DESeq2's results() would otherwise silently return only the
# last-vs-first-alphabetical coefficient (whichever level sorts first
# alphabetically, not necessarily the real control) -- the same class of
# silent-wrong-grouping bug PIPELINE_COLS (r/common.R) already guards against.
# Require an explicit choice instead of guessing. Optional for exactly 2
# levels: if set, still relevels for a correctly-signed fold-change; if
# unset, keeps the prior arbitrary-but-deterministic reference unchanged (no
# regression for existing 2-level runs).
control_condition <- trimws(opt[["control-condition"]])
if (length(condition_levels) > 2 && control_condition == "") {
    stop(sprintf(
        "[VERIFY FAIL] Condition column '%s' has %d levels (%s) but --control-condition is not set. Set CONTROL_CONDITION in my_project.sh to one of the levels above.",
        cond_col, length(condition_levels), paste(condition_levels, collapse = ", ")))
}
if (control_condition != "" && !(control_condition %in% condition_levels)) {
    stop(sprintf("[VERIFY FAIL] CONTROL_CONDITION='%s' is not among the detected levels: %s",
                 control_condition, paste(condition_levels, collapse = ", ")))
}

# ── Stage 7 submodules ────────────────────────────────────────────────────────
source(file.path(SCRIPT_DIR, "r", "stage7a_deseq2.R"))
source(file.path(SCRIPT_DIR, "r", "stage7b_monocle3.R"))
source(file.path(SCRIPT_DIR, "r", "stage7c_liana.R"))
source(file.path(SCRIPT_DIR, "r", "stage7d_nmf.R"))
source(file.path(SCRIPT_DIR, "r", "stage7e_gsea.R"))
source(file.path(SCRIPT_DIR, "r", "stage7f_milo.R"))
source(file.path(SCRIPT_DIR, "r", "stage7_summary.R"))

# =============================================================================
# Main
# =============================================================================
log_msg("Starting Stage 7 downstream analysis...")
log_msg("Threads: ", opt[["threads"]])
log_msg("Control condition: ", if (control_condition == "") "<not set>" else control_condition)

deseq2_out <- run_deseq2(counts_raw, obs, cond_col, ctype_col, control_condition, tables_dir, viz_dir,
                         threads = opt[["threads"]])
all_results  <- deseq2_out$results
contrasts_df <- deseq2_out$contrasts

run_monocle3(counts_raw,  obs, ctype_col, tables_dir, viz_dir, opt[["scanpy-dir"]])
run_liana(counts_lognorm, obs, ctype_col, tables_dir, viz_dir)
run_nmf(counts_lognorm,   tables_dir, k = 8)
if (exists("sce")) {
    run_milo(sce, obs, cond_col, ctype_col, control_condition, tables_dir, viz_dir)
} else {
    log_msg("Milo requires PCA coordinates from annotated.h5ad — skipping ",
            "(counts_matrix.csv fallback doesn't include them)")
}
run_gsea(all_results, contrasts_df, tables_dir, viz_dir, threads = opt[["threads"]])

make_downstream_summary(all_results, contrasts_df, viz_dir)

log_msg("All Stage 7 analyses complete.")
log_msg("Tables: ", tables_dir)
log_msg("Visualizations: ", viz_dir)

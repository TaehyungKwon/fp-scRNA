#!/usr/bin/env Rscript
# =============================================================================
# preprocess_qc.R — Stage 2: per-sample ambient RNA correction (2a: SoupX) +
# doublet detection (2b: scDblFinder), run once per sample before Stage 3
# (scanpy) combines them.
#
# Doublet detection must run per-sample, never on pooled/multi-batch data
# (artificial doublets are only meaningful within one batch) -- that's why
# this is a standalone per-sample stage rather than folded into Stage 3's
# already-combined AnnData.
# =============================================================================
suppressPackageStartupMessages({
    library(optparse)
    library(Matrix)
})

option_list <- list(
    make_option("--sample",       type = "character"),
    make_option("--filtered-dir", type = "character",
                help = "Cell Ranger filtered_feature_bc_matrix directory"),
    make_option("--raw-dir",      type = "character", default = NULL,
                help = "Cell Ranger raw_feature_bc_matrix directory, for SoupX's ambient profile. Omit to skip ambient RNA correction."),
    make_option("--out-dir",      type = "character",
                help = "Directory to write soupx_corrected_matrix/ and doublet_calls.tsv into")
)
opt <- parse_args(OptionParser(option_list = option_list))

.log_line <- function(tag, ...) {
    cat(sprintf("[%s] [%s]\t%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tag, paste0(...)))
}
log_msg  <- function(...) .log_line("INFO", ...)
log_warn <- function(...) .log_line("WARNING", ...)
options(error = function() {
    .log_line("ERROR", geterrmessage())
    quit(status = 1, save = "no")
})

dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(DropletUtils))

read_10x_deduped <- function(path) {
    sce <- read10xCounts(path, col.names = TRUE)
    mat <- as(counts(sce), "CsparseMatrix")
    rownames(mat) <- make.unique(rowData(sce)$Symbol)
    colnames(mat) <- colnames(sce)
    list(mat = mat, gene_id = rowData(sce)$ID, gene_symbol = rownames(mat))
}

log_msg("Reading filtered matrix: ", opt[["filtered-dir"]])
filt <- read_10x_deduped(opt[["filtered-dir"]])
corrected <- filt$mat  # falls through unchanged if ambient correction is skipped/unavailable
rho <- NA

if (!is.null(opt[["raw-dir"]])) {
    if (!requireNamespace("SoupX", quietly = TRUE) || !requireNamespace("scran", quietly = TRUE)) {
        log_warn("SoupX/scran not installed -- skipping ambient RNA correction")
    } else {
        suppressPackageStartupMessages({ library(SoupX); library(scran) })
        log_msg("Reading raw (unfiltered) matrix: ", opt[["raw-dir"]])
        raw <- read_10x_deduped(opt[["raw-dir"]])

        log_msg("Quick-clustering the filtered matrix for SoupX (SoupX's own results ",
                "aren't sensitive to cluster quality -- this is a cheap pre-clustering step, ",
                "not the pipeline's real clustering in Stage 5)")
        clusters <- tryCatch(
            as.character(scran::quickCluster(filt$mat)),
            error = function(e) {
                log_warn("quickCluster failed (", conditionMessage(e), ") -- using a single cluster")
                rep("1", ncol(filt$mat))
            }
        )
        names(clusters) <- colnames(filt$mat)

        log_msg("=== Stage 2a: SoupX -- estimating and correcting ambient RNA ===")
        channel <- SoupChannel(raw$mat, filt$mat, calcSoupProfile = TRUE)
        channel <- setClusters(channel, clusters)
        channel <- tryCatch(
            autoEstCont(channel, doPlot = FALSE),
            error = function(e) { log_warn("autoEstCont failed: ", conditionMessage(e)); NULL }
        )
        if (is.null(channel)) {
            log_warn("Ambient RNA correction failed -- using uncorrected counts")
        } else {
            corrected <- adjustCounts(channel, roundToInt = TRUE)
            rho <- channel$metaData$rho[1]
            log_msg("SoupX complete. Estimated contamination fraction: ", round(rho * 100, 2), "%")
        }
    }
} else {
    log_msg("No --raw-dir given -- skipping SoupX ambient RNA correction")
}

# Write the (possibly SoupX-corrected) matrix in the same flat 10x-mtx layout
# INPUT_MODE=matrix already expects, so scanpy_analysis.py reads it the same way.
out_matrix_dir <- file.path(opt[["out-dir"]], "soupx_corrected_matrix")
unlink(out_matrix_dir, recursive = TRUE)
write10xCounts(out_matrix_dir, corrected, gene.id = filt$gene_id, gene.symbol = filt$gene_symbol,
               version = "3", overwrite = TRUE)
log_msg("Corrected matrix written: ", out_matrix_dir)

# -- Stage 2b: Doublet detection (scDblFinder) -- always per-sample, see header note ---
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    log_warn("scDblFinder not installed -- skipping doublet detection")
} else {
    suppressPackageStartupMessages({ library(scDblFinder); library(SingleCellExperiment) })
    log_msg("=== Stage 2b: scDblFinder -- doublet detection ===")
    set.seed(123)
    sce_db <- scDblFinder(SingleCellExperiment(list(counts = corrected)))
    doublet_df <- data.frame(
        barcode            = colnames(sce_db),
        scDblFinder_score  = sce_db$scDblFinder.score,
        scDblFinder_class  = as.character(sce_db$scDblFinder.class)
    )
    out_doublet_path <- file.path(opt[["out-dir"]], "doublet_calls.tsv")
    write.table(doublet_df, out_doublet_path, sep = "\t", row.names = FALSE, quote = FALSE)
    n_doublet <- sum(doublet_df$scDblFinder_class == "doublet")
    log_msg("Doublet calls (", n_doublet, " of ", nrow(doublet_df), " flagged): ", out_doublet_path)
}

log_msg("Stage 2 complete: ", opt[["sample"]])

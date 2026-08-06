# =============================================================================
# r/common.R — shared infrastructure for downstream.R's Stage 7 submodules
# (r/stage7a_deseq2.R ... r/stage7f_milo.R, r/stage7_summary.R): logging,
# the shared visual style (mirrors py/common.py's PALETTE_QUALITATIVE/
# COLOR_* so Python- and R-generated figures read as one system), save_plot,
# output verification, the filename sanitizer, and condition-column
# auto-detection (also independently ported in differential_abundance.py --
# see PIPELINE_COLS's own comment below).
# =============================================================================

.log_line <- function(tag, ...) {
    cat(sprintf("[%s] [%s]\t%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tag, paste0(...)))
}
log_msg  <- function(...) .log_line("INFO", ...)     # kept as the name most call sites already use
log_info <- log_msg
log_warn  <- function(...) .log_line("WARNING", ...)
log_error <- function(...) .log_line("ERROR", ...)

# Uncaught stop()/errors otherwise print R's default "Error: message" or
# "Error in f(): message" with no tag -- route them through log_error first.
# Rscript (non-interactive) consults options("error") on any unhandled error.
options(error = function() {
    log_error(geterrmessage())
    quit(status = 1, save = "no")
})

# ── Shared visual style ──────────────────────────────────────────────────────
# Hex codes mirror py/common.py's PALETTE_QUALITATIVE / COLOR_* so Python-
# and R-generated figures read as one consistent visual system.
PALETTE_QUALITATIVE <- c(
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
)
COLOR_ACCENT <- "#D85A30"
COLOR_MUTED  <- "grey70"
FIGURE_DPI   <- 300

# Target page envelope for every summary figure: 180mm x 240mm (fits a
# double-column journal page), matching py/common.py's PAGE_W_IN/PAGE_H_IN.
PAGE_W_MM <- 180
PAGE_H_MM <- 240

# Real Arial isn't installed on Linux; Liberation Sans is the metrically-
# compatible open substitute (matches py/common.py's choice). R's plain
# pdf() device errors on unrecognized font family names ("invalid font
# type"), so save_plot() below uses the cairo_pdf device, which resolves
# fonts through the same fontconfig database `fc-list` reads.
# Additive theme layer -- append with `+ FIGURE_THEME` after each plot's own
# theme_classic(...)/etc. call so it only overrides text properties, not the
# whole theme.
FIGURE_THEME <- theme(
    text          = element_text(family = "Liberation Sans"),
    axis.text     = element_text(size = 6),
    legend.text   = element_text(size = 5),
    legend.title  = element_text(size = 6),
    legend.position = "bottom"
)

# Save a ggplot as PDF, sized in mm. `path`'s extension is ignored. Unlike
# matplotlib's bbox_inches="tight", ggsave does not auto-expand the canvas
# to fit overflowing artists, so requesting a size within the page envelope
# is sufficient to guarantee the output fits it.
save_plot <- function(plot, path, width_mm, height_mm, dpi = FIGURE_DPI) {
    if (width_mm > PAGE_W_MM + 1 || height_mm > PAGE_H_MM + 1) {
        log_warn("[VIZ] ", basename(path), ": requested size ", width_mm, "x", height_mm,
                "mm exceeds the ", PAGE_W_MM, "x", PAGE_H_MM, "mm page envelope")
    }
    out <- sub("\\.[^.]+$", ".pdf", path)
    ggsave(out, plot, width = width_mm, height = height_mm, units = "mm", dpi = dpi,
           device = cairo_pdf)
}

# Verify that a set of expected output files all exist after a stage completes.
# Logs a warning (does NOT stop) so optional modules that partially succeed
# are still reported clearly rather than silently failing.
verify_files <- function(paths, stage) {
    missing <- paths[!file.exists(paths)]
    if (length(missing) > 0) {
        log_warn("[VERIFY] ", stage, " — missing expected outputs: ",
                paste(basename(missing), collapse = ", "))
    } else {
        sizes <- vapply(paths, function(p) {
            sprintf("%s (%.1f KB)", basename(p), file.size(p) / 1024)
        }, character(1))
        log_msg("[VERIFY SUCCESS] ", stage, " OK — ", paste(sizes, collapse = ", "))
    }
}

# Shared filename sanitizer for both cell-type names and condition levels --
# used throughout Stage 7a/7e/7f's per-contrast output filenames.
safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# ── Condition-column auto-detection ───────────────────────────────────────────
# Technical columns produced by the pipeline — never candidates for "condition".
# Independently ported in differential_abundance.py's own PIPELINE_COLS -- the
# two copies must be kept in sync manually (this pipeline shares no code
# across languages by design, see codeshare/CLAUDE.md). Any new pipeline-
# generated obs column with 2-10 distinct values must be added to BOTH lists,
# or it can silently become the DEG/DA grouping variable instead of the real
# condition -- this actually happened once here when scDblFinder_class
# (singlet/doublet) was added without updating this list, and Stage 7 silently
# ran DESeq2/GSEA on "singlet vs doublet" instead of the real condition.
PIPELINE_COLS <- c(
    "sample", "cell_type", "leiden", "doublet_score", "predicted_doublet",
    "scDblFinder_score", "scDblFinder_class", "outlier", "mt_outlier",
    "n_genes_by_counts", "log1p_n_genes_by_counts",
    "total_counts", "log1p_total_counts", "pct_counts_in_top_20_genes",
    "total_counts_mt", "log1p_total_counts_mt", "pct_counts_mt",
    "total_counts_ribo", "log1p_total_counts_ribo", "pct_counts_ribo",
    "total_counts_hb", "log1p_total_counts_hb", "pct_counts_hb",
    "n_counts", "n_genes", "UMAP1", "UMAP2"
)

detect_condition_col <- function(obs, pipeline_cols) {
    # Candidates: columns from your sample_metadata.tsv that are NOT pipeline cols,
    # have 2–10 unique values, and vary across samples (not all cells same value).
    candidates <- obs |>
        dplyr::select(-dplyr::any_of(pipeline_cols)) |>
        dplyr::select(where(~ {
            n_unique <- dplyr::n_distinct(.)
            n_unique >= 2 && n_unique <= 10
        })) |>
        names()

    if (length(candidates) == 0)
        stop("Could not auto-detect a condition column. Check sample_metadata.tsv has a condition-like column.")

    if (length(candidates) > 1)
        log_msg("Multiple condition candidates found: [",
                paste(candidates, collapse = ", "),
                "] — using '", candidates[1], "'. Set manually if wrong.")

    log_msg("Condition column: '", candidates[1], "'")
    candidates[1]
}

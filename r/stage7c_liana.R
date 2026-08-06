# =============================================================================
# r/stage7c_liana.R — Stage 7c: Cell–cell communication (LIANA)
# =============================================================================
run_liana <- function(counts, obs, ctype_col, tables_dir, viz_dir) {
    if (!requireNamespace("liana", quietly = TRUE)) {
        log_msg("liana not available — skipping cell–cell communication")
        return(invisible(NULL))
    }
    suppressPackageStartupMessages({
        library(liana)
        library(SingleCellExperiment)
    })
    log_msg("=== Stage 7c: Cell–cell communication (LIANA) ===")

    sce <- SingleCellExperiment(assays = list(counts = counts, logcounts = log1p(counts)))
    colData(sce)[[ctype_col]] <- obs[[ctype_col]]

    liana_res <- liana_wrap(sce, idents_col = ctype_col)
    liana_res <- liana_aggregate(liana_res)
    results_path <- file.path(tables_dir, "liana_results.tsv")
    write_tsv(as.data.frame(liana_res), results_path)

    dotplot_path <- file.path(viz_dir, "07c_liana_dotplot.pdf")
    p <- liana_dotplot(liana_res, source_groups = unique(obs[[ctype_col]])[1:4],
                       target_groups = unique(obs[[ctype_col]])[1:4]) +
        FIGURE_THEME
    save_plot(p, dotplot_path, width_mm = 180, height_mm = 140)
    verify_files(c(results_path, dotplot_path), "Stage 7c (LIANA)")
    log_msg("Stage 7c complete.")
}

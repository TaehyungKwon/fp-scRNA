# =============================================================================
# r/stage7_summary.R — Stage 7 summary visualization (DEG counts per cell
# type/contrast). Cell-type composition per sample is already shown in
# py/reports.py's 05_clustering_summary.pdf (panel 4) -- not repeated here.
# =============================================================================
make_downstream_summary <- function(all_results, contrasts_df, viz_dir) {
    log_msg("=== Generating downstream summary visualization ===")

    if (length(all_results) == 0) {
        log_msg("[VIZ] No DEG results available for downstream summary — skipped")
        return(invisible(NULL))
    }

    deg_df <- contrasts_df |>
        dplyr::mutate(n_deg = vapply(key, function(k) {
            sum(all_results[[k]]$padj < 0.05 & abs(all_results[[k]]$log2FoldChange) > 1, na.rm = TRUE)
        }, integer(1)))

    # One file per non-control level -- same overflow-avoidance rationale as
    # umap_overview.pdf's split into umap_leiden.pdf/umap_cell_type.pdf:
    # cramming cell_type x contrast into one page stops being legible past a
    # couple of contrasts.
    for (lvl in unique(deg_df$level)) {
        lvl_df <- deg_df[deg_df$level == lvl, ] |> dplyr::arrange(dplyr::desc(n_deg))
        ctrl <- lvl_df$control[1]
        p <- ggplot(lvl_df, aes(x = reorder(cell_type, n_deg), y = n_deg)) +
            geom_bar(stat = "identity", fill = COLOR_ACCENT, width = 0.7) +
            geom_text(aes(label = n_deg), hjust = -0.2, size = 3) +
            coord_flip() +
            labs(title = sprintf("Sig. DEGs per cell type\n(%s vs %s, padj<0.05, |log2FC|>1)", lvl, ctrl),
                 x = NULL, y = "# DEGs") +
            theme_classic(base_size = 10) +
            FIGURE_THEME +
            expand_limits(y = max(lvl_df$n_deg) * 1.15)

        out_path <- file.path(viz_dir, paste0("07a_deg_summary__", safe_name(lvl), "_vs_", safe_name(ctrl), ".pdf"))
        save_plot(p, out_path, width_mm = 140, height_mm = 160)
        log_msg("[VIZ] Downstream summary: ", out_path)
    }
}

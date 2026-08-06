# =============================================================================
# r/stage7b_monocle3.R — Stage 7b: Trajectory analysis (Monocle3)
# =============================================================================
run_monocle3 <- function(counts, obs, ctype_col, tables_dir, viz_dir, scanpy_dir) {
    if (!requireNamespace("monocle3", quietly = TRUE)) {
        log_msg("monocle3 not available — skipping trajectory")
        return(invisible(NULL))
    }
    suppressPackageStartupMessages(library(monocle3))
    log_msg("=== Stage 7b: Trajectory (Monocle3) ===")

    umap_coords <- read_csv(
        file.path(scanpy_dir, "umap_coords.csv"),
        show_col_types = FALSE
    ) |> tibble::column_to_rownames("...1") |> as.matrix()

    gene_meta   <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
    cell_meta_m <- obs |> as.data.frame()
    rownames(cell_meta_m) <- colnames(counts)

    cds <- new_cell_data_set(
        counts,
        cell_metadata = cell_meta_m,
        gene_metadata = gene_meta
    )
    reducedDims(cds)[["UMAP"]] <- umap_coords[colnames(cds), ]

    cds <- cluster_cells(cds, reduction_method = "UMAP")
    cds <- learn_graph(cds)

    # Auto-pick root: cell with highest expression of a proliferation marker
    root_gene <- "MKI67"
    if (root_gene %in% rownames(counts)) {
        root_expr   <- counts[root_gene, ]
        root_cell   <- names(which.max(root_expr))
        root_pt     <- which(colnames(cds) == root_cell)
        cds <- order_cells(cds, root_cells = colnames(cds)[root_pt])
    } else {
        log_msg("  MKI67 not found; using interactive root (skipping auto-root)")
    }

    pseudotime_path <- file.path(viz_dir, "07b_pseudotime_umap.pdf")
    p <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
                    label_cell_groups = FALSE, graph_label_size = 3) +
        FIGURE_THEME
    save_plot(p, pseudotime_path, width_mm = 140, height_mm = 120)

    pt_path <- file.path(tables_dir, "pseudotime_values.tsv")
    pt_df <- data.frame(
        cell      = colnames(cds),
        pseudotime = pseudotime(cds),
        stringsAsFactors = FALSE
    )
    write_tsv(pt_df, pt_path)
    verify_files(c(pseudotime_path, pt_path), "Stage 7b (Monocle3)")
    log_msg("Stage 7b complete.")
}

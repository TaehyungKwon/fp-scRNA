# =============================================================================
# r/stage7e_gsea.R — Stage 7e: Gene Set Enrichment Analysis (GO:BP + KEGG +
# MSigDB Hallmark)
# =============================================================================
run_gsea <- function(all_results, contrasts_df, tables_dir, viz_dir, threads = 1, top_n = 5) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("msigdbr", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        log_msg("clusterProfiler/msigdbr/org.Hs.eg.db not available — skipping GSEA")
        return(invisible(NULL))
    }
    if (length(all_results) == 0) {
        log_msg("No DESeq2 results available — skipping GSEA")
        return(invisible(NULL))
    }
    suppressPackageStartupMessages(library(clusterProfiler))
    log_msg("=== Stage 7e: Gene Set Enrichment Analysis (GO:BP + KEGG + Hallmark) ===")

    gene_sets <- dplyr::bind_rows(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
            dplyr::mutate(collection = "Hallmark"),
        msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") |>
            dplyr::mutate(collection = "KEGG"),
        msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") |>
            dplyr::mutate(collection = "GO_BP")
    )
    t2g <- gene_sets |> dplyr::distinct(gs_name, gene_symbol)
    collection_map <- gene_sets |> dplyr::distinct(gs_name, collection)

    # Each (cell type x contrast)'s preranked GSEA test is independent and by
    # far the slowest part of Stage 7 -- run them concurrently instead of one
    # at a time. With >2 condition levels this now runs once per non-control
    # contrast per cell type (accepted ~(n_levels-1)x runtime cost).
    run_one <- function(full_key) {
        # clusterProfiler's GSEA(seed=TRUE) does set.seed(.Random.seed), which
        # requires .Random.seed to already exist -- not guaranteed in a freshly
        # forked mclapply worker after run_deseq2's own mclapply left the
        # parent unseeded. A real set.seed() call always initializes it.
        set.seed(42)
        res <- all_results[[full_key]]
        ranked <- res$stat
        names(ranked) <- res$gene
        ranked <- ranked[!is.na(ranked) & !duplicated(names(ranked))]
        ranked <- sort(ranked, decreasing = TRUE)
        if (length(ranked) < 10) return(NULL)

        gsea_res <- tryCatch(
            clusterProfiler::GSEA(ranked, TERM2GENE = t2g, pvalueCutoff = 0.05,
                                   minGSSize = 10, maxGSSize = 500, seed = TRUE, verbose = FALSE),
            error = function(e) { log_warn("  ", full_key, " — ", conditionMessage(e)); NULL }
        )
        if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
            log_msg("  ", full_key, ": no significant pathways")
            return(NULL)
        }
        df <- gsea_res@result |>
            dplyr::left_join(collection_map, by = c("ID" = "gs_name")) |>
            dplyr::arrange(p.adjust)
        row <- contrasts_df[contrasts_df$key == full_key, ]
        contrast_key <- paste0(safe_name(row$level), "_vs_", safe_name(row$control))
        write_tsv(df, file.path(tables_dir, paste0("GSEA_", safe_name(row$cell_type), "__", contrast_key, ".tsv")))
        df$cell_type <- row$cell_type
        df$level     <- row$level
        df$control   <- row$control
        log_msg("  ", full_key, ": ", nrow(df), " significant pathways (padj<0.05)")
        df
    }

    mc_cores <- max(1, min(threads, length(all_results)))
    gsea_list <- parallel::mclapply(names(all_results), run_one, mc.cores = mc_cores)
    names(gsea_list) <- names(all_results)
    all_gsea <- gsea_list[!vapply(gsea_list, is.null, logical(1))]

    gsea_tsvs <- list.files(tables_dir, pattern = "^GSEA_.*\\.tsv$", full.names = TRUE)
    log_msg("[VERIFY SUCCESS] GSEA — ", length(gsea_tsvs), " (cell type x contrast) table(s) written to ", tables_dir)

    if (length(all_gsea) == 0) {
        log_msg("[VIZ] No significant GSEA results across cell types — skipping 07e_gsea_summary")
        log_msg("Stage 7e complete.")
        return(invisible(NULL))
    }

    # One file per non-control level -- same overflow-avoidance rationale as
    # r/stage7_summary.R's make_downstream_summary()/umap_overview.pdf's
    # split: cell_type x contrast x collection stops fitting one legible
    # static page.
    all_gsea_df <- dplyr::bind_rows(all_gsea)
    for (lvl in unique(all_gsea_df$level)) {
        lvl_df <- all_gsea_df[all_gsea_df$level == lvl, ]
        ctrl <- lvl_df$control[1]
        top_terms <- lvl_df |>
            dplyr::group_by(cell_type) |>
            dplyr::slice_min(p.adjust, n = top_n, with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::mutate(Description = substr(Description, 1, 45))

        p <- ggplot(top_terms, aes(x = cell_type, y = reorder(Description, NES),
                                    size = -log10(p.adjust), color = NES)) +
            geom_point() +
            scale_color_gradient2(low = "#4C72B0", mid = "grey90", high = COLOR_ACCENT, midpoint = 0) +
            scale_size_continuous(range = c(1, 6) * 0.6) +  # 40% smaller than ggplot's own default (1,6)
            facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
            labs(title = sprintf("Top %d enriched pathways per cell type (%s vs %s)", top_n, lvl, ctrl),
                 x = NULL, y = NULL, size = "-log10(padj)", color = "NES") +
            theme_classic(base_size = 7) +
            FIGURE_THEME +
            # Overrides tuned for this figure's dense many-category layout --
            # FIGURE_THEME's 6pt axis.text would still be too large here.
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                  axis.text.y = element_text(size = 5.5),
                  strip.text.y = element_text(angle = 0, size = 7),
                  legend.box = "vertical",
                  panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.3))

        out_path <- file.path(viz_dir, paste0("07e_gsea_summary__", safe_name(lvl), "_vs_", safe_name(ctrl), ".pdf"))
        save_plot(p, out_path, width_mm = PAGE_W_MM, height_mm = PAGE_H_MM)
        log_msg("[VIZ] GSEA summary: ", out_path)
    }
    log_msg("Stage 7e complete.")
}

# =============================================================================
# r/stage7a_deseq2.R — Stage 7a: Pseudobulk DEG (DESeq2)
# =============================================================================
run_deseq2 <- function(counts, obs, cond_col, ctype_col, control_condition, tables_dir, viz_dir, threads = 1) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        log_msg("DESeq2 not available — skipping DEG analysis")
        return(invisible(list(results = list(), contrasts = NULL)))
    }
    suppressPackageStartupMessages(library(DESeq2))
    log_msg("=== Stage 7a: Pseudobulk DEG (DESeq2) ===")

    cell_types <- unique(obs[[ctype_col]])

    # Each cell type's DESeq2 fit is independent -- run them concurrently
    # (mc.cores workers, fork-based, Linux-only) instead of one at a time.
    # Returns a named list of per-contrast result data frames (one per
    # non-control level), or NULL if this cell type has no valid contrast.
    run_one <- function(ct) {
        ct_cells  <- which(obs[[ctype_col]] == ct)
        ct_obs    <- obs[ct_cells, ]
        ct_counts <- counts[, ct_cells, drop = FALSE]

        # Pseudobulk: aggregate counts per sample
        samples   <- unique(ct_obs$sample)
        if (length(samples) < 2) return(NULL)

        pb_counts <- vapply(samples, function(s) {
            idx <- which(ct_obs$sample == s)
            if (length(idx) == 0) return(rep(0L, nrow(ct_counts)))
            Matrix::rowSums(ct_counts[, idx, drop = FALSE])
        }, numeric(nrow(ct_counts)))
        colnames(pb_counts) <- samples

        # Sample-level metadata (one row per sample)
        pb_meta <- ct_obs |>
            dplyr::distinct(sample, .keep_all = TRUE) |>
            dplyr::filter(sample %in% samples) |>
            as.data.frame()
        rownames(pb_meta) <- pb_meta$sample
        pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]

        ct_levels <- as.character(unique(pb_meta[[cond_col]]))  # see condition_levels' comment in downstream.R
        if (length(ct_levels) < 2) return(NULL)

        # Relevel so the designated control is the DESeq2 reference. If this
        # cell type's own samples don't happen to include the control level
        # (e.g. a cell type entirely absent from the control sample), there's
        # no valid vs-control contrast to run for it -- skip rather than
        # silently falling back to a different reference.
        if (control_condition != "") {
            if (!(control_condition %in% ct_levels)) return(NULL)
            ref_level <- control_condition
        } else {
            ref_level <- ct_levels[1]  # 2-level, unset case -- unchanged prior behavior
        }
        pb_meta[[cond_col]] <- stats::relevel(factor(pb_meta[[cond_col]]), ref = ref_level)
        test_levels <- setdiff(levels(pb_meta[[cond_col]]), ref_level)

        dds <- tryCatch({
            dds <- DESeqDataSetFromMatrix(
                countData = round(pb_counts),
                colData   = pb_meta,
                design    = as.formula(paste("~", cond_col))
            )
            DESeq(dds, quiet = TRUE)
        }, error = function(e) { log_warn("  ", ct, " — ", conditionMessage(e)); NULL })
        if (is.null(dds)) return(NULL)

        out <- list()
        for (lvl in test_levels) {
            res <- tryCatch({
                results(dds, contrast = c(cond_col, lvl, ref_level), alpha = 0.05) |>
                    as.data.frame() |>
                    tibble::rownames_to_column("gene") |>
                    dplyr::filter(!is.na(padj)) |>
                    dplyr::arrange(padj)
            }, error = function(e) {
                log_warn("  ", ct, " (", lvl, " vs ", ref_level, ") — ", conditionMessage(e)); NULL
            })
            if (is.null(res)) next
            contrast_key <- paste0(safe_name(lvl), "_vs_", safe_name(ref_level))
            write_tsv(res, file.path(tables_dir, paste0("DEG_", safe_name(ct), "__", contrast_key, ".tsv")))
            log_msg("  ", ct, " (", lvl, " vs ", ref_level, "): ",
                    sum(res$padj < 0.05, na.rm = TRUE), " sig. DEGs")
            out[[contrast_key]] <- res
        }
        if (length(out) == 0) return(NULL)
        out
    }

    mc_cores <- max(1, min(threads, length(cell_types)))
    results_list <- parallel::mclapply(cell_types, run_one, mc.cores = mc_cores)
    names(results_list) <- cell_types
    results_list <- results_list[!vapply(results_list, is.null, logical(1))]

    # Flatten to a single list keyed by "<cell_type>||<contrast_key>", plus a
    # parallel contrasts data frame (cell_type, level, control, key) that
    # run_gsea()/make_downstream_summary()/run_milo() group and facet by.
    all_results <- list()
    contrasts_df <- data.frame(cell_type = character(), level = character(),
                               control = character(), key = character(),
                               stringsAsFactors = FALSE)
    for (ct in names(results_list)) {
        for (contrast_key in names(results_list[[ct]])) {
            full_key <- paste(ct, contrast_key, sep = "||")
            all_results[[full_key]] <- results_list[[ct]][[contrast_key]]
            parts <- strsplit(contrast_key, "_vs_")[[1]]
            contrasts_df <- rbind(contrasts_df, data.frame(
                cell_type = ct, level = parts[1], control = parts[2], key = full_key,
                stringsAsFactors = FALSE))
        }
    }

    # Volcano plot: single global pick (best cell_type x contrast pair by n DEGs).
    volcano_path <- file.path(viz_dir, "07a_volcano_top_celltype.pdf")
    if (length(all_results) > 0) {
        best_key <- names(which.max(sapply(all_results, function(r) sum(r$padj < 0.05, na.rm=TRUE))))
        best_row <- contrasts_df[contrasts_df$key == best_key, ]
        df <- all_results[[best_key]] |>
            dplyr::mutate(
                sig   = padj < 0.05 & abs(log2FoldChange) > 1,
                label = ifelse(sig & rank(padj) <= 15, gene, "")
            )
        p <- ggplot(df, aes(log2FoldChange, -log10(pvalue), color = sig)) +
            geom_point(alpha = 0.5, size = 0.8) +
            geom_text(aes(label = label), size = 2.5, vjust = -0.7,
                     show.legend = FALSE, check_overlap = TRUE) +
            scale_color_manual(values = c(COLOR_MUTED, COLOR_ACCENT)) +
            theme_classic(base_size = 12) +
            labs(title = sprintf("Volcano: %s (%s vs %s)", best_row$cell_type, best_row$level, best_row$control),
                 x = "log2 FC", y = "-log10 p") +
            FIGURE_THEME +
            theme(legend.position = "none")  # no legend here; overrides FIGURE_THEME's "bottom"
        save_plot(p, volcano_path, width_mm = 140, height_mm = 120)
    }

    deg_tsvs <- list.files(tables_dir, pattern = "^DEG_.*\\.tsv$", full.names = TRUE)
    log_msg("[VERIFY SUCCESS] DESeq2 — ", length(deg_tsvs), " (cell type x contrast) table(s) written to ", tables_dir)
    verify_files(volcano_path, "Stage 7a (DESeq2 volcano)")
    log_msg("Stage 7a complete.")
    invisible(list(results = all_results, contrasts = contrasts_df))
}

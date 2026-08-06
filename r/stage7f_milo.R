# =============================================================================
# r/stage7f_milo.R — Stage 7f: Differential abundance via kNN neighborhoods
# (Milo). Alternative differential-abundance method to Stage 6's scCODA --
# the two don't run adjacently since scCODA only needs obs_metadata.csv
# (runs right after Stage 5) while Milo needs PCA coordinates from the h5ad
# (runs here, inside Stage 7).
# =============================================================================
# Needs PCA coordinates -- only available via the annotated.h5ad/zellkonverter
# path (the counts_matrix.csv fallback near the top of downstream.R never
# builds an `sce` object at all). Callers must check `exists("sce")`
# themselves before calling this -- `sce` can't be a normal optional
# parameter here, since passing an unbound top-level variable as an argument
# raises before the function body ever runs, R's lazy evaluation notwithstanding.
run_milo <- function(sce, obs, cond_col, ctype_col, control_condition, tables_dir, viz_dir) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
        log_msg("miloR not available — skipping Milo differential abundance")
        return(invisible(NULL))
    }
    suppressPackageStartupMessages({
        library(miloR)
        library(SingleCellExperiment)
    })
    log_msg("=== Stage 7f: Differential abundance (Milo) ===")

    # Mirrors py/stage4_embed.py's own use_rep fallback: prefer the
    # Harmony-corrected embedding when present, since both land in the h5ad
    # via adata.obsm without any extra Python-side export.
    reduced_dim <- if ("X_pca_harmony" %in% reducedDimNames(sce)) "X_pca_harmony" else "X_pca"
    if (!(reduced_dim %in% reducedDimNames(sce))) {
        log_msg("Neither X_pca_harmony nor X_pca found in annotated.h5ad — skipping Milo")
        return(invisible(NULL))
    }
    log_msg("  Using reduced dim: ", reduced_dim)

    milo_obj <- Milo(sce)
    milo_obj <- buildGraph(milo_obj, k = 30, d = 30, reduced.dim = reduced_dim)
    milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = 30, refined = TRUE,
                           reduced_dims = reduced_dim)
    milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), samples = "sample")

    sample_meta <- obs |>
        dplyr::distinct(sample, .keep_all = TRUE) |>
        as.data.frame()
    rownames(sample_meta) <- sample_meta$sample
    sample_meta <- sample_meta[colnames(nhoodCounts(milo_obj)), , drop = FALSE]

    ct_levels <- as.character(unique(sample_meta[[cond_col]]))  # see condition_levels' comment in downstream.R
    ref_level <- if (control_condition != "" && control_condition %in% ct_levels) {
        control_condition
    } else {
        ct_levels[1]  # 2-level, unset case -- same arbitrary-but-deterministic reference as DESeq2
    }
    sample_meta[[cond_col]] <- stats::relevel(factor(sample_meta[[cond_col]]), ref = ref_level)
    test_levels <- setdiff(levels(sample_meta[[cond_col]]), ref_level)

    milo_obj <- calcNhoodDistance(milo_obj, d = 30, reduced.dim = reduced_dim)

    n_written <- 0
    for (lvl in test_levels) {
        contrast_str <- paste0(cond_col, lvl)
        da_res <- tryCatch({
            testNhoods(milo_obj, design = as.formula(paste("~", cond_col)),
                      design.df = sample_meta, reduced.dim = reduced_dim,
                      model.contrasts = contrast_str)
        }, error = function(e) {
            log_warn("  ", lvl, " vs ", ref_level, " — ", conditionMessage(e)); NULL
        })
        if (is.null(da_res)) next

        da_res <- annotateNhoods(milo_obj, da_res, coldata_col = ctype_col)
        contrast_key <- paste0(safe_name(lvl), "_vs_", safe_name(ref_level))
        out_path <- file.path(tables_dir, paste0("milo_DA_", contrast_key, ".tsv"))
        write_tsv(da_res, out_path)
        n_written <- n_written + 1
        n_sig <- sum(da_res$SpatialFDR < 0.1, na.rm = TRUE)
        log_msg("  ", lvl, " vs ", ref_level, ": ", n_sig, " DA neighborhoods (SpatialFDR<0.1)")

        if (n_sig == 0) {
            # plotDAbeeswarm() only colors neighborhoods with SpatialFDR <
            # alpha -- every other point gets logFC_color = NA, which ggplot2
            # renders invisibly. With 0 significant neighborhoods, that isn't
            # a smaller plot, it's a *blank* one: every point is transparent.
            # Same "skip and log why" rule run_gsea() already uses for 0
            # significant pathways, rather than shipping a technically-present
            # but uninformative PDF.
            log_msg("  [VIZ] 0 significant neighborhoods for ", contrast_key,
                    " — skipping beeswarm plot (would render blank)")
        } else {
            # ggplot objects build lazily -- ggplot2's scale/stat computation
            # (e.g. plotDAbeeswarm's color scale) doesn't actually run until
            # the plot is drawn/rendered, which happens inside save_plot()'s
            # ggsave() call, not at construction. The tryCatch must wrap both
            # steps together, or a rendering-time error escapes uncaught and
            # aborts the whole script instead of just skipping this one plot.
            tryCatch({
                p <- plotDAbeeswarm(da_res, group.by = ctype_col) + FIGURE_THEME
                beeswarm_path <- file.path(viz_dir, paste0("07f_milo_beeswarm__", contrast_key, ".pdf"))
                save_plot(p, beeswarm_path, width_mm = PAGE_W_MM, height_mm = PAGE_H_MM * 0.6)
            }, error = function(e) { log_warn("  beeswarm plot (", contrast_key, ") — ", conditionMessage(e)) })
        }
    }

    if (n_written == 0) {
        log_msg("[VIZ] No Milo DA results across any contrast — likely no residual degrees of ",
                "freedom for dispersion estimation (needs >=2 samples per condition)")
    } else {
        milo_tsvs <- list.files(tables_dir, pattern = "^milo_DA_.*\\.tsv$", full.names = TRUE)
        log_msg("[VERIFY SUCCESS] Milo — ", length(milo_tsvs), " contrast table(s) written to ", tables_dir)
    }
    log_msg("Stage 7f complete.")
}

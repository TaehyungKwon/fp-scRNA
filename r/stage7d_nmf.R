# =============================================================================
# r/stage7d_nmf.R — Stage 7d: Gene programs via NMF
# =============================================================================
run_nmf <- function(counts, tables_dir, k = 8) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
        log_msg("NMF not available — skipping gene program discovery")
        return(invisible(NULL))
    }
    suppressPackageStartupMessages(library(NMF))
    log_msg(sprintf("=== Stage 7d: Gene programs (NMF, k=%d) ===", k))

    # Subsample genes: top 2000 by variance (fast and sufficient)
    gene_var  <- apply(counts, 1, var)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
    mat       <- as.matrix(counts[top_genes, ])
    mat[mat < 0] <- 0  # NMF requires non-negative input

    res <- nmf(mat, rank = k, method = "snmf/r", nrun = 5, seed = 42, .options = "p4v")

    # W: gene × program loadings
    W <- basis(res)
    colnames(W) <- paste0("Program_", seq_len(k))
    loadings_path <- file.path(tables_dir, "nmf_gene_loadings.tsv")
    write_tsv(as.data.frame(W) |> tibble::rownames_to_column("gene"), loadings_path)

    # H: program × cell scores
    H <- coef(res)
    rownames(H) <- paste0("Program_", seq_len(k))
    scores_path <- file.path(tables_dir, "nmf_cell_scores.tsv")
    write_tsv(as.data.frame(t(H)) |> tibble::rownames_to_column("cell"), scores_path)

    # Top genes per program
    top_per_program <- apply(W, 2, function(w) {
        names(sort(w, decreasing = TRUE))[1:20]
    })
    top_path <- file.path(tables_dir, "nmf_top_genes_per_program.tsv")
    write_tsv(as.data.frame(top_per_program), top_path)

    verify_files(c(loadings_path, scores_path, top_path), "Stage 7d (NMF)")
    log_msg("Stage 7d complete. Programs saved: ", tables_dir)
}

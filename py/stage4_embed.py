"""
py/stage4_embed.py — Stage 4: Normalize & embed.
4a: normalize (normalization & HVG selection). 4b: embed (PCA + optional
Harmony batch correction + kNN + UMAP).
"""

import logging

import scanpy as sc
import anndata as ad

log = logging.getLogger(__name__)


# ── Stage 4a: Normalization & HVG selection ───────────────────────────────────
def normalize(adata: ad.AnnData, n_hvgs: int) -> ad.AnnData:
    log.info("Normalizing and selecting HVGs...")
    # Store raw counts before normalization
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # freeze log-normalized matrix for DEG later

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvgs,
        batch_key="sample",   # accounts for batch when selecting HVGs
        flavor="seurat_v3",
        layer="counts",
    )
    log.info(f"  HVGs selected: {adata.var['highly_variable'].sum()}")
    return adata


# ── Stage 4b: PCA + optional batch correction + kNN + UMAP ───────────────────
def embed(adata: ad.AnnData, n_pcs: int, n_neighbors: int) -> ad.AnnData:
    log.info("Scaling, PCA, batch correction, kNN, UMAP...")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack")

    use_rep = "X_pca"

    # Harmony batch correction: run when there are multiple samples to correct.
    # "sample" is always the batch key — it's set during load and is always present.
    n_samples = adata.obs["sample"].nunique()
    if n_samples > 1:
        try:
            import harmonypy as hm
            log.info(f"  Running Harmony batch correction ({n_samples} samples)...")
            ho = hm.run_harmony(
                adata.obsm["X_pca"],
                adata.obs,
                "sample",
                max_iter_harmony=20,
                verbose=False,
            )
            # harmonypy>=2.0's C++ backend returns Z_corr as (n_cells, n_pcs)
            # already -- the pre-2.0 pure-Python backend returned (n_pcs, n_cells).
            adata.obsm["X_pca_harmony"] = ho.Z_corr
            use_rep = "X_pca_harmony"
            log.info("  Harmony complete.")
        except ImportError:
            log.warning("  harmonypy not found — skipping batch correction (using raw PCA)")
    else:
        log.info("  Single sample detected — skipping batch correction.")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    sc.tl.umap(adata)
    log.info("  UMAP done.")
    return adata

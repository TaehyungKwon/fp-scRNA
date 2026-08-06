"""
py/stage5_cluster.py — Stage 5: Cluster & annotate.
5a: Leiden clustering + marker genes. 5b: cell type annotation (CellTypist,
optional; falls back to leiden cluster IDs).
"""

import logging
from pathlib import Path

import scanpy as sc
import anndata as ad

log = logging.getLogger(__name__)


def cluster_and_annotate(adata: ad.AnnData, resolution: float,
                          out_dir: Path) -> ad.AnnData:
    # ── Stage 5a: Leiden clustering + marker genes ────────────────────────────
    log.info(f"Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    log.info(f"  Clusters found: {n_clusters}")

    # Marker genes per cluster
    log.info("  Computing marker genes (Wilcoxon)...")
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=True,
        pts=True,           # fraction of cells expressing gene
    )

    # Save top markers to CSV
    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    markers_df.to_csv(out_dir / "cluster_markers.csv", index=False)
    log.info(f"  Markers saved: {out_dir / 'cluster_markers.csv'}")

    # ── Stage 5b: Cell type annotation (CellTypist, optional) ────────────────
    try:
        import celltypist
        from celltypist import models
        log.info("  Running CellTypist annotation...")
        model = models.Model.load(model="Immune_All_Low.pkl")
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        adata = predictions.to_adata()
        adata.obs["cell_type"] = adata.obs["majority_voting"]
        log.info("  CellTypist annotation complete.")
    except ImportError:
        log.warning("  celltypist not found — cell_type set to leiden cluster IDs")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)

    return adata

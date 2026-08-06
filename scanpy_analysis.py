#!/usr/bin/env python3
"""
scanpy_analysis.py — Stage 3 (Load & QC-filter) -> Stage 4 (Normalize &
embed) -> Stage 5 (Cluster & annotate). Entry point: parses CLI args
(py/common.py's parse_args()) and orchestrates the py/stage3_qc.py ->
py/stage4_embed.py -> py/stage5_cluster.py -> py/reports.py functions below,
in order. Minimal external dependencies: scanpy, anndata, scrublet,
harmonypy, matplotlib -- all heavy lifting done by scanpy's built-ins;
optional harmonypy for batch correction.
"""

import sys
import warnings
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc

from py.common import parse_args, verify_anndata, verify_outputs
from py.stage3_qc import load_samples, merge_metadata, run_qc
from py.stage4_embed import normalize, embed
from py.stage5_cluster import cluster_and_annotate
from py.reports import (
    make_cellranger_summary, make_qc_summary, make_qc_metrics_boxplot,
    make_clustering_summary, make_doublet_summary, make_plots,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

log = logging.getLogger(__name__)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    viz_dir = Path(args.viz_dir) if args.viz_dir else out_dir.parent / "summary_figures"
    viz_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Summary figures directory: {viz_dir}")

    tables_dir = Path(args.tables_dir) if args.tables_dir else out_dir.parent / "summary_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Summary tables directory: {tables_dir}")

    sc.settings.verbosity = 1
    sc.settings.n_jobs = args.threads

    # -- Stage 3a: Load & merge metadata ---------------------------------------
    adata = load_samples(args.matrix_dir, args.samples, args.input_mode, args.preprocess_qc_dir)
    verify_anndata(adata, "load_samples")
    cells_before_qc = adata.obs["sample"].value_counts().to_dict()

    adata = merge_metadata(adata, args.metadata_file)

    # -- Stage 3b: QC metrics & filtering ---------------------------------------
    adata, qc_metrics_before = run_qc(adata, args.min_genes, args.max_genes, args.max_mito,
                                       args.n_mads, args.n_mads_mt, args.remove_doublets)
    verify_anndata(adata, "QC filtering")

    # -- Stage 4: Normalize & embed ---------------------------------------------
    adata = normalize(adata, args.n_hvgs)
    adata = embed(adata, args.n_pcs, args.n_neighbors)

    # -- Stage 5: Cluster & annotate ---------------------------------------------
    adata = cluster_and_annotate(adata, args.resolution, out_dir)
    make_plots(adata, viz_dir)

    h5ad_path = out_dir / "annotated.h5ad"
    adata.write_h5ad(h5ad_path)
    log.info(f"AnnData saved: {h5ad_path}")

    # Also export obs table for R (cross-language data contract -- keep as CSV
    # in scanpy-dir; downstream.R reads it from exactly this path).
    obs_path = out_dir / "obs_metadata.csv"
    adata.obs.to_csv(obs_path)
    log.info(f"Metadata exported: {obs_path}")

    # Export UMAP coords for R (same contract as above)
    umap_df = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names
    )
    umap_df.to_csv(out_dir / "umap_coords.csv")

    verify_outputs(adata, out_dir)

    # -- Summary tables (human-facing TSV copies) ------------------------------
    log.info("Writing summary tables...")
    adata.obs.to_csv(tables_dir / "obs_metadata.tsv", sep="\t")
    umap_df.to_csv(tables_dir / "umap_coords.tsv", sep="\t")
    pd.read_csv(out_dir / "cluster_markers.csv").to_csv(
        tables_dir / "cluster_markers.tsv", sep="\t", index=False
    )

    # -- Summary visualizations -----------------------------------------------
    log.info("Generating summary visualizations...")
    if args.input_mode == "fastq":
        make_cellranger_summary(args.matrix_dir, args.samples, viz_dir, tables_dir)
    else:
        log.info("[VIZ] input-mode=matrix: no metrics_summary.csv available, "
                  "skipping 01_cellranger_summary")
    make_qc_summary(adata, cells_before_qc, viz_dir)
    make_qc_metrics_boxplot(qc_metrics_before, adata, viz_dir)
    make_clustering_summary(adata, viz_dir)
    make_doublet_summary(adata, viz_dir, tables_dir)
    log.info("Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # Otherwise this prints an untagged, unformatted Python traceback.
        log.exception(str(e))
        sys.exit(1)

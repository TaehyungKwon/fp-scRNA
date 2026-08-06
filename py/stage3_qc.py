"""
py/stage3_qc.py — Stage 3: Load & QC-filter.
3a: load_samples/merge_metadata (load & concatenate samples, merge sample
metadata). 3b: run_qc (QC metrics & MAD-based filtering, + scrublet fallback
doublet detection if Stage 2b's scDblFinder didn't run).
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import median_abs_deviation

from py.common import QC_METRICS

log = logging.getLogger(__name__)


def is_outlier(adata: ad.AnnData, metric: str, nmads: float) -> pd.Series:
    """MAD-based outlier flag: True where `metric` is more than `nmads` median
    absolute deviations from the median. Follows the sc-best-practices QC
    notebook's approach -- more permissive and less biased against small
    subpopulations than a single fixed min/max threshold."""
    m = adata.obs[metric]
    mad = median_abs_deviation(m)
    med = np.median(m)
    return (m < med - nmads * mad) | (med + nmads * mad < m)


# ── Stage 3a: Load & concatenate samples ─────────────────────────────────────
def load_samples(matrix_dir: str, samples: list[str], input_mode: str = "fastq",
                  preprocess_qc_dir: str | None = None) -> ad.AnnData:
    """
    input_mode="fastq":  <matrix_dir>/<sample>/outs/filtered_feature_bc_matrix/,
                         or .../outs/soupx_corrected_matrix/ if preprocess_qc.R's
                         ambient RNA correction ran for that sample (preferred
                         when present); .../outs/doublet_calls.tsv (scDblFinder,
                         also from preprocess_qc.R) is merged in when present.
    input_mode="matrix": <matrix_dir>/<sample>/  (flat -- matrix.mtx.gz,
                         features.tsv.gz, barcodes.tsv.gz directly inside) is
                         always the count matrix (SoupX never runs in this mode
                         -- no raw/unfiltered matrix available for its ambient
                         profile). scDblFinder still runs per-sample in this
                         mode (preprocess_qc.R), writing doublet_calls.tsv to
                         <preprocess_qc_dir>/<sample>/ -- merged in when present.
    """
    log.info(f"Loading feature-barcode matrices (input mode: {input_mode})...")
    adatas = {}
    for sample in samples:
        doublet_path = None
        if input_mode == "matrix":
            matrix_path = Path(matrix_dir) / sample
            if preprocess_qc_dir:
                doublet_path = Path(preprocess_qc_dir) / sample / "doublet_calls.tsv"
        else:
            sample_outs = Path(matrix_dir) / sample / "outs"
            soupx_path = sample_outs / "soupx_corrected_matrix"
            if soupx_path.exists():
                matrix_path = soupx_path
                doublet_path = sample_outs / "doublet_calls.tsv"
            else:
                matrix_path = sample_outs / "filtered_feature_bc_matrix"
        if not matrix_path.exists():
            raise FileNotFoundError(f"Matrix not found: {matrix_path}")
        # cache=False: scanpy's read-cache isn't invalidated when Cell Ranger
        # reruns at the same path, so it can silently serve stale data.
        adata = sc.read_10x_mtx(matrix_path, var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata.obs["sample"] = sample

        if doublet_path is not None and doublet_path.exists():
            doublets = pd.read_csv(doublet_path, sep="\t", index_col="barcode")
            adata.obs["scDblFinder_score"] = adata.obs_names.map(doublets["scDblFinder_score"])
            adata.obs["scDblFinder_class"] = adata.obs_names.map(doublets["scDblFinder_class"])
            n_doublets = (adata.obs["scDblFinder_class"] == "doublet").sum()
            soupx_note = "SoupX-corrected, " if matrix_path.name == "soupx_corrected_matrix" else ""
            log.info(f"  {sample}: {adata.n_obs} cells × {adata.n_vars} genes "
                     f"({soupx_note}{n_doublets} scDblFinder doublets)")
        else:
            log.info(f"  {sample}: {adata.n_obs} cells × {adata.n_vars} genes")
        adatas[sample] = adata

    combined = ad.concat(adatas, label="sample", merge="same")
    combined.obs_names_make_unique()
    log.info(f"Combined: {combined.n_obs} cells × {combined.n_vars} genes")
    return combined


def merge_metadata(adata: ad.AnnData, metadata_file: str) -> ad.AnnData:
    """
    Reads sample_metadata.tsv and left-joins it onto adata.obs by the 'sample'
    column.  Every column in the TSV (condition, sex, batch, timepoint, …)
    becomes available as a per-cell obs column automatically — no hardcoding.
    """
    log.info(f"Merging sample metadata from: {metadata_file}")
    meta = pd.read_csv(metadata_file, sep="\t", dtype=str)

    if "sample" not in meta.columns:
        raise ValueError("sample_metadata.tsv must contain a 'sample' column.")

    missing = set(adata.obs["sample"].unique()) - set(meta["sample"])
    if missing:
        raise ValueError(f"Samples in data but missing from metadata TSV: {missing}")

    meta = meta.set_index("sample")
    # Map each cell's sample to the metadata columns
    for col in meta.columns:
        adata.obs[col] = adata.obs["sample"].map(meta[col]).values
        log.info(f"  Added obs column: '{col}'")

    return adata


# ── Stage 3b: QC metrics & filtering ─────────────────────────────────────────
def run_qc(adata: ad.AnnData, min_genes: int, max_genes: int, max_mito: float,
           n_mads: float, n_mads_mt: float, remove_doublets: bool = False) -> tuple[ad.AnnData, pd.DataFrame]:
    log.info("Computing QC metrics...")

    # Annotate mitochondrial, ribosomal, and hemoglobin genes
    adata.var["mt"]   = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"]   = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],   # needed for pct_counts_in_top_20_genes, an outlier covariate below
        log1p=True,         # needed for log1p_total_counts / log1p_n_genes_by_counts below
        inplace=True,
    )
    qc_metrics_before = adata.obs[QC_METRICS].copy()  # snapshot for the before/after plot

    n_before = adata.n_obs
    log.info(f"  Cells before filtering: {n_before}")

    # ── Coarse absolute floor/ceiling ─────────────────────────────────────────
    # Removes obviously-broken barcodes before they can skew the median/MAD
    # statistics the outlier detection below relies on.
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)

    # ── MAD-based outlier detection (sc-best-practices QC notebook) ──────────
    # More permissive and less biased against small subpopulations than a
    # single fixed threshold -- see CLAUDE.md's "Filtering strategy" section.
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", n_mads)
        | is_outlier(adata, "log1p_n_genes_by_counts", n_mads)
        | is_outlier(adata, "pct_counts_in_top_20_genes", n_mads)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", n_mads_mt) | (adata.obs["pct_counts_mt"] > max_mito)
    log.info(f"  MAD outliers: {adata.obs['outlier'].sum()}, "
             f"mt outliers: {adata.obs['mt_outlier'].sum()} (>{max_mito}% or {n_mads_mt} MADs)")
    adata = adata[~(adata.obs["outlier"] | adata.obs["mt_outlier"])].copy()

    # ── Doublet detection ─────────────────────────────────────────────────────
    # Prefer scDblFinder calls from preprocess_qc.R (Stage 2b, run per-sample,
    # before combination -- doublet detection must never see pooled
    # multi-batch data); fall back to scrublet if that stage wasn't run.
    # Either way, doublets are flagged, not removed -- inspect
    # `predicted_doublet` during visualization/clustering rather than
    # dropping cells blind at this stage.
    if "scDblFinder_class" in adata.obs.columns:
        adata.obs["predicted_doublet"] = adata.obs["scDblFinder_class"] == "doublet"
        adata.obs["doublet_score"] = adata.obs["scDblFinder_score"]
        log.info(f"  Doublets flagged (scDblFinder): {adata.obs['predicted_doublet'].sum()}")
    else:
        try:
            import scrublet as scr
            scrub = scr.Scrublet(adata.X)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
            adata.obs["doublet_score"]     = doublet_scores
            adata.obs["predicted_doublet"] = predicted_doublets
            log.info(f"  Doublets flagged (scrublet fallback): {adata.obs['predicted_doublet'].sum()}")
        except ImportError:
            log.warning("  scrublet not found and no scDblFinder calls present — skipping doublet detection")
            adata.obs["predicted_doublet"] = False

    # Opt-in only (--remove-doublets / REMOVE_DOUBLETS=true) -- default keeps
    # doublets flagged-not-removed, per this pipeline's documented design.
    # Inspect 03b_doublet_summary.pdf before enabling this on a real dataset.
    if remove_doublets:
        n_doublets = int(adata.obs["predicted_doublet"].sum())
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        log.info(f"  Doublets removed (REMOVE_DOUBLETS=true): {n_doublets}")

    # ── Gene filter: keep genes expressed in ≥3 cells ─────────────────────────
    sc.pp.filter_genes(adata, min_cells=3)

    log.info(f"  Cells after filtering: {adata.n_obs} (removed {n_before - adata.n_obs})")
    log.info(f"  Genes retained: {adata.n_vars}")
    return adata, qc_metrics_before

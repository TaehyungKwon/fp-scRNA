"""
py/common.py — shared infrastructure for scanpy_analysis.py's Stage 3-5
submodules (py/stage3_qc.py, py/stage4_embed.py, py/stage5_cluster.py,
py/reports.py): logging setup, the shared visual style/page-envelope
constants (mirrored in differential_abundance.py and downstream.R so every
figure across all three languages reads as one system), generic plotting
helpers, cross-language contract verification, and the CLI.
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe on headless servers
import matplotlib.pyplot as plt

# ── Logging ───────────────────────────────────────────────────────────────────
# %(levelname)s is exactly "INFO"/"WARNING"/"ERROR" -- matches the
# [INFO]/[WARNING]/[ERROR] tag convention shared with run_pipeline.sh (bash)
# and downstream.R. Configured once here; every py/*.py submodule's own
# logging.getLogger(__name__) inherits this format.
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s]\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ── Shared visual style ────────────────────────────────────────────────────────
# Fixed hex codes (not a matplotlib colormap name) so downstream.R and
# differential_abundance.py can mirror the exact same values -- figures
# across all three languages then read as one system even though they're
# rendered by different plotting libraries.
PALETTE_QUALITATIVE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
]
COLOR_BEFORE = "#4C72B0"   # "before" / reference state
COLOR_AFTER  = "#55A868"   # "after" / good state
COLOR_WARN   = "#C44E52"   # low / bad
COLOR_MID    = "#DD8452"   # medium
COLOR_ACCENT = "#D85A30"   # highlight (DEGs, significant points)
COLOR_MUTED  = "0.7"       # de-emphasized background points -- matplotlib's grayscale-float
                            # string syntax; "0.7" is the equivalent of R/ggplot2's "grey70"
                            # (both mean 70% white), which isn't a valid matplotlib color name

FIGURE_DPI = 300

# Target page envelope for every summary figure: 180mm x 240mm (fits a
# double-column journal page). mm -> inch for matplotlib's figsize.
MM_PER_IN = 25.4
PAGE_W_IN = 180 / MM_PER_IN
PAGE_H_IN = 240 / MM_PER_IN

# Real Arial isn't installed on Linux; Liberation Sans is the metrically-
# compatible open substitute (what LibreOffice/most Linux distros ship as
# the Arial equivalent) -- request it directly rather than "Arial" first,
# which just logs a "not found" warning for every text element before
# falling back anyway. Sizes below match downstream.R's FIGURE_THEME.
plt.rcParams["font.family"]      = ["Liberation Sans", "sans-serif"]
plt.rcParams["axes.labelsize"]   = 6
plt.rcParams["xtick.labelsize"]  = 6
plt.rcParams["ytick.labelsize"]  = 6
plt.rcParams["legend.fontsize"]  = 5
AXIS_FONTSIZE   = 6
LEGEND_FONTSIZE = 5

# QC metrics tracked before/after filtering (see stage3_qc.py's run_qc() /
# reports.py's make_qc_metrics_boxplot()).
QC_METRICS = ["n_genes_by_counts", "total_counts", "pct_counts_mt",
              "pct_counts_ribo", "pct_counts_hb", "pct_counts_in_top_20_genes"]


def category_colors(categories) -> dict:
    """Stable name -> color mapping so a given category (e.g. a sample name)
    gets the same color in every figure, regardless of that figure's own
    plotting order."""
    cats = sorted(set(categories))
    return {c: PALETTE_QUALITATIVE[i % len(PALETTE_QUALITATIVE)] for i, c in enumerate(cats)}


def legend_below(ax, ncol: int | None = None, y_offset: float = -0.38) -> None:
    """Move an axis's legend below the plot instead of overlapping the data.
    `y_offset` defaults low enough to clear the 35 deg-rotated x-tick labels
    used throughout this module's per-sample bar charts."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    ncol = ncol or min(len(labels), 4)
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y_offset),
               ncol=ncol, fontsize=LEGEND_FONTSIZE, framealpha=0.9, frameon=False,
               handletextpad=0.3, columnspacing=0.8)


def grid_figure(n_panels: int, ncols: int = 2, row_h_in: float = 3.2, hspace: float = 1.0):
    """Figure with an `ncols`-column grid of axes sized to fit the shared
    180x240mm page envelope, with any unused trailing axes hidden.
    `hspace` is generous by default -- rows must reserve enough vertical
    gap for below-axis legends, or a top row's legend overlaps the row
    beneath it (plt.subplots packs rows tightly otherwise)."""
    nrows = -(-n_panels // ncols)  # ceil division
    h = min(row_h_in * nrows, PAGE_H_IN)
    fig, axes = plt.subplots(nrows, ncols, figsize=(PAGE_W_IN, h), squeeze=False,
                              gridspec_kw={"hspace": hspace, "wspace": 0.35})
    axes = axes.flatten()
    for ax in axes[n_panels:]:
        ax.axis("off")
    return fig, axes[:n_panels]


def save_fig(fig, out_path: Path, dpi: int = FIGURE_DPI) -> Path:
    """Save a matplotlib figure as PDF at `dpi` (affects rasterized elements
    like scatter points). `out_path`'s suffix is ignored. `bbox_inches="tight"`
    expands the exported canvas to include artists like below-axis legends
    that extend past the nominal figsize, so this also checks the final
    render against the shared 180x240mm page envelope and warns if exceeded."""
    out = out_path.with_suffix(".pdf")
    tight_bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    w_mm, h_mm = tight_bbox.width * MM_PER_IN, tight_bbox.height * MM_PER_IN
    if w_mm > 180 + 1 or h_mm > 240 + 1:
        log.warning(f"[VIZ] {out.name}: rendered size {w_mm:.0f}x{h_mm:.0f}mm "
                    f"exceeds the 180x240mm page envelope")
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    return out


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="scRNA-seq: Stage 3 (QC) -> Stage 4 (Normalize/Embed) -> Stage 5 (Cluster/Annotate)")
    p.add_argument("--matrix-dir",     required=True,
                   help="Root dir of feature-barcode matrices. Layout depends on "
                        "--input-mode: 'fastq' expects <matrix-dir>/<sample>/outs/"
                        "filtered_feature_bc_matrix/ (Cell Ranger's own layout); "
                        "'matrix' expects the flat <matrix-dir>/<sample>/ containing "
                        "matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz directly.")
    p.add_argument("--input-mode",     choices=["fastq", "matrix"], default="fastq",
                   help="'fastq': matrices were produced by this pipeline's own Cell "
                        "Ranger stage. 'matrix': pre-computed matrices from elsewhere "
                        "(skip Cell Ranger's outs/ nesting).")
    p.add_argument("--preprocess-qc-dir", type=str, default=None,
                   help="Stage 2 (preprocess_qc.R) output root, <dir>/<sample>/"
                        "soupx_corrected_matrix + doublet_calls.tsv. Only consulted in "
                        "--input-mode=matrix (fastq mode looks alongside --matrix-dir's "
                        "own outs/ subdirectory instead, unchanged).")
    p.add_argument("--out-dir",        required=True)
    p.add_argument("--samples",        required=True, nargs="+")
    p.add_argument("--min-genes",      type=int,   default=200)
    p.add_argument("--max-genes",      type=int,   default=6000)
    p.add_argument("--max-mito",       type=float, default=20.0)
    p.add_argument("--n-mads",         type=float, default=5.0,
                   help="MAD threshold for log1p_total_counts/log1p_n_genes_by_counts/"
                        "pct_counts_in_top_20_genes outlier detection")
    p.add_argument("--n-mads-mt",      type=float, default=3.0,
                   help="MAD threshold for pct_counts_mt outlier detection (tighter than --n-mads, "
                        "matching the sc-best-practices QC notebook)")
    p.add_argument("--remove-doublets", type=lambda s: s.lower() == "true", default=False,
                   help="Drop predicted_doublet==True cells during Stage 3b QC filtering "
                        "instead of just flagging them. Default false (flag-only, unchanged "
                        "behavior) -- inspect 03b_doublet_summary.pdf before enabling.")
    p.add_argument("--n-hvgs",         type=int,   default=3000)
    p.add_argument("--n-pcs",          type=int,   default=50)
    p.add_argument("--n-neighbors",    type=int,   default=15)
    p.add_argument("--resolution",     type=float, default=0.5)
    p.add_argument("--metadata-file",  type=str,   required=True,
                   help="Path to sample_metadata.tsv (tab-separated, must have 'sample' column)")
    p.add_argument("--viz-dir",        type=str,   default=None,
                   help="Directory for summary figure PDFs. Defaults to <out-dir>/../summary_figures")
    p.add_argument("--tables-dir",     type=str,   default=None,
                   help="Directory for summary data tables (TSV). Defaults to <out-dir>/../summary_tables")
    p.add_argument("--threads",        type=int,   default=1,
                   help="Threads for scanpy's own parallel dispatch (sc.settings.n_jobs) -- "
                        "e.g. neighbors/UMAP via pynndescent. BLAS threading (numpy/scipy) is "
                        "controlled separately via OMP_NUM_THREADS etc., set by run_pipeline.sh "
                        "before this process starts.")
    return p.parse_args()


# ── Verification helpers ──────────────────────────────────────────────────────
def verify_anndata(adata: ad.AnnData, stage: str) -> None:
    """Log shape and raise immediately if the AnnData is empty after a step."""
    log.info(f"[VERIFY SUCCESS] After {stage}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if adata.n_obs == 0:
        raise RuntimeError(
            f"[VERIFY FAIL] No cells remain after {stage}. "
            "Check QC thresholds (min/max genes, max mito) — they may be too strict."
        )
    if adata.n_vars == 0:
        raise RuntimeError(f"[VERIFY FAIL] No genes remain after {stage}.")


def verify_outputs(adata: ad.AnnData, out_dir: Path) -> None:
    """
    Confirm all expected output files exist on disk and that the AnnData
    object honours the cross-language contract required by downstream.R.
    Call this after all files have been written.
    """
    # ── File existence ────────────────────────────────────────────────────────
    expected = [
        "annotated.h5ad",
        "obs_metadata.csv",
        "umap_coords.csv",
        "cluster_markers.csv",
    ]
    missing_files = [f for f in expected if not (out_dir / f).exists()]
    if missing_files:
        raise FileNotFoundError(
            f"[VERIFY FAIL] Missing output files: {missing_files}"
        )

    # ── AnnData cross-language contract ──────────────────────────────────────
    errors = []
    if "counts" not in adata.layers:
        errors.append("adata.layers['counts'] missing — DESeq2 in downstream.R will fail")
    if adata.raw is None:
        errors.append("adata.raw not set — LIANA/NMF will use wrong matrix")
    if "cell_type" not in adata.obs.columns:
        errors.append("'cell_type' column missing from obs — downstream.R hardcodes this")
    if "X_umap" not in adata.obsm:
        errors.append("UMAP embedding missing from obsm — Monocle3 trajectory will fail")
    if errors:
        raise RuntimeError("[VERIFY FAIL] AnnData contract violations:\n  " + "\n  ".join(errors))

    # ── Summary ───────────────────────────────────────────────────────────────
    h5ad_mb = (out_dir / "annotated.h5ad").stat().st_size / 1e6
    log.info(
        f"[VERIFY SUCCESS] Output OK — "
        f"{adata.n_obs:,} cells, "
        f"{adata.obs['leiden'].nunique()} clusters, "
        f"{adata.obs['cell_type'].nunique()} cell types, "
        f"annotated.h5ad: {h5ad_mb:.1f} MB"
    )

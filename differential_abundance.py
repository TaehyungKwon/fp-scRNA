#!/usr/bin/env python3
"""
differential_abundance.py — Stage 6: Compositional differential abundance (scCODA)
Reads obs_metadata.csv written by scanpy_analysis.py; independent of Stage 7 (R)
so it can run right after Stage 3-5 without waiting on downstream.R.
Optional: skipped gracefully if sccoda isn't installed, or if the condition
column can't be auto-detected.
"""

import argparse
import re
import sys
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend -- safe on headless servers
import matplotlib.pyplot as plt

# ── Logging ───────────────────────────────────────────────────────────────────
# Matches scanpy_analysis.py's / downstream.R's [INFO]/[WARNING]/[ERROR] tags.
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s]\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Shared visual style ────────────────────────────────────────────────────────
# Duplicated from scanpy_analysis.py (same rationale as that file's own
# duplication note: fixed hex codes, not a colormap name, so R and Python
# figures read as one visual system despite being rendered by two libraries).
PALETTE_QUALITATIVE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
]
COLOR_ACCENT = "#D85A30"
COLOR_MUTED  = "#8C8C8C"

FIGURE_DPI = 300
MM_PER_IN = 25.4
PAGE_W_IN = 180 / MM_PER_IN
PAGE_H_IN = 240 / MM_PER_IN

plt.rcParams["font.family"]     = ["Liberation Sans", "sans-serif"]
plt.rcParams["axes.labelsize"]  = 6
plt.rcParams["xtick.labelsize"] = 6
plt.rcParams["ytick.labelsize"] = 6
plt.rcParams["legend.fontsize"] = 5
AXIS_FONTSIZE   = 6
LEGEND_FONTSIZE = 5

# ── Condition-column auto-detection ───────────────────────────────────────────
# Python port of downstream.R's PIPELINE_COLS/detect_condition_col(). The two
# copies must be kept in sync manually -- this pipeline shares no code across
# languages by design (see codeshare/CLAUDE.md), so any new obs column
# scanpy_analysis.py writes (QC metrics, outlier flags, doublet columns, ...)
# needs to be added to BOTH this list and downstream.R's PIPELINE_COLS, or it
# can silently become the condition/DA grouping variable instead of the real
# one (this already happened once in R with scDblFinder_class).
PIPELINE_COLS = {
    "sample", "cell_type", "leiden", "doublet_score", "predicted_doublet",
    "scDblFinder_score", "scDblFinder_class", "outlier", "mt_outlier",
    "n_genes_by_counts", "log1p_n_genes_by_counts",
    "total_counts", "log1p_total_counts", "pct_counts_in_top_20_genes",
    "total_counts_mt", "log1p_total_counts_mt", "pct_counts_mt",
    "total_counts_ribo", "log1p_total_counts_ribo", "pct_counts_ribo",
    "total_counts_hb", "log1p_total_counts_hb", "pct_counts_hb",
    "n_counts", "n_genes", "UMAP1", "UMAP2",
    "predicted_labels", "over_clustering", "majority_voting", "conf_score",
}


def detect_condition_col(obs: pd.DataFrame, pipeline_cols: set) -> str:
    candidates = []
    for col in obs.columns:
        if col in pipeline_cols or col.startswith("Unnamed"):
            continue
        n_unique = obs[col].nunique(dropna=True)
        if 2 <= n_unique <= 10:
            candidates.append(col)

    if not candidates:
        raise ValueError("Could not auto-detect a condition column. Check "
                          "sample_metadata.tsv has a condition-like column.")

    if len(candidates) > 1:
        log.info(f"Multiple condition candidates found: {candidates} — "
                 f"using '{candidates[0]}'. Set manually if wrong.")

    log.info(f"Condition column: '{candidates[0]}'")
    return candidates[0]


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Stage 6: scCODA compositional differential abundance")
    p.add_argument("--scanpy-dir",         required=True)
    p.add_argument("--tables-dir",         required=True)
    p.add_argument("--viz-dir",            required=True)
    p.add_argument("--control-condition",  default="",
                   help="Reference/baseline level of the auto-detected condition column. "
                        "Required if that column has >2 levels.")
    return p.parse_args()


def main():
    args = parse_args()
    scanpy_dir = Path(args.scanpy_dir)
    tables_dir = Path(args.tables_dir)
    viz_dir = Path(args.viz_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    viz_dir.mkdir(parents=True, exist_ok=True)
    control_condition = args.control_condition.strip()

    obs_path = scanpy_dir / "obs_metadata.csv"
    if not obs_path.exists():
        log.error(f"[VERIFY FAIL] {obs_path} not found — was Stage 3-5 completed?")
        sys.exit(1)
    obs = pd.read_csv(obs_path, index_col=0)

    if "cell_type" not in obs.columns:
        log.error("'cell_type' column not found in obs_metadata.csv. Was annotation completed?")
        sys.exit(1)

    cond_col = detect_condition_col(obs, PIPELINE_COLS)
    levels = sorted(obs[cond_col].dropna().unique().tolist())
    log.info(f"Conditions: {', '.join(levels)}")

    # Same hard-fail rule as downstream.R's Stage 7a (DESeq2): with >2 levels, an
    # arbitrary reference silently answers the wrong scientific question, so
    # require the user to say which level is control instead.
    if len(levels) > 2 and not control_condition:
        log.error(f"[VERIFY FAIL] Condition column '{cond_col}' has {len(levels)} levels "
                 f"({', '.join(levels)}) but CONTROL_CONDITION is not set. Set it in "
                 f"my_project.sh to one of the levels above.")
        sys.exit(1)
    if control_condition and control_condition not in levels:
        log.error(f"[VERIFY FAIL] CONTROL_CONDITION='{control_condition}' is not among the "
                 f"detected levels: {', '.join(levels)}")
        sys.exit(1)

    try:
        from sccoda.util import cell_composition_data as ccd
        from sccoda.util import comp_ana as ca
    except ImportError:
        log.warning("sccoda not available — skipping compositional differential abundance")
        sys.exit(0)

    # Sample x cell_type count matrix, one row per sample.
    counts = obs.groupby(["sample", "cell_type"]).size().unstack(fill_value=0)
    meta = obs.drop_duplicates("sample").set_index("sample")[[cond_col]]
    meta = meta.loc[counts.index]

    if control_condition:
        formula = f"C({cond_col}, Treatment('{control_condition}'))"
    else:
        formula = f"C({cond_col})"  # 2-level case, unset: patsy's default reference (first level)

    df = counts.join(meta)
    data = ccd.from_pandas(df, covariate_columns=[cond_col])
    log.info(f"Fitting scCODA model ({counts.shape[0]} samples x {counts.shape[1]} cell types, "
             f"formula: {formula})...")

    model = ca.CompositionalAnalysis(data, formula=formula, reference_cell_type="automatic")
    result = model.sample_hmc(verbose=False)

    _, effect_df = result.summary_prepare()
    effect_df = effect_df.reset_index()
    credible = result.credible_effects().rename("credible").reset_index()
    effect_df = effect_df.merge(credible, on=["Covariate", "Cell Type"])

    effect_df["condition_level"] = effect_df["Covariate"].str.extract(r"\[T\.(.+)\]$")
    # NOTE: "Final Parameter" (a CLR-space regression coefficient) and its
    # "hdi_3"/"hdi_97" are one consistent triple from the same posterior --
    # "log2_fold_change" is a *separately*-scaled derived quantity (log2 ratio
    # of expected sample composition) and its own credible interval isn't
    # reported by scCODA at all. Never plot hdi_3/hdi_97 as error bars around
    # log2_fold_change -- they don't bound it (final_parameter can and does
    # land outside its own log2_fold_change's neighborhood for shrunk/
    # non-credible effects). Keep both columns, but only final_parameter +
    # hdi_3/hdi_97 are a valid (point, interval) pair for plotting.
    out_df = effect_df.rename(columns={
        "Cell Type": "cell_type", "Final Parameter": "final_parameter",
        "log2-fold change": "log2_fold_change",
        "HDI 3%": "hdi_3", "HDI 97%": "hdi_97",
        "Inclusion probability": "inclusion_probability",
    })[["cell_type", "condition_level", "final_parameter", "hdi_3", "hdi_97",
        "log2_fold_change", "inclusion_probability", "credible"]]
    out_df.insert(2, "control_condition", control_condition or levels[0])
    out_df = out_df.sort_values(["condition_level", "cell_type"])

    out_path = tables_dir / "sccoda_effects.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    n_credible = int(out_df["credible"].sum())
    log.info(f"[VERIFY SUCCESS] scCODA — {len(out_df)} (cell_type x condition_level) effects "
             f"written to {out_path} ({n_credible} credible)")

    make_effects_plot(out_df, viz_dir / "06_sccoda_effects.pdf", control_condition or levels[0])

    log.info("Done.")


def make_effects_plot(df: pd.DataFrame, out_path: Path, control_label: str) -> None:
    levels = sorted(df["condition_level"].unique())
    n_levels = len(levels)
    fig, axes = plt.subplots(1, n_levels, figsize=(PAGE_W_IN, PAGE_H_IN * 0.45),
                              squeeze=False, sharey=True)
    axes = axes[0]

    cell_types = sorted(df["cell_type"].unique())
    y_pos = {ct: i for i, ct in enumerate(cell_types)}

    for ax, level in zip(axes, levels):
        sub = df[df["condition_level"] == level].set_index("cell_type").loc[cell_types]
        colors = [COLOR_ACCENT if c else COLOR_MUTED for c in sub["credible"]]
        # final_parameter + hdi_3/hdi_97 is the one internally-consistent
        # (point, interval) triple scCODA reports -- see the NOTE in main().
        lo = np.clip(sub["final_parameter"] - sub["hdi_3"], 0, None)
        hi = np.clip(sub["hdi_97"] - sub["final_parameter"], 0, None)
        ax.errorbar(sub["final_parameter"], [y_pos[c] for c in sub.index],
                    xerr=[lo, hi], fmt="none", ecolor="grey", elinewidth=0.6,
                    capsize=1.5, zorder=1)
        ax.scatter(sub["final_parameter"], [y_pos[c] for c in sub.index],
                   c=colors, s=10, zorder=2)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_title(f"{level}\nvs {control_label}", fontsize=7)
        ax.set_xlabel("effect (CLR coefficient)", fontsize=AXIS_FONTSIZE)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    axes[0].set_yticks(list(y_pos.values()))
    axes[0].set_yticklabels(list(y_pos.keys()), fontsize=5.5)
    fig.suptitle("scCODA compositional effects (credible = orange; whiskers = 94% HDI)",
                 fontsize=8, y=1.02)
    fig.tight_layout()
    out = out_path.with_suffix(".pdf")
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    log.info(f"[VIZ] scCODA effects plot: {out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log.exception(str(e))
        sys.exit(1)

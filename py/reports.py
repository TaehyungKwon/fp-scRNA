"""
py/reports.py — every summary figure scanpy_analysis.py writes: Stage 1's
Cell Ranger metrics (01_cellranger_summary.pdf), Stage 3b's QC/doublet
figures (03b_*.pdf), Stage 5's clustering/UMAP/dotplot figures
(05_clustering_summary.pdf, umap_*.pdf, dotplot__top_markers.pdf). Grouped
separately from the processing functions (py/stage3_qc.py etc.) since these
all read finished state and only write figures -- mirrors this file's own
original "Plots"/"Summary visualizations" section split.
"""

import csv
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from py.common import (
    save_fig, grid_figure, legend_below, category_colors, QC_METRICS,
    COLOR_BEFORE, COLOR_AFTER, COLOR_WARN, COLOR_MID, COLOR_MUTED,
    PAGE_W_IN, PAGE_H_IN, FIGURE_DPI, LEGEND_FONTSIZE,
)

log = logging.getLogger(__name__)


# ── Stage 5 plots (scanpy's own plotting functions) ───────────────────────────
def make_plots(adata: ad.AnnData, viz_dir: Path):
    """umap_overview.pdf (+ umap_leiden.pdf, umap_cell_type.pdf), dotplot__top_markers.pdf
    -- scanpy's own plotting functions, redirected to write straight into the shared
    summary_figures dir as PDF (scanpy picks the file format from the `save` string's
    extension; it also prepends its own plot-type prefix, hence the filenames below)."""
    log.info("Generating plots...")
    sc.settings.figdir = str(viz_dir)
    sc.settings.dpi = FIGURE_DPI

    # sc.pl.umap reads matplotlib's rcParams as the per-panel size (sc.settings.figsize
    # has no effect here). "sample" has a small, config-bounded category count, so it
    # shares this multi-panel grid with the two continuous QC metrics (colorbar) --
    # none of their legends/colorbars are wide enough to push a 2-column grid past
    # 180mm even with legend_loc enabled.
    old_figsize = plt.rcParams["figure.figsize"]
    plt.rcParams["figure.figsize"] = (PAGE_W_IN / 2 * 0.9, PAGE_H_IN / 3)
    sc.pl.umap(adata, color=["sample", "n_genes_by_counts", "pct_counts_mt"],
               ncols=2, show=False, legend_loc="right margin", save="_overview.pdf")
    plt.rcParams["figure.figsize"] = old_figsize

    # leiden/cell_type have unbounded, data-dependent category counts -- either
    # legend can be arbitrarily tall, and would blow a shared 2-column grid well
    # past the 180mm page envelope. Each gets its own single-panel file instead,
    # with the full page width available for its legend.
    for key in ("leiden", "cell_type"):
        # cell_type's legend text (CellTypist labels, e.g. "Tem/Effector helper
        # T cells") is far longer than leiden's bare cluster numbers -- shrink
        # the panel further and drop the legend font a point so the added
        # legend column doesn't push the page past 180mm.
        panel_w = PAGE_W_IN * (0.65 if key == "cell_type" else 0.95)
        legend_fs = LEGEND_FONTSIZE - 1 if key == "cell_type" else LEGEND_FONTSIZE
        plt.rcParams["figure.figsize"] = (panel_w, PAGE_H_IN * 0.55)
        sc.pl.umap(adata, color=key, show=False, legend_loc="right margin",
                   legend_fontsize=legend_fs, save=f"_{key}.pdf")
        plt.rcParams["figure.figsize"] = old_figsize

    # Dot plot: top 3 markers per cluster
    top_markers = (
        sc.get.rank_genes_groups_df(adata, group=None)
        .groupby("group")
        .head(3)["names"]
        .unique()
        .tolist()
    )
    # return_fig + .style(largest_dot=...) instead of save=, so dot size can be
    # scaled down (0.75 * 0.6 = 45% of scanpy's own DEFAULT_LARGEST_DOT=200 --
    # dots were still too large at the previous 75%, so shrunk a further 40%).
    dp = sc.pl.dotplot(adata, top_markers, groupby="cell_type", show=False, return_fig=True,
                       figsize=(PAGE_W_IN * 0.93, PAGE_H_IN * 0.6))
    dp.style(largest_dot=sc.pl.DotPlot.DEFAULT_LARGEST_DOT * 0.75 * 0.6)
    dp.savefig(str(viz_dir / "dotplot__top_markers.pdf"), dpi=FIGURE_DPI)
    log.info("  Plots saved.")


# ── Summary visualizations ────────────────────────────────────────────────────

def make_cellranger_summary(cellranger_dir: str, samples: list[str],
                             viz_dir: Path, tables_dir: Path) -> None:
    """
    01_cellranger_summary.pdf (2-column grid) + cellranger_metrics_summary.tsv
    One bar per sample for each key Cell Ranger metric parsed from
    metrics_summary.csv (estimated cells, median genes/cell, median UMI/cell,
    sequencing saturation). Handles column-name variation across CR versions
    by case-insensitive substring matching.
    """
    records: list[dict] = []
    for sample in samples:
        metrics_path = Path(cellranger_dir) / sample / "outs" / "metrics_summary.csv"
        if not metrics_path.exists():
            log.warning(f"  [VIZ] metrics_summary.csv not found for {sample} — skipped")
            continue
        with open(metrics_path, newline="") as fh:
            row = next(csv.DictReader(fh))
        # Normalise: strip whitespace/quotes, remove thousands commas
        norm = {
            k.strip().strip('"').lower(): v.strip().strip('"').replace(",", "")
            for k, v in row.items()
        }

        def _get(*keywords: str) -> float:
            for k, v in norm.items():
                if all(kw in k for kw in keywords):
                    try:
                        return float(v.rstrip("%"))
                    except ValueError:
                        pass
            return float("nan")

        records.append({
            "Sample":                sample,
            "Estimated Cells":       _get("estimated", "cell"),
            "Median Genes / Cell":   _get("median", "gene"),
            "Median UMI / Cell":     _get("median", "umi") or _get("median", "count"),
            "Seq. Saturation (%)":   _get("sequencing", "saturation"),
        })

    if not records:
        log.warning("[VIZ] No Cell Ranger metrics found; skipping 01_cellranger_summary")
        return

    df = pd.DataFrame(records).set_index("Sample")
    df.to_csv(tables_dir / "cellranger_metrics_summary.tsv", sep="\t")
    metrics = [c for c in df.columns if not df[c].isna().all()]
    cmap    = category_colors(df.index.tolist())
    colors  = [cmap[s] for s in df.index]

    fig, axes = grid_figure(len(metrics))

    for ax, metric in zip(axes, metrics):
        vals = df[metric].tolist()
        bars = ax.bar(df.index, vals, color=colors, edgecolor="white", linewidth=0.5)
        ax.set_title(metric, fontsize=10, fontweight="bold", pad=6)
        ax.set_ylabel(metric)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df.index, rotation=35, ha="right")
        ax.spines[["top", "right"]].set_visible(False)
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() * 1.01,
                        f"{val:,.0f}", ha="center", va="bottom", fontsize=7)

    fig.suptitle("Cell Ranger QC Summary — per sample", fontsize=13, fontweight="bold", y=1.02)
    out = save_fig(fig, viz_dir / "01_cellranger_summary")
    plt.close(fig)
    log.info(f"[VIZ] {out}")


def make_qc_summary(adata: ad.AnnData, cells_before_qc: dict[str, int],
                    viz_dir: Path) -> None:
    """
    03b_qc_summary.pdf (2-column grid)
    Panel 1: grouped bar — cells per sample before vs after QC.
    Panel 2: % cells retained per sample (dashed line at 80 % as reference).
    Per-metric before/after distributions live in 03b_qc_metrics_boxplot.pdf.
    """
    samples = list(cells_before_qc.keys())
    after   = adata.obs["sample"].value_counts()
    n_before = np.array([cells_before_qc.get(s, 0) for s in samples], dtype=float)
    n_after  = np.array([after.get(s, 0)            for s in samples], dtype=float)
    pct_kept = np.where(n_before > 0, 100 * n_after / n_before, 0.0)

    fig, axes = grid_figure(2)
    x = np.arange(len(samples))
    w = 0.38

    # Panel 1 — before / after
    axes[0].bar(x - w / 2, n_before, w, label="Before QC", color=COLOR_BEFORE, edgecolor="white")
    axes[0].bar(x + w / 2, n_after,  w, label="After QC",  color=COLOR_AFTER, edgecolor="white")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(samples, rotation=35, ha="right")
    axes[0].set_ylabel("Cell count")
    legend_below(axes[0], ncol=2)
    axes[0].set_title("Cells per sample:\nbefore vs after QC", fontweight="bold")
    axes[0].spines[["top", "right"]].set_visible(False)

    # Panel 2 — retention rate
    bar_colors = [COLOR_WARN if p < 50 else COLOR_MID if p < 80 else COLOR_AFTER for p in pct_kept]
    axes[1].bar(samples, pct_kept, color=bar_colors, edgecolor="white")
    axes[1].axhline(80, ls="--", color="gray", lw=1, label="80 % reference")
    axes[1].set_ylim(0, 110)
    axes[1].set_ylabel("% cells retained")
    legend_below(axes[1], ncol=1)
    axes[1].set_xticks(range(len(samples)))
    axes[1].set_xticklabels(samples, rotation=35, ha="right")
    axes[1].set_title("QC retention rate\nper sample", fontweight="bold")
    axes[1].spines[["top", "right"]].set_visible(False)
    for i, (s, p) in enumerate(zip(samples, pct_kept)):
        axes[1].text(i, p + 1.5, f"{p:.0f}%", ha="center", va="bottom", fontsize=7)

    fig.suptitle("QC Summary", fontsize=13, fontweight="bold")
    out = save_fig(fig, viz_dir / "03b_qc_summary")
    plt.close(fig)
    log.info(f"[VIZ] {out}")


def make_qc_metrics_boxplot(qc_before: pd.DataFrame, adata: ad.AnnData, viz_dir: Path) -> None:
    """
    03b_qc_metrics_boxplot.pdf (2-column grid)
    One panel per QC metric (QC_METRICS): before vs after filtering,
    side by side. Boxes show the bulk distribution (quartiles + whiskers);
    only points beyond the whiskers are drawn as individual dots -- the
    full per-cell distribution isn't plotted, just its outliers.
    """
    fig, axes = grid_figure(len(QC_METRICS), row_h_in=2.6)
    flier_style = dict(marker="o", markersize=1.5, markeredgewidth=0,
                       markerfacecolor=COLOR_WARN, alpha=0.5)

    for ax, metric in zip(axes, QC_METRICS):
        before = qc_before[metric].dropna().values
        after  = adata.obs[metric].dropna().values
        bp = ax.boxplot([before, after], tick_labels=["Before", "After"],
                        patch_artist=True, widths=0.55, flierprops=flier_style)
        for patch, color in zip(bp["boxes"], [COLOR_BEFORE, COLOR_AFTER]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for median in bp["medians"]:
            median.set_color("black")
        ax.set_title(metric, fontsize=9, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("QC metrics: before vs after filtering", fontsize=13, fontweight="bold")
    out = save_fig(fig, viz_dir / "03b_qc_metrics_boxplot")
    plt.close(fig)
    log.info(f"[VIZ] {out}")


def make_clustering_summary(adata: ad.AnnData, viz_dir: Path) -> None:
    """
    05_clustering_summary.pdf (2-column grid)
    Panel 1-3: UMAP coloured by leiden cluster, cell type, and sample.
    Panel 4: stacked bar — cell-type composition (%) per sample.
    """
    umap = adata.obsm["X_umap"]
    u1, u2 = umap[:, 0], umap[:, 1]

    fig, (ax1, ax2, ax3, ax4) = grid_figure(4, row_h_in=4.15)

    def _scatter(ax, labels: list[str], title: str) -> None:
        cmap = category_colors(labels)
        cats = sorted(cmap)
        for cat in cats:
            mask = np.array(labels) == cat
            ax.scatter(u1[mask], u2[mask], s=1.5, alpha=0.4,
                       color=cmap[cat], label=cat, rasterized=True)
        ax.set_title(title, fontweight="bold", fontsize=9)
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.spines[["top", "right"]].set_visible(False)
        ncol = max(1, len(cats) // 6)
        handles, lbls = ax.get_legend_handles_labels()
        ax.legend(handles, lbls, markerscale=5, fontsize=LEGEND_FONTSIZE,
                  loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=ncol,
                  framealpha=0.4, handletextpad=0.2, columnspacing=0.5)

    _scatter(ax1, adata.obs["leiden"].tolist(),    "Leiden clusters")
    _scatter(ax2, adata.obs["cell_type"].tolist(), "Cell types")
    _scatter(ax3, adata.obs["sample"].tolist(),    "Sample")

    # Panel 4 — stacked proportion bar
    prop = (
        adata.obs
        .groupby(["sample", "cell_type"], observed=True)
        .size()
        .unstack(fill_value=0)
        .apply(lambda r: r / r.sum() * 100, axis=1)
    )
    ct_cmap = category_colors(prop.columns.tolist())
    bottoms = np.zeros(len(prop))
    for ct in prop.columns:
        vals = prop[ct].values
        ax4.bar(prop.index, vals, bottom=bottoms,
                color=ct_cmap[ct], label=ct, width=0.6)
        bottoms += vals
    ax4.set_ylabel("% cells")
    ax4.set_title("Cell-type composition\nper sample", fontweight="bold", fontsize=9)
    ax4.set_xticks(range(len(prop)))
    ax4.set_xticklabels(prop.index, rotation=35, ha="right")
    ax4.spines[["top", "right"]].set_visible(False)
    legend_below(ax4, ncol=max(1, len(prop.columns) // 6))

    fig.suptitle("Clustering & Annotation Summary", fontsize=13, fontweight="bold")
    out = save_fig(fig, viz_dir / "05_clustering_summary")
    plt.close(fig)
    log.info(f"[VIZ] {out}")


def make_doublet_summary(adata: ad.AnnData, viz_dir: Path, tables_dir: Path) -> None:
    """
    03b_doublet_summary.pdf (2-column grid) + doublet_summary.tsv
    Panel 1: doublet_score distribution per sample (boxplot, outliers only --
             same style as 03b_qc_metrics_boxplot.pdf).
    Panel 2: doublet fraction (%) per sample.
    Panel 3: UMAP colored by predicted_doublet (singlet vs doublet).
    Doublets are flagged only by default (REMOVE_DOUBLETS=false, see
    stage3_qc.py's run_qc()) -- inspect this figure before deciding whether
    to enable removal.
    """
    if "doublet_score" not in adata.obs.columns:
        log.info("[VIZ] doublet_score not available (no scDblFinder/scrublet ran) "
                 "— skipping 03b_doublet_summary")
        return

    samples = sorted(adata.obs["sample"].unique())
    fig, axes = grid_figure(3, row_h_in=3.2)

    # Panel 1 -- doublet_score distribution per sample
    flier_style = dict(marker="o", markersize=1.5, markeredgewidth=0,
                       markerfacecolor=COLOR_WARN, alpha=0.5)
    scores_by_sample = [adata.obs.loc[adata.obs["sample"] == s, "doublet_score"].dropna().values
                       for s in samples]
    bp = axes[0].boxplot(scores_by_sample, tick_labels=samples, patch_artist=True,
                         widths=0.55, flierprops=flier_style)
    for patch in bp["boxes"]:
        patch.set_facecolor(COLOR_MID)
    axes[0].set_xticklabels(samples, rotation=35, ha="right")
    axes[0].set_ylabel("doublet_score")
    axes[0].set_title("Doublet score\nper sample", fontweight="bold")
    axes[0].spines[["top", "right"]].set_visible(False)

    # Panel 2 -- doublet fraction per sample
    frac = (adata.obs.groupby("sample")["predicted_doublet"].mean() * 100).reindex(samples)
    axes[1].bar(samples, frac.values, color=COLOR_WARN, edgecolor="white")
    axes[1].set_ylabel("% predicted doublets")
    axes[1].set_xticks(range(len(samples)))
    axes[1].set_xticklabels(samples, rotation=35, ha="right")
    axes[1].set_title("Doublet fraction\nper sample", fontweight="bold")
    axes[1].spines[["top", "right"]].set_visible(False)
    for i, v in enumerate(frac.values):
        axes[1].text(i, v + max(frac.values) * 0.03, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)

    # Panel 3 -- UMAP colored by predicted_doublet
    umap = adata.obsm["X_umap"]
    u1, u2 = umap[:, 0], umap[:, 1]
    is_doublet = adata.obs["predicted_doublet"].values.astype(bool)
    axes[2].scatter(u1[~is_doublet], u2[~is_doublet], s=1.5, alpha=0.3,
                    color=COLOR_MUTED, label="singlet", rasterized=True)
    axes[2].scatter(u1[is_doublet], u2[is_doublet], s=3, alpha=0.8,
                    color=COLOR_WARN, label="doublet", rasterized=True)
    axes[2].set_title("UMAP: predicted doublets", fontweight="bold", fontsize=9)
    axes[2].set_xlabel("UMAP 1"); axes[2].set_ylabel("UMAP 2")
    axes[2].spines[["top", "right"]].set_visible(False)
    legend_below(axes[2], ncol=2)

    fig.suptitle("Doublet Detection Summary", fontsize=13, fontweight="bold")
    out = save_fig(fig, viz_dir / "03b_doublet_summary")
    plt.close(fig)
    log.info(f"[VIZ] {out}")

    # Companion table, same shape/spirit as cellranger_metrics_summary.tsv
    summary_df = pd.DataFrame({
        "sample":       samples,
        "n_cells":      [int((adata.obs["sample"] == s).sum()) for s in samples],
        "n_doublets":   [int(adata.obs.loc[adata.obs["sample"] == s, "predicted_doublet"].sum())
                         for s in samples],
        "pct_doublets": frac.values,
    })
    summary_path = tables_dir / "doublet_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    log.info(f"  Doublet summary table: {summary_path}")

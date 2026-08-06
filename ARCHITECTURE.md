# ARCHITECTURE.md — scrna_pipeline

Contributor-facing reference for how this pipeline is put together: stages,
cross-language data contracts, configuration, and package choices. For
setup/quick-start, see README.md. For agent-specific operational notes
(known gotchas, environment quirks), see CLAUDE.md.

---

## Stages

Bash orchestrates seven stages across three languages. Numbers are
sequential top-level steps; letters mark multiple programs/sub-outputs
within one conceptual stage (e.g. Stage 2's two techniques, or Stage 7's six
independent optional analyses).

| Stage | Script | Purpose |
|-------|--------|---------|
| **1** | `run_pipeline.sh` → `cellranger count` (skipped when `INPUT_MODE=matrix`) | fastq → matrix (fastq mode only) |
| **2** | `preprocess_qc.R` (env `scrna_r`), per sample, both modes | Ambient RNA & doublet correction |
| 　2a | " | SoupX ambient correction (fastq-only — needs Cell Ranger's raw matrix) |
| 　2b | " | scDblFinder doublet detection (both modes) |
| **3** | `scanpy_analysis.py` → `py/stage3_qc.py` (env `scrna_py`) | Load & QC-filter |
| 　3a | " | Load & concatenate samples, merge sample metadata |
| 　3b | " | QC metrics & MAD-based filtering (+ scrublet fallback if 2b didn't run) |
| **4** | `scanpy_analysis.py` → `py/stage4_embed.py` | Normalize & embed |
| 　4a | " | Normalization & HVG selection |
| 　4b | " | PCA + Harmony batch correction (optional) + kNN + UMAP |
| **5** | `scanpy_analysis.py` → `py/stage5_cluster.py` | Cluster & annotate |
| 　5a | " | Leiden clustering + marker genes |
| 　5b | " | Cell type annotation (CellTypist, optional; else leiden cluster ID) |
| **6** | `differential_abundance.py` (env `scrna_py`) | Compositional differential abundance (scCODA) — independent, runs right after Stage 5 |
| **7** | `downstream.R` → `r/stage7*.R` (env `scrna_r`) | Downstream analysis, each letter an independent optional module |
| 　7a | `r/stage7a_deseq2.R` | Pseudobulk DEG (DESeq2) |
| 　7b | `r/stage7b_monocle3.R` | Trajectory (Monocle3, optional) |
| 　7c | `r/stage7c_liana.R` | Cell–cell communication (LIANA, optional) |
| 　7d | `r/stage7d_nmf.R` | Gene programs (NMF, optional) |
| 　7e | `r/stage7e_gsea.R` | GSEA (GO:BP + KEGG + Hallmark) |
| 　7f | `r/stage7f_milo.R` | Differential abundance via kNN neighborhoods (Milo) |

Stage 6 (scCODA) and Stage 7f (Milo) are the pipeline's two *alternative*
differential-abundance methods, testing the same underlying question
("did cell-type/neighborhood abundance change with condition?") via
different statistics — they don't run adjacently because of a data
dependency, not because they're unrelated: scCODA only needs
`obs_metadata.csv` (already written by Stage 5), so it runs immediately
after Stage 5, independent of Stage 7 (R); Milo needs PCA coordinates from
the annotated `.h5ad`, so it runs inside Stage 7 instead.

Every stage's Python/R implementation is split into a thin entry point
(`scanpy_analysis.py`, `downstream.R`) plus per-stage submodules
(`py/*.py`, `r/*.R`) that the entry point imports/sources — so a heavy
library (scanpy, DESeq2, a Bioconductor package) still loads once per
process, not once per stage, while each stage's code stays independently
reviewable. `preprocess_qc.R` and `differential_abundance.py` are each
already a single-purpose standalone script (Stage 2 and Stage 6
respectively), invoked as their own subprocess by `run_pipeline.sh`, so
they aren't split further.

Every stage's **human-facing** outputs (figures, tables) land in two flat,
shared directories instead of being scattered per-stage —
`${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}/` (PDF only) and
`${SUMMARY_TABLES_OUT}/${PROJECT_NAME}/` (TSV only), both siblings of
`PROJECT_OUT` rather than nested inside it. The 12 numbered summary figures
(`01_cellranger_summary.pdf` … `07f_milo_beeswarm__*.pdf`) are prefixed with
the stage that produced them; `umap_overview.pdf`/`umap_leiden.pdf`/
`umap_cell_type.pdf`/`dotplot__top_markers.pdf` (Stage 5) are auto-named by
scanpy's own plotting functions and left unprefixed. Table filenames
(`DEG_*.tsv`, `GSEA_*.tsv`, `milo_DA_*.tsv`, `sccoda_effects.tsv`, …)
already self-describe their analysis and aren't stage-prefixed.

---

## Cross-language data contracts

Stages communicate through files in `$CELLRANGER_OUT`/`$PREPROCESS_QC_OUT`/
`$PROJECT_OUT`, not shared code — this pipeline shares no code across
languages by design.

- **Two count layers.** `normalize()` (`py/stage4_embed.py`) stores raw
  integer counts in `adata.layers["counts"]` before `normalize_total`/
  `log1p`, then freezes log-norm in `adata.raw`. `downstream.R` reads
  `layers["counts"]` → assay `"counts"` (for DESeq2) and `adata.X` → assay
  `"X"` (for LIANA/NMF). Removing the `counts` layer hard-fails Stage 7a.
- **`cell_type` is a fixed contract column.** Always written by
  `cluster_and_annotate()` (`py/stage5_cluster.py`) — CellTypist labels if
  available, else leiden IDs as strings. `downstream.R` hardcodes
  `ctype_col <- "cell_type"`.
- **`sample` is the universal batch key.** Set at load time
  (`py/stage3_qc.py`'s `load_samples()`); used for HVG batch-awareness and
  as the Harmony key.
- **Condition column is auto-detected, independently, in both R and
  Python.** `r/common.R`'s `detect_condition_col()`/`PIPELINE_COLS` and
  `differential_abundance.py`'s port of the same logic both pick the first
  `obs` column not in their (separately maintained) `PIPELINE_COLS` list
  with 2–10 distinct values. **Any new pipeline-generated `obs` column with
  2–10 distinct values must be added to BOTH copies of `PIPELINE_COLS`, or
  it can silently become the DEG/DA grouping variable instead of the real
  condition.** Both lists must stay in sync with every `obs` column
  `scanpy_analysis.py` writes (QC metrics, `outlier`/`mt_outlier`, doublet
  columns).
- **`sample_metadata.tsv`'s `sample` column drives the sample list.**
  `run_pipeline.sh` reads it directly (no separate `SAMPLES=()` array to
  keep in sync); each value must match a subdirectory in `FASTQ_DIR`
  (fastq mode) or `MATRIX_DIR` (matrix mode). Duplicate `sample` values
  abort the run immediately.

---

## Multi-condition analysis (`CONTROL_CONDITION`)

`CONTROL_CONDITION` names the reference/baseline level of the
auto-detected condition column. It's consumed identically by Stage 6
(scCODA), Stage 7a (DESeq2), and Stage 7f (Milo):

- **Required if the condition column has more than 2 levels.** Without it,
  DESeq2's `results()` (and edgeR inside Milo's `testNhoods()`) would
  otherwise silently return only the last-vs-first-alphabetical
  coefficient — e.g. a 4-level `condition` column with values
  `cocktail`/`control`/`plastic5`/`plastic50` would silently compare
  `plastic50 vs cocktail`, never touching the real control.
  `differential_abundance.py` and `downstream.R` each independently
  enforce this and abort (exit 1) with an error naming the detected levels
  if it's unset — whichever stage runs first (Stage 6 runs before Stage 7)
  surfaces the error via `run_pipeline.sh`'s `set -o pipefail`.
- **Optional for exactly 2 levels** — only affects which level
  fold-changes are computed against; leaving it unset keeps prior behavior
  unchanged (an arbitrary but deterministic reference).
- **Every non-control level gets its own vs-control contrast**, named
  `<level>_vs_<CONTROL_CONDITION>` (sanitized like cell-type names,
  `[^A-Za-z0-9_]` → `_`): `DEG_<celltype>__<level>_vs_<control>.tsv`,
  `GSEA_<celltype>__<level>_vs_<control>.tsv`,
  `milo_DA_<level>_vs_<control>.tsv`. GSEA runs once per (cell type ×
  contrast) — accepted ~(n_levels−1)× runtime cost over the 2-level case.
- **Figure-splitting rule**: when a summary figure can't legibly hold cell
  type × contrast × collection on one static page, it splits into one file
  per non-control level: `07a_deg_summary__<level>_vs_<control>.pdf`,
  `07e_gsea_summary__<level>_vs_<control>.pdf`,
  `07f_milo_beeswarm__<level>_vs_<control>.pdf`.
  `07a_volcano_top_celltype.pdf` stays a single file (global best
  cell_type × contrast pair by DEG count). **scCODA is the one
  exception** — its whole point is one joint model across every condition
  level at once, so `sccoda_effects.tsv`/`06_sccoda_effects.pdf` stay
  singular, faceted internally by level instead.
- Milo needs replicated samples per condition just like DESeq2 does — with
  1 sample per condition, every contrast reports 0 residual degrees of
  freedom for dispersion estimation, and 0 `milo_DA_*.tsv` tables get
  written. scCODA (Bayesian) can still produce output in that case, just
  with limited statistical confidence.

---

## `PREPROCESSING_ONLY`

Set `true` to stop the pipeline after Stage 5 (Cell Ranger through QC,
normalization, embedding, clustering, and annotation), before Stage 6
(scCODA) and Stage 7 (R downstream: DEG/GSEA/trajectory/etc.) — useful for
reviewing clustering/annotation quality, or handing off the annotated
dataset, before committing to the statistical analysis stages. Default
`false` (run the full pipeline).

---

## config.sh reference

| Variable | Required | Description |
|----------|----------|-------------|
| `WORKDIR` | Yes (env var) | Root directory; export before running |
| `INPUT_MODE` | No | `"fastq"` (default) or `"matrix"` |
| `MATRIX_DIR` | Only if `INPUT_MODE=matrix` | Flat per-sample dir of pre-computed matrices |
| `FASTQ_DIR` | Only if `INPUT_MODE=fastq` | Dir with one subdir per sample containing `.fastq.gz` files |
| `GENOME_REF` | Only if `INPUT_MODE=fastq` | Cell Ranger genome reference directory |
| `OUT_DIR` | Yes | Base directory. `CELLRANGER_OUT`, `PREPROCESS_QC_OUT`, `PROJECT_OUT`, `SUMMARY_FIGURES_OUT`, `SUMMARY_TABLES_OUT` are computed automatically from it |
| `CELLRANGER_OUT` | No (auto-derived) | Shared across projects; per-sample Cell Ranger subdirs land here |
| `PREPROCESS_QC_OUT` | No (auto-derived) | Shared across projects; Stage 2 output for matrix-mode samples lands here (fastq mode uses `CELLRANGER_OUT`) |
| `PROJECT_OUT` | No (auto-derived) | Project dir = `${PROJECT_OUT}/${PROJECT_NAME}` for scanpy/downstream/logs |
| `SUMMARY_FIGURES_OUT` / `SUMMARY_TABLES_OUT` | No (auto-derived) | Figures/tables dir = `${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}` / `${SUMMARY_TABLES_OUT}/${PROJECT_NAME}` |
| `METADATA_FILE` | Yes | TSV with `sample` column + any condition columns; drives the sample list directly |
| `CONTROL_CONDITION` | Only if the condition column has >2 levels | Reference/baseline level for Stage 6/7a/7f — see "Multi-condition analysis" above |
| `CELLRANGER_DIR` | No | Dir containing `cellranger` binary; leave `""` if on PATH |
| `CONDA_BIN` | No | Dir containing `conda` binary; leave `""` if on PATH |
| `THREADS` / `MEM_GB` | No | Pipeline-wide TOTAL resource budget |
| `CELLRANGER_THREADS_PER_SAMPLE` / `CELLRANGER_MEM_GB_PER_SAMPLE` | No | Per-sample caps for Stage 1's parallel fan-out (default 16 / 64) |
| `MIN_GENES`, `MAX_GENES`, `MAX_MITO` | No | Coarse pre-filter floor/ceiling, applied before MAD-based outlier detection in Stage 3b |
| `REMOVE_DOUBLETS` | No | `false` (default) flags doublets only; `true` drops them during Stage 3b |
| `N_HVGS`, `N_PCS`, `N_NEIGHBORS`, `LEIDEN_RES` | No | Stage 4/5 clustering parameters |
| `PYTHON_ENV` / `R_ENV` | No | Conda env names (default: `scrna_py` / `scrna_r`) |
| `SETUP_ENVS` | No | Set `true` on first run to create conda envs |
| `SKIP_CELLRANGER` | No | fastq mode only: set `true` to skip Stage 1 when matrices already exist |
| `SKIP_AMBIENT_QC` | No | set `true` to skip Stage 2 entirely (both modes) |
| `PREPROCESSING_ONLY` | No | set `true` to stop after Stage 5 — see above |

---

## Package rationale

### Python (Scanpy stack, env `scrna_py`)

| Package      | Role                        | Optional? | Why minimal |
|--------------|-----------------------------|-----------|-------------|
| scanpy       | Core: Stage 3-5 (QC, norm, PCA, UMAP, clustering) | No | Replaces Seurat, Monocle2, many others |
| anndata      | Data container              | No        | Required by scanpy |
| numpy/pandas | Numerics + tables           | No        | Standard lib |
| scipy        | Sparse matrix, stats        | No        | scanpy dependency |
| matplotlib   | Static plots                | No        | No seaborn needed |
| leidenalg    | Stage 5a: Leiden clustering | No*       | scanpy calls it internally |
| scrublet     | Stage 3b: doublet detection fallback | Yes | Gracefully skipped if scDblFinder (Stage 2b) already ran |
| harmonypy    | Stage 4b: batch correction  | Yes       | Falls back to uncorrected PCA |
| celltypist   | Stage 5b: automated annotation | Yes    | Falls back to leiden cluster IDs |
| sccoda       | Stage 6: compositional differential abundance | Yes | Pulls in TensorFlow + tensorflow-probability (its Bayesian backend) — heavier than the rest of this env, but confirmed not to conflict with the `anndata<0.11` pin `zellkonverter` needs |

### R (Downstream, env `scrna_r`)

| Package                                  | Role                              | Optional? |
|-------------------------------------------|-----------------------------------|-----------|
| SoupX / scran                             | Stage 2a: ambient RNA correction  | Yes       |
| scDblFinder                               | Stage 2b: doublet detection       | Yes       |
| DropletUtils                              | Stage 2: 10x mtx read/write       | Yes       |
| DESeq2                                    | Stage 7a: pseudobulk DEG          | Yes       |
| zellkonverter                             | Read `.h5ad` in R                 | Yes       |
| clusterProfiler / msigdbr / org.Hs.eg.db  | Stage 7e: GSEA (GO:BP + KEGG + Hallmark) | Yes |
| miloR                                     | Stage 7f: kNN-neighborhood differential abundance | Yes |
| monocle3                                  | Stage 7b: trajectory / pseudotime | Yes       |
| liana                                     | Stage 7c: ligand–receptor communication | Yes |
| NMF                                       | Stage 7d: gene program discovery  | Yes       |
| ggplot2/dplyr                             | Plots + data wrangling            | No        |

`monocle3`, `liana`, and `NMF` are **intentionally not in `r_env.yaml`**
(conda installs are unreliable) — install them manually via
`install.packages`/`pak` after env creation. All other R Stage 7 packages
*are* in `r_env.yaml` but every Stage 7 module is still wrapped in a
`requireNamespace()` check — the script continues with the next module if
a package is missing, so you can install only what you need.

---

## FASTQ filename convention (`cellranger --sample`)

Cell Ranger's `--sample` argument must match the **filename prefix** (the
portion before `_S<N>_L<N>_`), not the folder name. The pipeline
auto-detects this from the first `.fastq.gz` file in each sample directory:

```
Donor1/5k_Human_Donor1_PBMC_3p_gem-x_GEX_S1_L001_R1_001.fastq.gz
         └──────── detected prefix ─────────────────┘
```

No manual mapping needed — just make sure each `sample` value in
`METADATA_FILE` matches its folder name.

---

## Figure typography

Every figure in both languages uses **Liberation Sans** (real Arial isn't
installed on Linux; Liberation Sans is the metrically-compatible open
substitute most Linux distros ship as its equivalent). Axis text is 6pt,
legend text 5pt, legends placed below the plot — `py/common.py`
(`AXIS_FONTSIZE`/`LEGEND_FONTSIZE` constants, `legend_below()` helper) and
`r/common.R`'s `FIGURE_THEME` (append `+ FIGURE_THEME` after each plot's
own `theme_classic(...)`/etc. call). R's plain `pdf()` device errors on
named font families — `save_plot()` (`r/common.R`) uses `device =
cairo_pdf`, which resolves fonts through the same fontconfig database
`fc-list` reads.

# CLAUDE.md — scrna_pipeline

Portable, Linux-native 10X Genomics scRNA-seq pipeline.  
Three-language orchestration: **Bash** (orchestrator) → **Python/Scanpy** (QC + clustering) → **R** (downstream analysis).

---

## Running the pipeline

See README.md's Quick Start for the standard invocation. To run a stage in
isolation, call `scanpy_analysis.py`/`downstream.R` directly with `conda run
-n <env>` — pass `--help` to either for the current flag set (they change
often enough that an example here would go stale; see `run_pipeline.sh`'s
`run_scanpy()`/`run_r_downstream()` for the exact flags a full run uses).

---

## Portability: how tool paths work

The pipeline does **not** require cellranger or conda to be on the system PATH.  
Set `CELLRANGER_DIR` and `CONDA_BIN` in your `config.sh` and the orchestrator injects them into PATH at runtime:

```bash
# In config.sh:
CELLRANGER_DIR="/opt/cellranger-10.0.0"   # dir containing the cellranger binary
CONDA_BIN="/home/user/miniforge3/bin"      # dir containing the conda binary
```

Leave either empty (`""`) to rely on the system PATH instead.  
Nothing is written to `~/.bashrc` or any system file.

**`python3` and `Rscript` are intentionally NOT checked in the outer PATH** — they live inside the `scrna_py` / `scrna_r` conda envs and are called via `conda run -n <env>`. The `check_deps` function probes them that way.

---

## Architecture: the cross-language data contract

Stages communicate through files in `$CELLRANGER_OUT` / `$PROJECT_OUT`, not shared code:

| Stage | Script | Output |
|-------|--------|--------|
| 1 | `run_pipeline.sh` → `cellranger count` (skipped when `INPUT_MODE=matrix`) | `cellranger/<sample>/outs/filtered_feature_bc_matrix/`, `.../raw_feature_bc_matrix/` |
| 2 (2a/2b) | `preprocess_qc.R` (env `scrna_r`), per sample, both modes (2a: SoupX correction, fastq-only; 2b: scDblFinder, both modes) | fastq: `.../outs/soupx_corrected_matrix/`, `.../outs/doublet_calls.tsv`. matrix: `${PREPROCESS_QC_OUT}/<sample>/doublet_calls.tsv` |
| 3–5 | `scanpy_analysis.py` (env `scrna_py`) → `py/stage3_qc.py`, `py/stage4_embed.py`, `py/stage5_cluster.py`, `py/reports.py` | `scanpy/annotated.h5ad` (+ CSV side-cars used only by Stage 6/7) |
| 6 | `differential_abundance.py` (env `scrna_py`) — scCODA, independent of Stage 7 | `sccoda_effects.tsv`, `06_sccoda_effects.pdf` |
| 7 (7a–7f) | `downstream.R` (env `scrna_r`) → `r/stage7a_deseq2.R` … `r/stage7f_milo.R`, `r/stage7_summary.R` | in-memory only — no per-module subdirectories |

`scanpy_analysis.py` and `downstream.R` are thin entry points (CLI parsing,
orchestration) over these per-stage submodules — `py/common.py`/`r/common.R`
hold the infrastructure shared across all of them (logging, the visual
style/page-envelope constants, verification helpers). `preprocess_qc.R` and
`differential_abundance.py` are each already single-purpose standalone
scripts (Stage 2 and Stage 6 respectively) and aren't split further.

Every stage's **human-facing** outputs (figures, tables) land in two flat,
shared directories instead of being scattered per-module —
`${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}/` (PDF only) and
`${SUMMARY_TABLES_OUT}/${PROJECT_NAME}/` (TSV only), both siblings of
`PROJECT_OUT` under `OUT_DIR` rather than nested inside it — so every
project's figures/tables are browsable side-by-side without descending into
each project's own directory; see README.md for the current file list. The
12 numbered summary figures (`01_cellranger_summary.pdf` …
`07f_milo_beeswarm__*.pdf`) are prefixed with the stage that produced them.
Every figure is sized to fit a 180mm x 240mm page (`PAGE_W_IN`/`PAGE_H_IN`
in `py/common.py`, `PAGE_W_MM`/`PAGE_H_MM` in `r/common.R`) and multi-panel
figures use a 2-column grid (`grid_figure()` in Python; no direct R
equivalent needed since `save_plot()`'s fixed mm size already prevents
overflow). Both languages share one hardcoded hex palette
(`PALETTE_QUALITATIVE` / `COLOR_*`) so a given sample or cell type gets the
same color everywhere.

### Starting from existing feature matrices (`INPUT_MODE=matrix`)

Set `INPUT_MODE="matrix"` and `MATRIX_DIR` in config.sh to skip Cell Ranger
entirely — useful when matrices were already produced elsewhere (a different
pipeline, a GEO/public dataset download, etc.). Layout is **flat**, unlike
Cell Ranger's own nested `outs/filtered_feature_bc_matrix/`:

```
MATRIX_DIR/<sample>/matrix.mtx.gz
MATRIX_DIR/<sample>/features.tsv.gz
MATRIX_DIR/<sample>/barcodes.tsv.gz
```

`<sample>` must match a `sample` value in `METADATA_FILE`, same as fastq mode.
This bypasses `FASTQ_DIR`, `GENOME_REF`, `CELLRANGER_DIR`, and `SKIP_CELLRANGER` — none of
those are read when `INPUT_MODE=matrix`. `run_pipeline.sh` verifies all three
matrix files exist per sample before invoking Stage 3 (`verify_matrix_inputs`),
mirroring `verify_cellranger_sample` for fastq mode. `py/stage3_qc.py`'s
`load_samples()` takes `--input-mode` to pick the right path layout; the
Cell Ranger metrics-summary visualization (`01_cellranger_summary.pdf`) is
skipped in matrix mode since `metrics_summary.csv` doesn't exist for
externally-sourced matrices.

### QC: ambient RNA correction, MAD-based filtering, doublets

Follows the sc-best-practices QC notebook's approach rather than fixed
thresholds:

- **Ambient RNA correction (2a: SoupX)** and **doublet detection (2b:
  scDblFinder)** run in `preprocess_qc.R`, once per sample, *before* Stage 3
  combines samples — doublet detection is meaningless on pooled multi-batch
  data, so this can't be folded into Stage 3's already-combined AnnData.
  **SoupX (2a) is fastq-mode only** — it needs Cell Ranger's raw (unfiltered)
  matrix for its ambient profile (`SoupChannel(raw, filtered)`), which matrix
  mode never has. **scDblFinder (2b) has no such requirement and runs in both
  modes** — `preprocess_qc.R` itself needed no changes for this, since it
  already skips SoupX gracefully (logging why) whenever `--raw-dir` is
  omitted and runs scDblFinder unconditionally afterward; only
  `run_pipeline.sh`'s Stage 2 gating and `py/stage3_qc.py`'s `load_samples()`
  needed updating. Skip the whole stage with `SKIP_AMBIENT_QC=true`.
  - fastq mode: output lands in `${CELLRANGER_OUT}/<sample>/outs/`, same as
    before. `load_samples()` prefers `outs/soupx_corrected_matrix/` over
    `filtered_feature_bc_matrix/` when present, and merges
    `outs/doublet_calls.tsv` into `.obs` by barcode.
  - matrix mode: output lands in `${PREPROCESS_QC_OUT}/<sample>/` instead —
    a new OUT_DIR-derived dir, shared across projects and keyed by sample
    like `CELLRANGER_OUT`, kept separate from `MATRIX_DIR` since that's
    user-provided external data the pipeline shouldn't write into. Only
    `doublet_calls.tsv` is meaningful here (the `soupx_corrected_matrix/`
    written alongside it is an uncorrected passthrough copy, since SoupX
    never ran) — `py/stage3_qc.py` still reads the count matrix straight
    from `MATRIX_DIR`, and only looks in `PREPROCESS_QC_OUT` for
    `doublet_calls.tsv`.
- **Cell filtering is MAD-based (`is_outlier()` in `py/stage3_qc.py`), not
  fixed min/max thresholds.** A cell is dropped if `log1p_total_counts`,
  `log1p_n_genes_by_counts`, or `pct_counts_in_top_20_genes` is more than
  `--n-mads` (default 5) median absolute deviations from the median, or if
  `pct_counts_mt` is more than `--n-mads-mt` (default 3) MADs from the median
  *or* exceeds `--max-mito`. `--min-genes`/`--max-genes` still run first as a
  coarse floor/ceiling, so obviously-broken barcodes don't skew the
  median/MAD statistics themselves — they're a pre-filter now, not the
  primary filtering mechanism.
- **Doublets are flagged, not removed, by default.** `predicted_doublet`/`doublet_score`
  stay in `.obs` (from scDblFinder if `preprocess_qc.R` ran, else a scrublet
  fallback) for inspection during clustering/visualization. Inspect
  `03b_doublet_summary.pdf` (score distribution + fraction per sample + UMAP
  colored by `predicted_doublet`) and `doublet_summary.tsv`, then set
  `REMOVE_DOUBLETS=true` in config.sh/my_project.sh if you want them dropped
  during Stage 3b QC filtering (same place MAD-based outliers are dropped, in
  `py/stage3_qc.py`'s `run_qc()`) — default stays `false` (flag-only,
  unchanged prior behavior).

### Key contract details (do not break these)

- **Two count layers.** `py/stage4_embed.py`'s `normalize()` stores raw integer counts in `adata.layers["counts"]` before `normalize_total`/`log1p`, then freezes log-norm in `adata.raw`. `downstream.R` reads `layers["counts"]` → assay `"counts"` (for DESeq2) and `adata.X` → assay `"X"` (for LIANA/NMF). Removing the `counts` layer hard-fails Stage 7a.
- **`cell_type` is a fixed contract column.** Always written by `py/stage5_cluster.py`'s `cluster_and_annotate()` (CellTypist labels if available, else leiden IDs as strings). `downstream.R` hardcodes `ctype_col <- "cell_type"`.
- **`sample` is the universal batch key.** Set at load time; used for HVG batch-awareness and as the Harmony key.
- **Condition column is auto-detected, independently, in both R and Python.** `r/common.R`'s `detect_condition_col()`/`PIPELINE_COLS` and `differential_abundance.py`'s port of the same logic both pick the first `obs` column not in their (separately maintained) `PIPELINE_COLS` list with 2–10 distinct values. **Any new pipeline-generated `obs` column with 2–10 distinct values must be added to BOTH copies of `PIPELINE_COLS`, or it can silently become the DEG/DA grouping variable instead of the real condition** — this actually happened once in R when `scDblFinder_class` (singlet/doublet) was added without updating the list, and Stage 7 silently ran DESeq2/GSEA on "singlet vs doublet" instead of the real experimental condition. Both lists must stay in sync with every `obs` column `scanpy_analysis.py` writes (QC metrics, `outlier`/`mt_outlier`, doublet columns).
- **`SAMPLES` is derived from `METADATA_FILE`'s `sample` column, not a config.sh array.** `run_pipeline.sh` reads the `sample` column directly (`awk` against the TSV header, in file order) right after config loads — there's nothing left in config.sh to fall out of sync with `METADATA_FILE`. Each value must still match a subdirectory in `FASTQ_DIR`/`MATRIX_DIR`. Duplicate `sample` values abort the run immediately (`die`).

---

## Multi-condition analysis (`CONTROL_CONDITION`)

`CONTROL_CONDITION` names the reference/baseline level of the auto-detected
condition column. It's consumed identically by Stage 6 (scCODA), Stage 7a
(DESeq2), and Stage 7f (Milo):

- **Required if the condition column has more than 2 levels.** Without it,
  DESeq2's `results()` (and edgeR inside Milo's `testNhoods()`) would otherwise
  silently return only the last-vs-first-alphabetical coefficient — e.g. a
  4-level `condition` column with values `cocktail`/`control`/`plastic5`/
  `plastic50` would silently compare `plastic50 vs cocktail`, never touching
  the real control. `differential_abundance.py` and `downstream.R` each
  independently enforce this and abort (exit 1) with an error naming the
  detected levels if it's unset — whichever stage runs first (Stage 6 runs
  before Stage 7) surfaces the error via `run_pipeline.sh`'s `set -o pipefail`.
- **Optional for exactly 2 levels** — only affects which level fold-changes
  are computed against; leaving it unset keeps prior behavior unchanged (an
  arbitrary but deterministic reference).
- **Every non-control level gets its own vs-control contrast**, named
  `<level>_vs_<CONTROL_CONDITION>` (sanitized like cell-type names,
  `[^A-Za-z0-9_]` → `_`): `DEG_<celltype>__<level>_vs_<control>.tsv`,
  `GSEA_<celltype>__<level>_vs_<control>.tsv`,
  `milo_DA_<level>_vs_<control>.tsv`. GSEA runs once per (cell type × contrast)
  — accepted ~(n_levels−1)× runtime cost over the 2-level case.
- **Figure-splitting rule**: when a summary figure can't legibly hold
  cell-type × contrast × collection on one static page, it splits into one
  file per non-control level — same rationale as `umap_overview.pdf`'s split
  into `umap_leiden.pdf`/`umap_cell_type.pdf` to avoid legend overflow.
  `07a_deg_summary__<level>_vs_<control>.pdf`,
  `07e_gsea_summary__<level>_vs_<control>.pdf`,
  `07f_milo_beeswarm__<level>_vs_<control>.pdf` each split this way;
  `07a_volcano_top_celltype.pdf` stays a single file (global best
  cell_type × contrast pair by DEG count). **scCODA is the one exception** —
  its whole point is one joint model across every condition level at once
  (`C(condition, Treatment('<control>'))` in a single `sample_hmc()` fit), so
  `sccoda_effects.tsv`/`06_sccoda_effects.pdf` stay singular, faceted
  internally by level instead.
- **Milo needs replicated samples per condition just like DESeq2 does** —
  `testNhoods()` calls edgeR's `estimateDisp()` internally, which needs
  residual degrees of freedom the same way DESeq2 does. With 1 sample per
  condition, every contrast logs "No residual df: setting dispersion to NA"
  and 0 `milo_DA_*.tsv` tables get written — not a bug, the same
  experimental-design limit as DESeq2, just surfacing a second time. scCODA
  (Bayesian) can still produce output in this case, just with limited
  statistical confidence.
- **Stage 6 (`differential_abundance.py`, scCODA)** reads only
  `obs_metadata.csv` (already written by Stage 3-5) — independent of Stage 7,
  runs right after it, doesn't wait on `downstream.R`. `sccoda`'s reference
  cell type is auto-selected (its own dispersion-based heuristic) — no
  separate config var for it.
- **Stage 7f (Milo, in `r/stage7f_milo.R`)** needs PCA coordinates, so it only
  runs on the `annotated.h5ad`/zellkonverter path — the `counts_matrix.csv`
  fallback near the top of `downstream.R` never builds an `sce` object, so
  Milo is skipped with a log line in that case. Builds its kNN graph on
  `reducedDim(sce, "X_pca_harmony")` if present else `"X_pca"` (mirrors
  `py/stage4_embed.py`'s own `use_rep` fallback) — both land in the h5ad via
  `adata.obsm` automatically, no extra Python-side export needed.

---

## config.sh reference

| Variable | Required | Description |
|----------|----------|-------------|
| `WORKDIR` | Yes (env var) | Root directory; export before running |
| `INPUT_MODE` | No | `"fastq"` (default) or `"matrix"` — see below |
| `MATRIX_DIR` | Only if `INPUT_MODE=matrix` | Flat per-sample dir of pre-computed matrices |
| `FASTQ_DIR` | Only if `INPUT_MODE=fastq` | Dir with one subdir per sample containing `.fastq.gz` files |
| `GENOME_REF` | Only if `INPUT_MODE=fastq` | Cell Ranger genome reference directory |
| `OUT_DIR` | Yes | Base directory. `CELLRANGER_OUT`, `PREPROCESS_QC_OUT`, `PROJECT_OUT`, `SUMMARY_FIGURES_OUT`, `SUMMARY_TABLES_OUT` are computed automatically inside `run_pipeline.sh` from `OUT_DIR` (`${OUT_DIR}/out_cellranger`, `out_preprocess_qc`, `out_project`, `summary_figures`, `summary_tables`) — config.sh/my_project.sh don't set them at all; only add one explicitly if you need it somewhere other than under `OUT_DIR` |
| `CELLRANGER_OUT` | No (auto-derived from `OUT_DIR`) | Shared across projects; per-sample Cell Ranger subdirs land here |
| `PREPROCESS_QC_OUT` | No (auto-derived from `OUT_DIR`) | Shared across projects; Stage 2 (scDblFinder) output for matrix-mode samples lands here (fastq mode still uses `CELLRANGER_OUT`) |
| `PROJECT_OUT` | No (auto-derived from `OUT_DIR`) | Project dir = `${PROJECT_OUT}/${PROJECT_NAME}` for scanpy/downstream/logs (machine-facing only) |
| `SUMMARY_FIGURES_OUT` / `SUMMARY_TABLES_OUT` | No (auto-derived from `OUT_DIR`) | Figures/tables dir = `${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}` / `${SUMMARY_TABLES_OUT}/${PROJECT_NAME}` — shared across projects, siblings of `PROJECT_OUT` rather than nested inside it |
| `METADATA_FILE` | Yes | TSV with `sample` column + any condition columns. `run_pipeline.sh` derives the sample list from this column directly — not a config.sh setting, and there is no `SAMPLES` variable to set |
| `CONTROL_CONDITION` | Only if the condition column has >2 levels | Reference/baseline level for DEG (DESeq2), GSEA, Milo, and scCODA — see "Multi-condition analysis" above |
| `CELLRANGER_DIR` | No | Dir containing `cellranger` binary; leave `""` if on PATH |
| `CONDA_BIN` | No | Dir containing `conda` binary; leave `""` if on PATH |
| `THREADS` / `MEM_GB` | No | Pipeline-wide TOTAL resource budget — Stage 1 divides it across concurrent Cell Ranger jobs (see below); Stage 3-5/7 use `THREADS` directly (BLAS env vars, `sc.settings.n_jobs`, R `mclapply`) |
| `CELLRANGER_THREADS_PER_SAMPLE` / `CELLRANGER_MEM_GB_PER_SAMPLE` | No | Per-sample caps for Stage 1's parallel fan-out (default 16 / 64) — `floor(THREADS/this)` and `floor(MEM_GB/this)` bound how many samples run concurrently |
| `MIN_GENES`, `MAX_GENES`, `MAX_MITO` | No | Coarse pre-filter floor/ceiling, applied before MAD-based outlier detection (see "QC" above); `MAX_MITO` also doubles as `pct_counts_mt`'s hard cap |
| `REMOVE_DOUBLETS` | No | `false` (default) flags doublets only; `true` drops `predicted_doublet` cells during Stage 3b QC — inspect `03b_doublet_summary.pdf` first |
| `N_HVGS`, `N_PCS`, `N_NEIGHBORS`, `LEIDEN_RES` | No | Stage 4/5 clustering parameters |
| `PYTHON_ENV` / `R_ENV` | No | Conda env names (default: `scrna_py` / `scrna_r`) |
| `SETUP_ENVS` | No | Set `true` on first run to create conda envs |
| `SKIP_CELLRANGER` | No | fastq mode only: set `true` to skip Stage 1 when matrices already exist in `CELLRANGER_OUT`'s nested layout (different from `INPUT_MODE=matrix`'s flat layout) |
| `SKIP_AMBIENT_QC` | No | set `true` to skip Stage 2 entirely (both modes) — SoupX (fastq only) + scDblFinder (both modes) |
| `PREPROCESSING_ONLY` | No | set `true` to stop after Stage 5 (Cell Ranger through clustering/annotation), skipping Stage 6 (scCODA) and Stage 7 (R downstream) |

---

## FASTQ filename convention (cellranger --sample)

Cell Ranger's `--sample` argument must match the **filename prefix** (the portion before `_S<N>_L<N>_`), not the folder name. The pipeline auto-detects this from the first `.fastq.gz` file in each sample directory:

```
Donor1/5k_Human_Donor1_PBMC_3p_gem-x_GEX_S1_L001_R1_001.fastq.gz
         └──────── detected prefix ─────────────────┘
```

No manual mapping needed — just make sure each `sample` value in `METADATA_FILE` matches its folder name.

---

## Graceful degradation — do not "fix" it

Every optional dependency is wrapped so a missing package downgrades rather than crashes:

- **Python**: `scrublet` (doublets), `harmonypy` (batch correction → raw PCA fallback), `celltypist` (annotation → leiden IDs fallback), `sccoda` (Stage 6 differential abundance) — all `try/except ImportError`. `sccoda` is heavier than the rest of this env (pulls in full TensorFlow + tensorflow-probability as its Bayesian backend, plus `tf-keras` to fix an import error the base package has without it) — confirmed not to conflict with the `anndata<0.11` pin zellkonverter needs (see `envs/python_env.yaml`'s comment). `pertpy`'s modern scCODA reimplementation was considered instead but transitively upgrades `anndata` past 0.11, breaking Stage 7's `readH5AD` — avoid it.
- **R**: Every Stage 7 module opens with `requireNamespace(...)` and returns early if the package is absent. `monocle3`, `liana`, and `NMF` are **intentionally not in `r_env.yaml`** (conda installs are unreliable) — install them manually via `install.packages`/`pak` after env creation. `clusterProfiler`/`msigdbr`/`org.Hs.eg.db` (Stage 7e GSEA), `SoupX`/`scDblFinder`/`DropletUtils`/`scran`/`scater` (Stage 2), and `miloR` (Stage 7f differential abundance) *are* in `r_env.yaml` but still guarded the same way — `preprocess_qc.R` falls back to uncorrected counts if `SoupX`/`scran` are missing, and skips doublet detection if `scDblFinder` is missing.

When adding features, guard new optional deps the same way and log a warning rather than aborting.

---

## Known gotchas

- `WORKDIR` must be **exported** before `source config.sh` runs — it is not defined inside config.sh. Always use `export WORKDIR=... && bash run_pipeline.sh` or set it in your shell profile.
- `run_pipeline.sh` runs under `set -euo pipefail`. Any unset variable or failed command will abort the run.
- `SETUP_ENVS` and `SKIP_CELLRANGER` set on the command line take precedence over values in config.sh (they are captured before the file is sourced).
- Cell Ranger is re-entrant: a sample whose `filtered_feature_bc_matrix` directory already exists is skipped automatically.
- The full pipeline run log is one file, shared across projects under `${OUT_DIR}/logs/pipeline_${PROJECT_NAME}_<timestamp>.log` (not nested under `PROJECT_OUT`) — every `log()`/`warn()`/`die()` call and every stage's own stdout/stderr (`2>&1 | tee -a "$LOG_FILE"`) funnels into it. Per-sample Cell Ranger + Stage 2 logs are separate, fine-grained files that live alongside their own shared output instead: `${CELLRANGER_OUT}/logs/` (fastq mode) or `${PREPROCESS_QC_OUT}/logs/` (matrix mode). Figures/tables land under `${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}/` and `${SUMMARY_TABLES_OUT}/${PROJECT_NAME}/`, not under `PROJECT_OUT`.
- `conda run -n <env>` does not activate the environment interactively — it is a subprocess call. Scripts that rely on conda shell functions (`conda activate`, `deactivate`) will not work inside the pipeline.
- `scanpy_analysis.py`/`downstream.R` are thin entry points over `py/*.py`/`r/*.R` submodules (one file per stage) — when editing a stage's logic, edit its submodule (e.g. `py/stage5_cluster.py`, `r/stage7a_deseq2.R`), not the entry point. `py/common.py`/`r/common.R` hold the infrastructure shared across submodules (logging, visual style, verification helpers); adding a new stage submodule means adding one `import`/`source()` line to the entry point, not restructuring it.
- `downstream.R` self-locates its own directory via `commandArgs()`'s `--file=` argument (R has no Python-style "directory of the running script" built-in) so it can `source()` `r/*.R` by absolute path regardless of the caller's working directory — don't replace this with a relative `source("r/common.R")`, which breaks the moment `run_pipeline.sh` is invoked from anywhere other than this directory.

---

## Logging convention

Every stdout/stderr line across all three languages starts with exactly one
of `[INFO]\t`, `[WARNING]\t`, `[ERROR]\t` (tab-separated from the message),
in addition to a leading timestamp:

- **Bash**: `log()` → `[INFO]`, `warn()` → `[WARNING]`, `die()` → `[ERROR]`
  (and exits 1). `plog()` is a file-only `[INFO]` variant for parallel
  workers, to avoid interleaved terminal output when several samples log
  concurrently — errors still go through `warn()`/`die()` so they surface
  immediately.
- **Python**: stdlib `logging`, format `"[%(asctime)s] [%(levelname)s]\t%(message)s"`.
  `%(levelname)s` is already exactly `INFO`/`WARNING`/`ERROR`. Uncaught
  exceptions are caught at the `if __name__ == "__main__":` boundary and
  logged via `log.exception(...)` (tag + traceback) before `sys.exit(1)` —
  otherwise they'd print a bare untagged Python traceback.
- **R**: `log_msg()`/`log_info()` → `[INFO]`, `log_warn()` → `[WARNING]`,
  `log_error()` → `[ERROR]`. `options(error = ...)` installs a global handler
  so uncaught `stop()` calls are logged via `log_error()` before `quit(status=1)`
  — otherwise they'd print R's default untagged `Error: ...`.

`[VERIFY SUCCESS]` / `[VERIFY WARN]` / `[VERIFY FAIL]` markers from the
`/verify`-style verification helpers are message-body text, not a substitute
for the level tag — e.g. a verification failure is `[ERROR]\t[VERIFY FAIL] ...`.

---

## Figure typography

Every figure in both languages uses **Liberation Sans** (real Arial isn't
installed on Linux; Liberation Sans is the metrically-compatible open
substitute most Linux distros ship as its equivalent — visually
indistinguishable, and requesting "Arial" directly just logs a per-text-
element "not found" warning before falling back anyway). Axis text is 6pt,
legend text 5pt, legends placed below the plot — `plt.rcParams` in
`py/common.py` (`AXIS_FONTSIZE`/`LEGEND_FONTSIZE` constants,
`legend_below()` helper) and `FIGURE_THEME` in `r/common.R` (append
`+ FIGURE_THEME` after each plot's own `theme_classic(...)`/etc. call).
**R's plain `pdf()` device errors on named font families** ("invalid font
type") — `save_plot()` (`r/common.R`) uses `device = cairo_pdf`, which
resolves fonts through the same fontconfig database `fc-list` reads.

---

## Development environment note

This pipeline was originally developed/tested on WSL (Ubuntu userland on Windows) and is targeted to run in production on an Ubuntu remote server. Script edits made from other environments (e.g. macOS) are done as plain text/syntax-checked only (`bash -n`, `python3 -m py_compile`) — actual execution (conda envs, `cellranger`, real data) should be verified on the Ubuntu target, since this repo is a portability layer between machines, not a dev environment itself.

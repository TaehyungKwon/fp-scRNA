Finn-Perkins lab reference pipeline for scRNA-seq data analysis
=======
# 10X Genomics scRNA-seq Pipeline

Portable, Linux-native scRNA-seq pipeline using Bash, Python (Scanpy), and R (DESeq2 + friends).  
Designed for minimal external dependencies — every optional package degrades gracefully if absent.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the full stage-by-stage reference
(data contracts, config.sh reference, package rationale). See
[CLAUDE.md](CLAUDE.md) for agent-specific operational notes.

---

## Directory layout

```
scrna_pipeline/
├── run_pipeline.sh          # Main orchestrator
├── config.sh                # All user-editable parameters
├── preprocess_qc.R          # Stage 2: ambient RNA correction (SoupX) + doublet detection (scDblFinder)
├── scanpy_analysis.py       # Stage 3-5 entry point: QC, normalization, clustering
├── py/                      # Stage 3-5 submodules (imported by scanpy_analysis.py)
│   ├── common.py            #   shared infra: logging, visual style, verify, CLI
│   ├── stage3_qc.py         #   3a: load/merge metadata, 3b: QC filtering
│   ├── stage4_embed.py      #   4a: normalize/HVG, 4b: PCA/Harmony/UMAP
│   ├── stage5_cluster.py    #   5a: Leiden, 5b: CellTypist annotation
│   └── reports.py           #   every summary figure
├── differential_abundance.py  # Stage 6: compositional differential abundance (scCODA)
├── downstream.R             # Stage 7 entry point: DEG, trajectory, LR, NMF, GSEA, Milo
├── r/                       # Stage 7 submodules (sourced by downstream.R)
│   ├── common.R             #   shared infra: logging, visual style, condition detection
│   ├── stage7a_deseq2.R
│   ├── stage7b_monocle3.R
│   ├── stage7c_liana.R
│   ├── stage7d_nmf.R
│   ├── stage7e_gsea.R
│   ├── stage7f_milo.R
│   └── stage7_summary.R
└── envs/
    ├── python_env.yaml      # Conda env for Python
    └── r_env.yaml           # Conda env for R
```

---

## Quick start

```bash
# 1. Copy and edit config.sh for your project
cp config.sh my_project.sh
nano my_project.sh

# 2. Create conda environments (first run only)
SETUP_ENVS=true bash run_pipeline.sh my_project.sh

# 3-1. Run the full pipeline
bash run_pipeline.sh my_project.sh

# 3-2. Skip Cell Ranger if matrices already exist (same nested outs/ layout)
SKIP_CELLRANGER=true bash run_pipeline.sh my_project.sh

# 3-3. Start from pre-computed matrices from elsewhere (flat layout)
    # Set in my_project.sh: INPUT_MODE="matrix"  MATRIX_DIR="/path/to/matrices"
    #   MATRIX_DIR/<sample>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
bash run_pipeline.sh my_project.sh
```

---

## Output structure

Every human-facing output (figures, result tables) lands under two flat,
shared base directories regardless of which stage produced it —
`${SUMMARY_FIGURES_OUT}/<PROJECT_NAME>/` (PDF only, each sized to fit a
180mm x 240mm page) and `${SUMMARY_TABLES_OUT}/<PROJECT_NAME>/` (TSV only).
Both derive from `OUT_DIR` by default and are kept out of `PROJECT_OUT` so
every project's figures/tables are browsable side-by-side without
descending into each project's own directory. `${OUT_DIR}/logs/` holds the
full pipeline run log; `${CELLRANGER_OUT}`/`${PREPROCESS_QC_OUT}` hold
per-sample Stage 1/2 output; `${PROJECT_OUT}/<PROJECT_NAME>/scanpy/` holds
the machine-facing `.h5ad`/CSV cross-language contract.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the exact figure/table filenames
per stage and the `<level>_vs_<control>` contrast-naming convention.

---

## Package rationale

See [ARCHITECTURE.md](ARCHITECTURE.md#package-rationale) for the full
per-package table (role, optional?, why). In short: every optional
dependency (`scrublet`, `harmonypy`, `celltypist`, `sccoda` in Python;
`SoupX`, `scDblFinder`, `DESeq2`, `monocle3`, `liana`, `NMF`, `miloR`,
GSEA's `clusterProfiler`/`msigdbr`/`org.Hs.eg.db` in R) degrades gracefully
if missing — the corresponding stage is skipped with a log line rather than
crashing the run.

---

## Portability notes

- Tested on Ubuntu 20.04+ and CentOS 7+ (x86_64).
- All Python/R dependencies managed through conda — no system-level `apt`/`yum` required beyond conda itself.
- `cellranger` must be installed separately (download from 10x Genomics website and add to `$PATH`).
- `bash >= 4.0` required (`associative arrays`, `set -euo pipefail`).
- Script is re-entrant: already-completed Cell Ranger samples are skipped automatically.

---

## Tuning parameters

| Parameter       | Default | Notes |
|----------------|---------|-------|
| `MIN_GENES`    | 200     | Coarse pre-filter floor, applied before MAD-based outlier detection |
| `MAX_GENES`    | 6000    | Coarse pre-filter ceiling, same |
| `MAX_MITO`     | 20%     | Hard cap on `pct_counts_mt`, combined with the 3-MAD check below; use 10–15% for neurons (naturally high) |
| `N_MADS`       | 5       | MAD threshold for total_counts/n_genes/pct_counts_in_top_20_genes outliers -- 5 is permissive, matches sc-best-practices |
| `N_MADS_MT`    | 3       | MAD threshold for `pct_counts_mt` specifically -- tighter than `N_MADS` |
| `N_HVGS`       | 3000    | 2000 for fast runs, 5000 for discovery |
| `N_PCS`        | 50      | Check elbow plot; 30–50 typical |
| `LEIDEN_RES`   | 0.5     | 0.3 = coarse, 1.2 = fine-grained |
| `N_NEIGHBORS`  | 15      | Lower for rarer cell types |
>>>>>>> e77e0e1 (Initial public release)

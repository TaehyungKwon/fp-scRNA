#!/usr/bin/env bash
# =============================================================================
# 10X Genomics scRNA-seq Pipeline -- Main Orchestrator
# Compatible: Linux x86_64, bash >= 4.0
# Dependencies: cellranger, conda/mamba (manages Python + R envs)
# Usage: bash run_pipeline.sh [config.sh]
# =============================================================================
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/config.sh}"

# -- Preserve caller-supplied flags before config.sh can overwrite them -------
_PRE_SETUP_ENVS="${SETUP_ENVS:-}"
_PRE_SKIP_CELLRANGER="${SKIP_CELLRANGER:-}"
_PRE_SKIP_AMBIENT_QC="${SKIP_AMBIENT_QC:-}"
_PRE_PREPROCESSING_ONLY="${PREPROCESSING_ONLY:-}"

# -- Load config -------------------------------------------------------------
if [[ ! -f "$CONFIG" ]]; then
    echo -e "[ERROR]\tConfig file not found: $CONFIG" >&2
    exit 1
fi
source "$CONFIG"

[[ -n "$_PRE_SETUP_ENVS"        ]] && SETUP_ENVS="$_PRE_SETUP_ENVS"
[[ -n "$_PRE_SKIP_CELLRANGER"   ]] && SKIP_CELLRANGER="$_PRE_SKIP_CELLRANGER"
[[ -n "$_PRE_SKIP_AMBIENT_QC"   ]] && SKIP_AMBIENT_QC="$_PRE_SKIP_AMBIENT_QC"
[[ -n "$_PRE_PREPROCESSING_ONLY" ]] && PREPROCESSING_ONLY="$_PRE_PREPROCESSING_ONLY"

# INPUT_MODE defaults to "fastq" so config.sh files written before this option
# existed keep working unchanged.
INPUT_MODE="${INPUT_MODE:-fastq}"
[[ "$INPUT_MODE" == "fastq" || "$INPUT_MODE" == "matrix" ]] \
    || { echo -e "[ERROR]\tINPUT_MODE must be 'fastq' or 'matrix' (got: $INPUT_MODE)" >&2; exit 1; }

# Per-sample resource caps for the parallel Cell Ranger fan-out (Stage 1).
# THREADS/MEM_GB are the pipeline-wide budget; these control how many samples
# run concurrently within it. Defaulted here so config.sh files written
# before this option existed keep working unchanged.
CELLRANGER_THREADS_PER_SAMPLE="${CELLRANGER_THREADS_PER_SAMPLE:-16}"
CELLRANGER_MEM_GB_PER_SAMPLE="${CELLRANGER_MEM_GB_PER_SAMPLE:-64}"

# MAD-based QC outlier thresholds (see scanpy_analysis.py's is_outlier()).
# Defaulted here so config.sh files written before this option existed keep
# working unchanged.
N_MADS="${N_MADS:-5}"
N_MADS_MT="${N_MADS_MT:-3}"

# Doublets are flagged (predicted_doublet/doublet_score in .obs), not removed,
# by default -- inspect 03b_doublet_summary.pdf before setting this true.
REMOVE_DOUBLETS="${REMOVE_DOUBLETS:-false}"

# All five OUT_DIR-derived dirs default from OUT_DIR -- config.sh/my_project.sh
# only need to set OUT_DIR itself; override any of these individually here
# only if you need one somewhere other than under OUT_DIR.
CELLRANGER_OUT="${CELLRANGER_OUT:-${OUT_DIR}/out_cellranger}"
PROJECT_OUT="${PROJECT_OUT:-${OUT_DIR}/out_project}"
SUMMARY_FIGURES_OUT="${SUMMARY_FIGURES_OUT:-${OUT_DIR}/summary_figures}"
SUMMARY_TABLES_OUT="${SUMMARY_TABLES_OUT:-${OUT_DIR}/summary_tables}"
# Stage 2 output, shared across projects and keyed by sample, same rationale
# as CELLRANGER_OUT -- kept separate from MATRIX_DIR (matrix mode) since that's
# user-provided external data the pipeline shouldn't write into.
PREPROCESS_QC_OUT="${PREPROCESS_QC_OUT:-${OUT_DIR}/out_preprocess_qc}"

# Stage 2 (2b: scDblFinder doublet detection, both modes; 2a: SoupX ambient
# RNA correction, fastq mode only -- it needs Cell Ranger's raw matrix, which
# matrix mode never has). Set true to skip Stage 2 entirely.
SKIP_AMBIENT_QC="${SKIP_AMBIENT_QC:-false}"

# Stop the pipeline after Stage 5 (Cell Ranger through QC, normalization,
# embedding, clustering, and annotation), before Stage 6 (scCODA) and
# Stage 7 (R downstream: DEG/GSEA/trajectory/etc.) -- useful for reviewing
# clustering/annotation quality (or handing off the annotated dataset)
# before committing to the statistical analysis stages. Default false (run
# the full pipeline, unchanged behavior).
PREPROCESSING_ONLY="${PREPROCESSING_ONLY:-false}"

# Reference/baseline level of the auto-detected condition column for DEG
# (DESeq2), GSEA, Milo, and scCODA. Required if that column has >2 levels --
# Stage 6 (scCODA) and Stage 7a/7f (DESeq2/Milo) enforce this themselves and
# abort with a clear error rather than silently picking an arbitrary
# reference. Optional (may be left "") for exactly 2 levels.
CONTROL_CONDITION="${CONTROL_CONDITION:-}"

# -- PATH bootstrap ----------------------------------------------------------
# If CELLRANGER_DIR / CONDA_BIN are set in config.sh and not already on PATH,
# inject them. Users set these in config.sh rather than relying on system-wide
# installations, keeping the pipeline portable across environments.
[[ -n "${CELLRANGER_DIR:-}" && -d "$CELLRANGER_DIR" ]] && export PATH="${CELLRANGER_DIR}:${PATH}"
[[ -n "${CONDA_BIN:-}"      && -d "$CONDA_BIN"      ]] && export PATH="${CONDA_BIN}:${PATH}"

# -- Project directory -------------------------------------------------------
# PROJECT_DIR is where machine-facing per-project outputs (scanpy, downstream,
# logs) land. Human-facing outputs (figures, tables) go under the shared
# SUMMARY_FIGURES_OUT/SUMMARY_TABLES_OUT dirs instead, each with its own
# PROJECT_NAME subdir -- kept out of PROJECT_DIR so they're browsable across
# projects without descending into each project's own directory. Cell Ranger
# outputs go to CELLRANGER_OUT, which is shared across projects.
[[ -z "${PROJECT_NAME:-}" ]] && { echo -e "[ERROR]\tPROJECT_NAME is not set in config.sh" >&2; exit 1; }
PROJECT_DIR="${PROJECT_OUT}/${PROJECT_NAME}"
SUMMARY_FIGURES_DIR="${SUMMARY_FIGURES_OUT}/${PROJECT_NAME}"
SUMMARY_TABLES_DIR="${SUMMARY_TABLES_OUT}/${PROJECT_NAME}"
mkdir -p "$PROJECT_DIR" "$SUMMARY_FIGURES_DIR" "$SUMMARY_TABLES_DIR"
[[ "$INPUT_MODE" == "fastq" ]] && mkdir -p "$CELLRANGER_OUT"
[[ "$SKIP_AMBIENT_QC" != "true" ]] && mkdir -p "$PREPROCESS_QC_OUT"

# -- Logging -----------------------------------------------------------------
# The full pipeline run log is one file, shared across projects under
# OUT_DIR/logs (not nested in PROJECT_DIR) -- every log()/warn()/die() call
# and every stage's `2>&1 | tee -a "$LOG_FILE"` output funnels into it, so
# PROJECT_NAME is baked into the filename to tell concurrent/past runs of
# different projects apart in that shared directory. Stage 1/2 per-sample
# logs still live alongside their shared, per-sample output instead
# (CELLRANGER_OUT in fastq mode; PREPROCESS_QC_OUT in matrix mode, since
# there's no Cell Ranger output to colocate with there) -- CR_LOG_DIR is
# always defined, since Stage 2's _run_preprocess_qc_one() writes into it
# in both modes.
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$LOG_DIR"
if [[ "$INPUT_MODE" == "fastq" ]]; then
    CR_LOG_DIR="${CELLRANGER_OUT}/logs"
else
    CR_LOG_DIR="${PREPROCESS_QC_OUT}/logs"
fi
mkdir -p "$CR_LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_${PROJECT_NAME}_$(date +%Y%m%d_%H%M%S).log"
export LOG_FILE  # log()/warn()/die()/plog() run inside parallel xargs -P subshells too

log()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]\t$*"    | tee -a "$LOG_FILE"; }
warn() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [WARNING]\t$*" | tee -a "$LOG_FILE"; }
die()  { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR]\t$*"   | tee -a "$LOG_FILE"; exit 1; }
# File-only progress line for parallel workers -- avoids interleaved terminal
# output when many samples log concurrently. Warnings/errors still use
# warn()/die() so they surface immediately.
plog() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]\t$*" >> "$LOG_FILE"; }
export -f log warn die plog

# -- Samples (derived from METADATA_FILE) ------------------------------------
# No manual SAMPLES=() array in config.sh -- every sample ID comes straight
# from METADATA_FILE's "sample" column instead, in the order rows appear
# there. Each value must still match a subdirectory name in FASTQ_DIR (fastq
# mode) or MATRIX_DIR (matrix mode) -- same requirement as before, just no
# longer duplicated by hand in two places that could drift out of sync.
[[ -n "${METADATA_FILE:-}" ]] || die "METADATA_FILE is not set in config.sh"
[[ -f "$METADATA_FILE" ]] || die "Metadata file not found: $METADATA_FILE"

sample_col=$(head -1 "$METADATA_FILE" | awk -F'\t' '{for (i=1;i<=NF;i++) if ($i=="sample") print i}')
[[ -n "$sample_col" ]] || die "METADATA_FILE has no 'sample' column: $METADATA_FILE"

mapfile -t SAMPLES < <(awk -F'\t' -v col="$sample_col" 'NR>1 && $col!="" {print $col}' "$METADATA_FILE")
(( ${#SAMPLES[@]} > 0 )) || die "No sample rows found in METADATA_FILE: $METADATA_FILE"

dup_samples=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d)
[[ -z "$dup_samples" ]] || die "Duplicate sample name(s) in METADATA_FILE: $(echo "$dup_samples" | tr '\n' ' ')"

# -- Resource-aware parallel fan-out (see Stage 1 below for why) -------------
compute_parallel_jobs() {
    local threads_per_job="$1" mem_gb_per_job="$2"
    local by_threads=$(( THREADS / threads_per_job ))
    local by_mem=$(( MEM_GB / mem_gb_per_job ))
    local jobs=$(( by_threads < by_mem ? by_threads : by_mem ))
    (( jobs < 1 )) && jobs=1
    (( jobs > ${#SAMPLES[@]} )) && jobs=${#SAMPLES[@]}
    echo "$jobs"
}

# -- Dependency checks -------------------------------------------------------
check_conda() {
    # Minimal check before env setup -- only conda itself must be on PATH.
    command -v conda &>/dev/null \
        || die "conda not found. Set CONDA_BIN in config.sh or install miniconda/mambaforge."
    log "conda found: $(conda --version)"
}

check_deps() {
    # Full check run after conda envs exist.
    # python3 and Rscript live inside their conda envs, not the outer PATH,
    # so we probe them via 'conda run' rather than 'command -v'.
    log "Checking dependencies..."
    local missing=()
    if [[ "$INPUT_MODE" == "fastq" ]]; then
        command -v cellranger &>/dev/null \
            || missing+=("cellranger (set CELLRANGER_DIR in config.sh)")
    fi
    conda run -n "$PYTHON_ENV" python3 --version &>/dev/null \
        || missing+=("python3 in conda env '$PYTHON_ENV' (run with SETUP_ENVS=true first)")
    conda run -n "$R_ENV" Rscript --version &>/dev/null \
        || missing+=("Rscript in conda env '$R_ENV' (run with SETUP_ENVS=true first)")
    [[ ${#missing[@]} -gt 0 ]] && die "Missing: ${missing[*]}"
    log "All dependencies found."
}

# -- Verification helpers ----------------------------------------------------

# Confirm Cell Ranger finished cleanly for one sample and surface key metrics.
verify_cellranger_sample() {
    local sample="$1"
    local outs_dir="${CELLRANGER_OUT}/${sample}/outs"
    local matrix_dir="${outs_dir}/filtered_feature_bc_matrix"
    local metrics="${outs_dir}/metrics_summary.csv"

    # Required matrix files
    local missing=()
    for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        [[ -f "${matrix_dir}/${f}" ]] || missing+=("$f")
    done
    [[ ${#missing[@]} -gt 0 ]] \
        && die "[VERIFY FAIL] $sample — Cell Ranger matrix incomplete, missing: ${missing[*]}"

    # Column names vary by CR version (matched case-insensitively on substrings).
    # Values are quoted with thousands-separator commas (e.g. "5,710"), which a
    # naive IFS=',' split would misalign -- shell out to Python's csv module
    # instead (matches make_cellranger_summary()'s parsing in scanpy_analysis.py).
    if [[ -f "$metrics" ]]; then
        local parsed
        parsed=$(conda run -n "$PYTHON_ENV" python3 -c "
import csv
with open('$metrics', newline='') as fh:
    row = next(csv.DictReader(fh))
cells = median_genes = sat = '?'
for k, v in row.items():
    kl = k.strip().lower()
    if 'estimated' in kl and 'cell' in kl: cells = v
    if 'median' in kl and 'gene' in kl: median_genes = v
    if 'sequencing' in kl and 'saturation' in kl: sat = v
print(f'{cells}\t{median_genes}\t{sat}')
" 2>/dev/null)
        IFS=$'\t' read -r cells median_genes sat <<< "$parsed"
        log "  [VERIFY SUCCESS] $sample — estimated cells: ${cells:-?}, median genes/cell: ${median_genes:-?}, sequencing saturation: ${sat:-?}"
    fi

    log "  [VERIFY SUCCESS] $sample: Cell Ranger output OK"
}

# Confirm a pre-computed feature matrix exists for one sample (INPUT_MODE=matrix).
# Layout is flat: MATRIX_DIR/<sample>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
# -- unlike Cell Ranger's nested <sample>/outs/filtered_feature_bc_matrix/.
verify_matrix_sample() {
    local sample="$1"
    local matrix_dir="${MATRIX_DIR}/${sample}"

    local missing=()
    for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        [[ -f "${matrix_dir}/${f}" ]] || missing+=("$f")
    done
    [[ ${#missing[@]} -gt 0 ]] \
        && die "[VERIFY FAIL] $sample — matrix incomplete in ${matrix_dir}, missing: ${missing[*]}"

    log "  [VERIFY SUCCESS] $sample: matrix files OK (${matrix_dir})"
}

verify_matrix_inputs() {
    log "=== Input mode: matrix (skipping Cell Ranger) ==="
    [[ -n "${MATRIX_DIR:-}" ]] || die "MATRIX_DIR must be set in config.sh when INPUT_MODE=matrix"
    [[ -d "$MATRIX_DIR" ]] || die "MATRIX_DIR not found: $MATRIX_DIR"
    for sample in "${SAMPLES[@]}"; do
        verify_matrix_sample "$sample"
    done
    log "All matrix inputs verified."
}

# Confirm the two conda envs can actually import their critical packages.
verify_envs() {
    log "Verifying conda environments..."
    # conda run consumes its own stdin before passing control to the subprocess,
    # so a heredoc fed to "python3 -" is silently swallowed — python3 sees EOF,
    # exits 0, and prints nothing.  Write the check to a temp file instead.
    local py_check
    py_check=$(mktemp /tmp/scrna_pycheck_XXXXXX.py)
    printf 'import importlib.metadata as _m\nprint(f"""  [VERIFY SUCCESS] scrna_py: scanpy={_m.version(\"scanpy\")}, anndata={_m.version(\"anndata\")}""")\n' \
        > "$py_check"
    conda run -n "$PYTHON_ENV" python3 "$py_check" 2>&1 | tee -a "$LOG_FILE" \
        || { rm -f "$py_check"; die "[VERIFY FAIL] scrna_py env broken — could not import scanpy/anndata"; }
    rm -f "$py_check"

    conda run -n "$R_ENV" Rscript -e \
        'cat("  [VERIFY SUCCESS] scrna_r:", R.version$version.string, "\n")' \
        2>&1 | tee -a "$LOG_FILE" \
        || die "[VERIFY FAIL] scrna_r env broken — Rscript failed"

    log "Environment verification passed."
}

# Confirm scanpy stage outputs exist on disk.
verify_scanpy_outputs() {
    local scanpy_dir="${PROJECT_DIR}/scanpy"
    local missing=()
    for f in annotated.h5ad obs_metadata.csv umap_coords.csv cluster_markers.csv; do
        [[ -f "${scanpy_dir}/${f}" ]] || missing+=("$f")
    done
    [[ ${#missing[@]} -gt 0 ]] \
        && die "[VERIFY FAIL] Scanpy stage incomplete — missing outputs: ${missing[*]}"
    local h5ad_bytes
    h5ad_bytes=$(wc -c < "${scanpy_dir}/annotated.h5ad")
    log "[VERIFY SUCCESS] Scanpy outputs OK (annotated.h5ad: ${h5ad_bytes} bytes)"
}

# Confirm R downstream produced at least the DEG summary tables (DESeq2 is the
# only Stage 7 module always installed; monocle3/liana/NMF/GSEA degrade
# gracefully and may legitimately produce nothing).
verify_r_outputs() {
    local tables_dir="$SUMMARY_TABLES_DIR"
    local n_deg
    n_deg=$(find "$tables_dir" -maxdepth 1 -name 'DEG_*.tsv' 2>/dev/null | wc -l)
    [[ "$n_deg" -eq 0 ]] \
        && die "[VERIFY FAIL] Downstream R stage produced no DEG_*.tsv tables in ${tables_dir}"
    log "[VERIFY SUCCESS] Downstream R outputs OK (${n_deg} DEG table(s) in ${tables_dir})"
}

# -- Stage 1: Cell Ranger (parallel across samples) --------------------------
# cellranger count's own multithreading stops scaling well before THREADS is
# exhausted, so fan out across samples via xargs -P instead, each capped at
# CELLRANGER_THREADS_PER_SAMPLE/CELLRANGER_MEM_GB_PER_SAMPLE.
_run_cellranger_one() {
    local sample="$1"
    local sample_out="${CELLRANGER_OUT}/${sample}"
    if [[ -d "${sample_out}/outs/filtered_feature_bc_matrix" ]]; then
        plog "  Skipping $sample (already complete)"
        return 0
    fi
    plog "  Processing sample: $sample"

    # Auto-detect the FASTQ sample prefix from the first matching file.
    # cellranger --sample must equal the portion of the filename before
    # _S<number>_L<number>_ -- it is NOT the same as the folder name.
    local fastq_dir="${FASTQ_DIR}/${sample}"
    [[ -d "$fastq_dir" ]] || die "FASTQ directory not found: $fastq_dir"
    local fastq_prefix
    fastq_prefix=$(ls "${fastq_dir}/"*.fastq.gz 2>/dev/null \
        | head -1 | xargs -I{} basename {} \
        | sed 's/_S[0-9]\+_L[0-9]\+_.*//')
    [[ -z "$fastq_prefix" ]] && die "No .fastq.gz files found in ${fastq_dir}/"
    plog "  FASTQ sample prefix ($sample): $fastq_prefix"

    cellranger count \
        --id="$sample" \
        --create-bam true \
        --transcriptome="$GENOME_REF" \
        --fastqs="$fastq_dir" \
        --sample="$fastq_prefix" \
        --localcores="$CELLRANGER_THREADS_PER_SAMPLE" \
        --localmem="$CELLRANGER_MEM_GB_PER_SAMPLE" \
        --output-dir="$sample_out" \
        >> "${CR_LOG_DIR}/${sample}_cellranger.log" 2>&1
    plog "  Done: $sample"
}
export -f _run_cellranger_one

run_cellranger() {
    log "=== Stage 1: Cell Ranger count ==="
    mkdir -p "$CELLRANGER_OUT"

    local jobs
    jobs=$(compute_parallel_jobs "$CELLRANGER_THREADS_PER_SAMPLE" "$CELLRANGER_MEM_GB_PER_SAMPLE")
    log "  Parallel Cell Ranger jobs: $jobs (of ${#SAMPLES[@]} samples; ${CELLRANGER_THREADS_PER_SAMPLE} threads / ${CELLRANGER_MEM_GB_PER_SAMPLE} GB each, within a ${THREADS}-thread / ${MEM_GB}GB total budget)"
    export CELLRANGER_OUT FASTQ_DIR GENOME_REF CR_LOG_DIR CELLRANGER_THREADS_PER_SAMPLE CELLRANGER_MEM_GB_PER_SAMPLE
    printf '%s\n' "${SAMPLES[@]}" \
        | xargs -I{} -P "$jobs" bash -c 'set -euo pipefail; _run_cellranger_one "$@"' _ {}

    for sample in "${SAMPLES[@]}"; do
        verify_cellranger_sample "$sample"
    done
    log "Stage 1 complete."
}

# -- Stage 2: Ambient RNA correction (2a: SoupX) + doublet detection (2b: scDblFinder) -
# Also per-sample, parallel, like Stage 1 -- and must run before Stage 3
# combines samples, since doublet detection is meaningless on pooled
# multi-batch data. scDblFinder (2b) needs no raw/unfiltered matrix, so it
# runs in BOTH modes; SoupX (2a) needs Cell Ranger's raw matrix for its
# ambient profile, so --raw-dir is only passed in fastq mode --
# preprocess_qc.R itself already skips SoupX gracefully (logging why)
# whenever --raw-dir is omitted, so no changes were needed there. Skip this
# whole stage with SKIP_AMBIENT_QC=true.
_run_preprocess_qc_one() {
    local sample="$1"
    local filtered_dir out_dir raw_dir_args=()
    if [[ "$INPUT_MODE" == "matrix" ]]; then
        filtered_dir="${MATRIX_DIR}/${sample}"
        out_dir="${PREPROCESS_QC_OUT}/${sample}"
    else
        filtered_dir="${CELLRANGER_OUT}/${sample}/outs/filtered_feature_bc_matrix"
        out_dir="${CELLRANGER_OUT}/${sample}/outs"
        raw_dir_args=(--raw-dir "${CELLRANGER_OUT}/${sample}/outs/raw_feature_bc_matrix")
    fi
    if [[ -d "${out_dir}/soupx_corrected_matrix" ]]; then
        plog "  Skipping $sample (already complete)"
        return 0
    fi
    plog "  Processing sample: $sample"
    conda run -n "$R_ENV" Rscript "${SCRIPT_DIR}/preprocess_qc.R" \
        --sample "$sample" \
        --filtered-dir "$filtered_dir" \
        "${raw_dir_args[@]}" \
        --out-dir "$out_dir" \
        >> "${CR_LOG_DIR}/${sample}_preprocess_qc.log" 2>&1
    plog "  Done: $sample"
}
export -f _run_preprocess_qc_one

run_preprocess_qc() {
    log "=== Stage 2: Ambient RNA correction (2a: SoupX) + doublet detection (2b: scDblFinder) ==="
    local jobs
    jobs=$(compute_parallel_jobs "$CELLRANGER_THREADS_PER_SAMPLE" "$CELLRANGER_MEM_GB_PER_SAMPLE")
    log "  Parallel preprocess_qc jobs: $jobs"
    export CELLRANGER_OUT PREPROCESS_QC_OUT MATRIX_DIR INPUT_MODE R_ENV SCRIPT_DIR CR_LOG_DIR
    printf '%s\n' "${SAMPLES[@]}" \
        | xargs -I{} -P "$jobs" bash -c 'set -euo pipefail; _run_preprocess_qc_one "$@"' _ {}
    log "Stage 2 complete."
}

# -- Stage 3-5: Python (Scanpy) ----------------------------------------------
run_scanpy() {
    log "=== Stages 3-5: QC -> Normalization -> Clustering (Scanpy) ==="

    local matrix_source
    if [[ "$INPUT_MODE" == "matrix" ]]; then
        matrix_source="$MATRIX_DIR"
    else
        matrix_source="$CELLRANGER_OUT"
    fi

    # Set BEFORE the process starts (not inside the script) so numpy/scipy's
    # BLAS backend picks it up regardless of when it gets imported relative to
    # argument parsing -- env vars set after import may not take effect.
    export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
    conda run -n "$PYTHON_ENV" python3 \
        "${SCRIPT_DIR}/scanpy_analysis.py" \
        --input-mode        "$INPUT_MODE" \
        --matrix-dir        "$matrix_source" \
        --preprocess-qc-dir "$PREPROCESS_QC_OUT" \
        --out-dir        "${PROJECT_DIR}/scanpy" \
        --viz-dir        "$SUMMARY_FIGURES_DIR" \
        --tables-dir     "$SUMMARY_TABLES_DIR" \
        --samples        "${SAMPLES[@]}" \
        --metadata-file  "$METADATA_FILE" \
        --min-genes      "$MIN_GENES" \
        --max-genes      "$MAX_GENES" \
        --max-mito       "$MAX_MITO" \
        --n-mads         "$N_MADS" \
        --n-mads-mt      "$N_MADS_MT" \
        --remove-doublets "$REMOVE_DOUBLETS" \
        --n-hvgs         "$N_HVGS" \
        --n-pcs          "$N_PCS" \
        --n-neighbors    "$N_NEIGHBORS" \
        --resolution     "$LEIDEN_RES" \
        --threads        "$THREADS" \
        2>&1 | tee -a "$LOG_FILE"
    verify_scanpy_outputs
    log "Stages 3-5 complete."
}

# -- Stage 6: Python (scCODA compositional differential abundance) ----------
# Independent of Stage 7 (R) -- reads only obs_metadata.csv, already written
# by Stage 3-5, so it doesn't wait on downstream.R. A missing CONTROL_CONDITION
# on a >2-level condition column aborts the whole pipeline here (via
# set -o pipefail), same hard-fail rule Stage 7a's DESeq2 enforces -- whichever
# stage runs first surfaces the error.
run_sccoda() {
    log "=== Stage 6: Compositional differential abundance (scCODA) ==="
    conda run -n "$PYTHON_ENV" python3 \
        "${SCRIPT_DIR}/differential_abundance.py" \
        --scanpy-dir        "${PROJECT_DIR}/scanpy" \
        --viz-dir           "$SUMMARY_FIGURES_DIR" \
        --tables-dir        "$SUMMARY_TABLES_DIR" \
        --control-condition "$CONTROL_CONDITION" \
        2>&1 | tee -a "$LOG_FILE"
    log "Stage 6 complete."
}

# -- Stage 7: R (downstream) -------------------------------------------------
run_r_downstream() {
    log "=== Stage 7: Downstream analysis (R) ==="
    export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
    conda run -n "$R_ENV" Rscript \
        "${SCRIPT_DIR}/downstream.R" \
        --scanpy-dir        "${PROJECT_DIR}/scanpy" \
        --out-dir           "${PROJECT_DIR}/downstream" \
        --viz-dir           "$SUMMARY_FIGURES_DIR" \
        --tables-dir        "$SUMMARY_TABLES_DIR" \
        --threads           "$THREADS" \
        --control-condition "$CONTROL_CONDITION" \
        2>&1 | tee -a "$LOG_FILE"
    verify_r_outputs
    log "Stage 7 complete."
}

# -- Environment setup (first run) -------------------------------------------
setup_envs() {
    log "Setting up conda environments: $PYTHON_ENV"
    conda env create -f "${SCRIPT_DIR}/envs/python_env.yaml" --name "$PYTHON_ENV" -q || true
    log "Setting up conda environments: $R_ENV"
    conda env create -f "${SCRIPT_DIR}/envs/r_env.yaml"      --name "$R_ENV"      -q || true
    log "Environments ready."
    verify_envs
}

# -- Main --------------------------------------------------------------------
main() {
    log "Pipeline started."
    log "Config:           $CONFIG"
    log "Log file:         $LOG_FILE"
    log "Project name:     $PROJECT_NAME"
    log "Project dir:      $PROJECT_DIR"
    log "Summary figures:  $SUMMARY_FIGURES_DIR"
    log "Summary tables:   $SUMMARY_TABLES_DIR"
    log "Input mode:       $INPUT_MODE"
    # IFS is $'\n\t' script-wide, so ${SAMPLES[*]} alone would join with
    # newlines here -- force a comma-space join for a readable one-line log.
    # (array-join with IFS only ever uses IFS's first character, hence the sed.)
    log "Samples (${#SAMPLES[@]}, from METADATA_FILE): $(IFS=,; echo "${SAMPLES[*]}" | sed 's/,/, /g')"
    [[ "$INPUT_MODE" == "fastq" ]] && log "Cell Ranger dir:  $CELLRANGER_OUT"
    [[ "$INPUT_MODE" == "matrix" ]] && log "Matrix dir:       $MATRIX_DIR"
    [[ "$SKIP_AMBIENT_QC" != "true" ]] && log "Preprocess QC dir: $PREPROCESS_QC_OUT"
    log "Control condition: ${CONTROL_CONDITION:-<not set>}"
    log "REMOVE_DOUBLETS=${REMOVE_DOUBLETS}"
    log "SETUP_ENVS=${SETUP_ENVS}"
    [[ "$INPUT_MODE" == "fastq" ]] && log "SKIP_CELLRANGER=${SKIP_CELLRANGER}"
    log "SKIP_AMBIENT_QC=${SKIP_AMBIENT_QC}"
    log "PREPROCESSING_ONLY=${PREPROCESSING_ONLY}"

    check_conda

    [[ "$SETUP_ENVS" == "true" ]] && setup_envs

    check_deps

    if [[ "$INPUT_MODE" == "matrix" ]]; then
        verify_matrix_inputs
    else
        [[ "$SKIP_CELLRANGER" != "true" ]] && run_cellranger
    fi
    [[ "$SKIP_AMBIENT_QC" != "true" ]] && run_preprocess_qc
    run_scanpy

    if [[ "$PREPROCESSING_ONLY" == "true" ]]; then
        log "PREPROCESSING_ONLY=true -- stopping after Stage 5 (Cell Ranger through clustering/annotation), before Stage 6 (scCODA) and Stage 7 (R downstream)."
        log "Project results:  $PROJECT_DIR"
        log "Summary figures:  $SUMMARY_FIGURES_DIR"
        log "Summary tables:   $SUMMARY_TABLES_DIR"
        log "Log file:         $LOG_FILE"
        return 0
    fi

    run_sccoda
    run_r_downstream

    log "Pipeline finished successfully."
    log "Project results:  $PROJECT_DIR"
    log "Summary figures:  $SUMMARY_FIGURES_DIR"
    log "Summary tables:   $SUMMARY_TABLES_DIR"
    log "Log file:         $LOG_FILE"
}

main "$@"

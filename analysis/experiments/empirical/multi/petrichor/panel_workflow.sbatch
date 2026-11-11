#!/bin/bash
# =====================================================================
# SLURM array task for the panel-conditioned diagnostic workflow.
# One array task = one line of jobs.tsv = one (panel, stage) unit.
# Fractional-node request (16 cores / modest mem / tight time) so tasks
# backfill onto `mix` nodes immediately instead of waiting for a whole node.
#
# Submit via submit_panel_workflow.sh (sets RUNDIR, NPERM, NREF, etc.).
# Expects in the environment / job script dir:
#   RUNDIR   self-contained run dir on the work FS
#   jobs.tsv with columns: panel<TAB>stage   (1 line per array task)
# =====================================================================
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=panelwf
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log

set -euo pipefail
module load R/4.4.1
export R_LIBS_USER=/path/to/work/Rlibs
# pin BLAS so process-parallel mclapply forks don't oversubscribe cores
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

RUNDIR="${RUNDIR:?set RUNDIR}"
cd "$RUNDIR"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" jobs.tsv)
PANEL=$(echo "$LINE" | cut -f1)
STAGE=$(echo "$LINE" | cut -f2)

export DATA_DIR="$RUNDIR/data"
export LIB_PANEL="$RUNDIR/lib_panel_workflow.R"
export LIB_ORTHO="$RUNDIR/lib_ortho.R"
export POLYDX_SRC="${POLYDX_SRC:-$RUNDIR/polyDx_R}"
export CV_CORES="${SLURM_CPUS_PER_TASK}"
# settings (overridable by the submit script via --export)
export PWF_NPERM="${PWF_NPERM:-30}"
export PWF_NREF="${PWF_NREF:-30}"
export PWF_SEED="${PWF_SEED:-2026}"
export PWF_MAXTR="${PWF_MAXTR:-8}"
export PWF_REPS="${PWF_REPS:-5}"
export PWF_USE_CACHE="${PWF_USE_CACHE:-1}"
export PWF_MARKER_QC="${PWF_MARKER_QC:-1}"
export PWF_MAF_MIN="${PWF_MAF_MIN:-0.05}"
export PWF_CALLRATE_MIN="${PWF_CALLRATE_MIN:-0.8}"

echo "[$(date '+%F %T')] task=${SLURM_ARRAY_TASK_ID} panel=${PANEL} stage=${STAGE} cores=${SLURM_CPUS_PER_TASK}"
Rscript "$RUNDIR/run_panel_workflow.R" "$PANEL" "$STAGE" "$RUNDIR/out"
echo "[$(date '+%F %T')] DONE task=${SLURM_ARRAY_TASK_ID}"

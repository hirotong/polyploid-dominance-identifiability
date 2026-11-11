#!/bin/bash
# Re-run the two large blueberry panels (rows 3=blueberry_gc, 4=blueberry2) that
# OOM-killed at 48G/32-forks. n~1800 kernels make each BGLR fork memory-heavy, so
# give a high-mem allocation and a moderate fork count. Numerics unchanged.
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=predcv
#SBATCH --partition=defq
#SBATCH --array=3,4
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --output=logs/predcv_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1
export CV_CORES=24
export R_LIBS_USER=/path/to/work/Rlibs
module load R/4.4.1

RUNDIR=/path/to/work/code/pred_cv_run
cd "$RUNDIR"
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" datasets.tsv)
key=$(echo "$row" | cut -f1); dos=$(echo "$row" | cut -f2); phe=$(echo "$row" | cut -f3); pl=$(echo "$row" | cut -f4)
echo "TASK ${SLURM_ARRAY_TASK_ID} key=$key ploidy=$pl host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=256G start=$(date)"
Rscript run_pred_cv.R "data/$dos" "data/$phe" "out/$key" "$pl" auto 5
echo "DONE $key $(date)"

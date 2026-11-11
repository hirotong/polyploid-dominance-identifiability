#!/bin/bash
# Cross-crop prediction-CV array: A vs A+D vs orthogonal(A+D2+..+Dm) per dataset.
# One array task per row of datasets.tsv. Submit from the run root (RUNDIR below).
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=predcv
#SBATCH --partition=defq
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/predcv_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1
export CV_CORES=8
export R_LIBS_USER=/path/to/work/Rlibs
module load R/4.4.1

RUNDIR=/path/to/work/code/pred_cv_run
cd "$RUNDIR"
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" datasets.tsv)
key=$(echo "$row" | cut -f1); dos=$(echo "$row" | cut -f2); phe=$(echo "$row" | cut -f3); pl=$(echo "$row" | cut -f4)
echo "TASK ${SLURM_ARRAY_TASK_ID} key=$key ploidy=$pl host=$(hostname) start=$(date)"
Rscript run_pred_cv.R "data/$dos" "data/$phe" "out/$key" "$pl" auto 5
echo "DONE $key $(date)"

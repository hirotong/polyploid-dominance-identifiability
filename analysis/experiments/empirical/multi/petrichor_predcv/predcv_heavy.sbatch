#!/bin/bash
# Re-run ONLY the heavy panels (rows 3=blueberry_gc, 4=blueberry2, 5=miscanthus of
# datasets.tsv) with the unit-parallel driver on full 64-core nodes. The other 5
# panels already finished with identical numerics and are not re-run.
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=predcv
#SBATCH --partition=defq
#SBATCH --array=3,4,5
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/predcv_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1
export CV_CORES=32
export R_LIBS_USER=/path/to/work/Rlibs
module load R/4.4.1

RUNDIR=/path/to/work/code/pred_cv_run
cd "$RUNDIR"
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" datasets.tsv)
key=$(echo "$row" | cut -f1); dos=$(echo "$row" | cut -f2); phe=$(echo "$row" | cut -f3); pl=$(echo "$row" | cut -f4)
echo "TASK ${SLURM_ARRAY_TASK_ID} key=$key ploidy=$pl host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} start=$(date)"
Rscript run_pred_cv.R "data/$dos" "data/$phe" "out/$key" "$pl" auto 5
echo "DONE $key $(date)"

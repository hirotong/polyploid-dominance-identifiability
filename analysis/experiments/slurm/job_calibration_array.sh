#!/bin/bash
# Petrichor SLURM ARRAY for the calibration-curve sweep — ONE (structure x fraction)
# cell per task (3 structures x 6 true-D2 fractions = 18 cells). Each task writes
# experiments/simulations/calibration/out/calib_cell_<NN>.csv; concatenate afterwards.
#
# Sized to BACKFILL: a fraction of a node (16 cores) + tight time/mem so tasks start
# on partially-used 'mix' nodes immediately instead of waiting for whole free nodes.
# PREREQ: an `r_env` conda env with AlphaSimR + BGLR must exist on petrichor (home
# dirs are NOT shared with krypton). Build once, e.g.:
#   module load miniforge3/24.9.0 && conda create -y -n r_env -c conda-forge \
#     r-base=4.4 r-alphasimr r-bglr r-data.table && \
#   # (AGHmatrix only needed for the optional Slater/Endelman supplement)
# The primary grid is small enough to run on krypton (run_primary_krypton.sh); this
# array is mainly for the HEAVY secondary stress grid (large n/M) — bump --mem/--time
# and point at the secondary script when that phase runs.
#SBATCH --job-name=calib_arr
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=defq
#SBATCH --array=1-18
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=experiments/simulations/calibration/calib_arr_%A_%a.log
source /etc/profile 2>/dev/null
module load miniforge3/24.9.0
cd "${CALIB_DIR:-/path/to/work/code/calibration}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NC=${SLURM_CPUS_PER_TASK:-16} REPS=${REPS:-100}
conda run -n r_env Rscript calibration_sim.R out
echo "CALIB_ARR cell ${SLURM_ARRAY_TASK_ID} EXIT $?"

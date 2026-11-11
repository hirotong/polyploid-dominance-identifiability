#!/bin/bash
# Petrichor driver for the SECONDARY calibration stress grid. One script,
# submitted once per N-tier with per-tier --array/--cpus-per-task/--mem/--time and
# --export=ALL,CALIB_TIER=..,REGIME_GRAIN=..,CHUNK=.. (plus optional CALIB_NSET,
# CALIB_MSET, CALIB_LDSET, CALIB_SSET, CALIB_ORTH, CALIB_ORTH_MAXN, REPS_* passed
# through the shell env -- NOT inline in --export, sbatch splits commas).
#
#   grain=cell    -> task = one (N,M,LD,s) cell; mclapply over reps; BLAS pinned 1
#   grain=cellrep -> task = one cell-rep; threaded BLAS (large-N eigendecomposition)
#
#SBATCH --job-name=calibsec
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/path/to/work/code/calibration_secondary/logs/calibsec_%x_%A_%a.log
source /etc/profile 2>/dev/null
module load miniforge3/24.9.0
cd /path/to/work/code/calibration_secondary
mkdir -p logs experiments/simulations/calibration/out
GRAIN=${REGIME_GRAIN:-cell}
if [ "$GRAIN" = "cellrep" ]; then
  export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} NC=1
else
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NC=${SLURM_CPUS_PER_TASK}
fi
export CHUNK=${CHUNK:-1} CALIB_OUTDIR=experiments/simulations/calibration/out
echo "CALIBSEC tier=${CALIB_TIER} grain=${GRAIN} arr=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} chunk=${CHUNK} $(date)"
mamba run -n r_env Rscript experiments/simulations/calibration/calibration_secondary.R
echo "CALIBSEC tier=${CALIB_TIER} arr=${SLURM_ARRAY_TASK_ID} rc=$? $(date)"

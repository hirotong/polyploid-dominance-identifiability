#!/bin/bash
# Generic driver for the EXPANDED regime sweep. One script, submitted once per
# N-tier with per-tier --array/--cpus-per-task/--mem/--time/--partition and
# --export=ALL,REGIME_TIER=..,REGIME_GRAIN=..,CHUNK=..  (mirrors job_calib.sh).
#
#   grain=cell    -> task = 1 (N,M,LD,dr) cell; mclapply over reps; BLAS pinned to 1
#   grain=cellrep -> task = 1 cell-rep; threaded BLAS (n^3 eigendecomp uses all cores)
#
#SBATCH --job-name=regexp
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/path/to/work/code/regime_expanded/logs/regexp_t%x_%A_%a.log
source /etc/profile 2>/dev/null
module load miniforge3/24.9.0
cd /path/to/work/code/regime_expanded
mkdir -p logs experiments/simulations/regime_map/out/expanded
GRAIN=${REGIME_GRAIN:-cell}
if [ "$GRAIN" = "cellrep" ]; then          # one rep/task: parallelize the eigendecomposition via BLAS
  export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} NC=1
else                                        # many reps/task via mclapply: pin BLAS, fork across cores
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NC=${SLURM_CPUS_PER_TASK}
fi
export CHUNK=${CHUNK:-1}
echo "REGEXP tier=${REGIME_TIER} grain=${GRAIN} arr=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} chunk=${CHUNK} $(date)"
mamba run -n r_env Rscript experiments/simulations/regime_map/regime_map_expanded.R
echo "REGEXP tier=${REGIME_TIER} arr=${SLURM_ARRAY_TASK_ID} rc=$? $(date)"

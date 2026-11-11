#!/bin/bash
# =====================================================================
# Petrichor SLURM array: updog dosage calling, 32 chunks x 16 cores = 512 cores.
# Submit:  sbatch updog_array.sbatch
# Reads/writes the shared work path (visible on krypton).
# =====================================================================
#SBATCH --job-name=bb_updog
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=defq
#SBATCH --array=1-32
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/work/code/bb_dosage/logs/chunk_%a.out
set -e
module load R/4.4.1
export R_LIBS_USER=/path/to/work/Rlibs
export WD=/path/to/work/code/bb_dosage
export NCHUNK=32
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript $WD/updog_chunk.R

#!/bin/bash
# Petrichor SLURM job: pure-architecture simulation behind Figure 1, with the
# alternative dominance kernels (Slater, Endelman) fitted on the same replicate
# traits (supplement table tab:robust). 3 structures x 5 truths x REPS replicates,
# four BGLR fits per replicate. BLAS pinned to one thread per mclapply worker.
# Outputs go to the gitignored runs/ tree; submit from the repository root after
#   mkdir -p runs/poly_ortho_rerun
#SBATCH --job-name=poly_ortho_alt
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/path/to/work/code/ar_alfalfa_haplotype_breeding/runs/poly_ortho_rerun/poly_ortho_%j.log
ROOT=/path/to/work
cd "$ROOT/code/ar_alfalfa_haplotype_breeding" || exit 1
export TMPDIR=/path/to/scratch/tmp; mkdir -p "$TMPDIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NC=${SLURM_CPUS_PER_TASK:-16} REPS=${REPS:-50}
OUT=${OUT:-runs/poly_ortho_rerun/out}
mkdir -p "$OUT"
git rev-parse HEAD > "$(dirname "$OUT")/CODE_COMMIT_${SLURM_JOB_ID:-local}"
echo "POLY_ORTHO_ALT start $(date) REPS=$REPS NC=$NC OUT=$OUT commit=$(git rev-parse --short HEAD)"
"$ROOT/conda/envs/r-shared/bin/Rscript" experiments/simulations/poly_ortho/poly_ortho_framework.R "$OUT"
echo "POLY_ORTHO_ALT EXIT $? $(date)"

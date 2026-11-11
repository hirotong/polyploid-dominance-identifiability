#!/bin/bash
# Re-run Miscanthus restricted to a single-site (Chuncheon) CORE-5 trait set, for
# BOTH the variance partition (ortho) and the prediction CV, so the Fig.4 Miscanthus
# point is internally consistent (size/x/y over the SAME 5 traits) and no longer
# mixes 5 environments / 113 near-duplicate trait-columns into a genotype-only test.
#
# Core-5 = the traits with complete cross-site coverage in the source dataset
# (Njuguna et al., IDB-6439486), taken at the primary Chuncheon trial:
#   Biomass.Yield.dry.weight, Culm.Length, Culm.per.area, HD1, Total.Culm
#
# Writes to the out/miscanthus_core5* prefix (does NOT clobber the committed
# 113-trait _PA.csv / old 5-trait _ortho.csv until verified + integrated locally).
# Run dir = pred_cv_run on the shared work path (has data/, run_pred_cv.R,
# experiments/simulations/poly_ortho/lib_ortho{,_mp}.R). run_ortho_dataset.R is copied in
# alongside (same relative source path for lib_ortho.R).
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=misc_core5
#SBATCH --partition=defq
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/misc_core5_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1
export CV_CORES=${SLURM_CPUS_PER_TASK:-16}
export R_LIBS_USER=/path/to/work/Rlibs
module load R/4.4.1

RUNDIR=/path/to/work/code/pred_cv_run
cd "$RUNDIR"

TRAITS="Chuncheon.Biomass.Yield.dry.weight,Chuncheon.Culm.Length,Chuncheon.Culm.per.area,Chuncheon.HD1,Chuncheon.Total.Culm"
echo "HOST=$(hostname) CPUS=${SLURM_CPUS_PER_TASK} START=$(date)"
echo "TRAITS=$TRAITS"

# 1) variance partition + diagnostics + 30-perm null (do_cv=0); perm-null is on TR[1]=Biomass.Yield
echo "=== ORTHO (variance partition) ==="
Rscript run_ortho_dataset.R data/miscanthus_dosage.rds data/miscanthus_pheno.rds out/miscanthus_core5 "$TRAITS" 0

# 2) paired 5-fold x 5-rep prediction CV (A / A+D / orthogonal), ploidy 4
echo "=== PREDICTION CV ==="
Rscript run_pred_cv.R data/miscanthus_dosage.rds data/miscanthus_pheno.rds out/miscanthus_core5 4 "$TRAITS" 5

echo "DONE $(date)"

#!/bin/bash
# =====================================================================
# Submit the estimator / prior sensitivity run (estimator_sensitivity.R) as one
# SLURM array, one task per panel, from a self-contained run dir on the work FS.
#
# Run on the petrichor LOGIN node:
#   bash submit_estimator_sensitivity.sh [full|smoke]
#     full  -> all 11 empirical panels + the simulated positive control
#     smoke -> rose only, short chains, 3 permutations
#   PANELS_OVERRIDE="rose potato" overrides the panel list.
# Outputs: $RUNDIR/out/<panel>__estsens.csv, <panel>__remlprofile.csv
# =====================================================================
set -euo pipefail
MODE="${1:-full}"
WORK=/path/to/work/code
PROJECT="${PROJECT:-$WORK/ar_alfalfa_haplotype_breeding}"
SRC_DATA="${SRC_DATA:-$PROJECT/data/processed/panel_inputs}"
RUNDIR="${RUNDIR:-$PROJECT/runs/estimator_sensitivity}"
mkdir -p "$RUNDIR/data" "$RUNDIR/out" "$RUNDIR/polyDx_R"

PANELS=(alfalfa potato potato2 blueberry blueberry2 blueberry_gc panicum rose miscanthus sweetpotato sugarcane poscontrol)
SMOKE=0; TIME=04:00:00
if [ "$MODE" = "smoke" ]; then PANELS=(rose); SMOKE=1; TIME=00:20:00; fi
if [ -n "${PANELS_OVERRIDE:-}" ]; then read -r -a PANELS <<< "$PANELS_OVERRIDE"; fi

for p in "${PANELS[@]}"; do
  for suf in dosage pheno; do cp -f "$SRC_DATA/${p}_${suf}.rds" "$RUNDIR/data/"; done
done
cp -f "$SRC_DATA/sugarcane_marker_ploidy.rds" "$RUNDIR/data/" 2>/dev/null || true
cp -f "$PROJECT/experiments/empirical/multi/estimator_sensitivity.R" "$RUNDIR/"
cp -f "$PROJECT/experiments/empirical/multi/lib_reml.R"              "$RUNDIR/"
cp -f "$PROJECT/experiments/empirical/multi/lib_panel_workflow.R"    "$RUNDIR/"
cp -f "$PROJECT/experiments/simulations/poly_ortho/lib_ortho.R"      "$RUNDIR/"
rm -f "$RUNDIR/polyDx_R/"*.R
cp -f "$PROJECT/software/polyDx/R/"*.R "$RUNDIR/polyDx_R/"
( cd "$PROJECT" && git rev-parse HEAD ) > "$RUNDIR/CODE_COMMIT"
md5sum "$RUNDIR"/data/*.rds > "$RUNDIR/DATA_MD5"

printf "%s\n" "${PANELS[@]}" > "$RUNDIR/panels.txt"
N=$(wc -l < "$RUNDIR/panels.txt")
RLIB=/path/to/work/Rlibs
cd "$RUNDIR"
JID=$(sbatch --parsable --account=YOUR_ACCOUNT --job-name=estsens \
  --cpus-per-task=16 --mem=64G --time="$TIME" --output="$RUNDIR/estsens_%A_%a.log" --array=1-"$N" \
  --wrap="module load R/4.4.1; export R_LIBS_USER=$RLIB OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 DATA_DIR=$RUNDIR/data LIB_PANEL=$RUNDIR/lib_panel_workflow.R LIB_ORTHO=$RUNDIR/lib_ortho.R LIB_REML=$RUNDIR/lib_reml.R POLYDX_SRC=$RUNDIR/polyDx_R PWF_MARKER_QC=1 PWF_MAF_MIN=0.05 PWF_CALLRATE_MIN=0.8 PWF_SEED=2026 EST_SMOKE=$SMOKE; P=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $RUNDIR/panels.txt); Rscript $RUNDIR/estimator_sensitivity.R \$P $RUNDIR/out")
echo "estsens array: $JID ($N panels, mode=$MODE, rundir=$RUNDIR)"
echo "$JID" > "$RUNDIR/.jobid"

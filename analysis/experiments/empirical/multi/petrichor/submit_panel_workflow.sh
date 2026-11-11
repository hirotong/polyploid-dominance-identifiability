#!/bin/bash
# =====================================================================
# Submit the panel-conditioned diagnostic workflow as two SLURM arrays
# from the self-contained run dir on the work FS (code already staged there):
#   (1) build  -- cache each panel's kernels         (afterok dependency for (2))
#   (2) main   -- one task per (panel, stage) of jobs.tsv
#
# Prereq (done once by the operator): RUNDIR holds the code files
#   lib_ortho.R lib_ortho_mp.R lib_panel_workflow.R run_panel_workflow.R panel_workflow.sbatch
#   polyDx_R/ (software/polyDx/R/*.R, staged as plain sources -- see below)
#
# Run on the petrichor LOGIN node:
#   bash submit_panel_workflow.sh [full|affected|smoke]
#     full  -> NPERM=NREF=30, all 11 panels  (default)
#     affected -> NPERM=NREF=30, panels whose marker set changes under QC
#     smoke -> NPERM=NREF=10, rose,sugarcane only
#   PANELS_OVERRIDE="poscontrol" RUNDIR=... bash submit_panel_workflow.sh full
#     -> the simulated positive control only (build its inputs first with
#        build_poscontrol_panel.R <lib_ortho.R> $SRC_DATA)
# Outputs: $RUNDIR/out/<panel>__<stage>.csv
# =====================================================================
set -euo pipefail
MODE="${1:-full}"

WORK=/path/to/work/code
PROJECT="${PROJECT:-$WORK/ar_alfalfa_haplotype_breeding}"
SRC_DATA="${SRC_DATA:-$PROJECT/data/processed/panel_inputs}"
RUNDIR="${RUNDIR:-$PROJECT/runs/panel_cond_workflow}"
mkdir -p "$RUNDIR/data" "$RUNDIR/out"

# --- data: copy all panels' dosage/pheno (+group, +marker_ploidy) -----------
ALL_PANELS=(alfalfa potato potato2 blueberry blueberry2 blueberry_gc panicum rose miscanthus sweetpotato sugarcane)
for p in "${ALL_PANELS[@]}"; do
  for suf in dosage pheno group; do
    f="$SRC_DATA/${p}_${suf}.rds"; [ -f "$f" ] && cp -n "$f" "$RUNDIR/data/" || true
  done
done
cp -n "$SRC_DATA/sugarcane_marker_ploidy.rds" "$RUNDIR/data/" 2>/dev/null || true

# --- code: stage run scripts + libs into RUNDIR (self-contained) ------------
cp -f "$PROJECT/experiments/empirical/multi/run_panel_workflow.R"            "$RUNDIR/"
cp -f "$PROJECT/experiments/empirical/multi/lib_panel_workflow.R"            "$RUNDIR/"
cp -f "$PROJECT/experiments/empirical/multi/petrichor/panel_workflow.sbatch" "$RUNDIR/"
cp -f "$PROJECT/experiments/simulations/poly_ortho/lib_ortho.R"              "$RUNDIR/"
cp -f "$PROJECT/experiments/simulations/poly_ortho/lib_ortho_mp.R"           "$RUNDIR/"

# polyDx package sources (the workflow implementation). Staged as plain sources so
# the run dir stays self-contained and needs no R CMD INSTALL on the compute node.
# Clear stale sources first so a file removed from software/polyDx/R/ since the
# last stage doesn't linger in $RUNDIR/polyDx_R and get sourced anyway.
mkdir -p "$RUNDIR/polyDx_R"
rm -f "$RUNDIR/polyDx_R/"*.R
cp -f "$PROJECT/software/polyDx/R/"*.R "$RUNDIR/polyDx_R/"

# --- mode-specific settings + panel set -------------------------------------
if [ "$MODE" = "smoke" ]; then
  NPERM=10; NREF=10; MAXTR=4; REPS=3; PANELS=(rose sugarcane)
elif [ "$MODE" = "affected" ]; then
  NPERM=30; NREF=30; MAXTR=8; REPS=5; PANELS=(alfalfa potato blueberry_gc rose miscanthus sugarcane)
else
  NPERM=30; NREF=30; MAXTR=8; REPS=5; PANELS=("${ALL_PANELS[@]}")
fi
# optional explicit panel subset (space-separated), e.g. PANELS_OVERRIDE="potato2 sweetpotato"
if [ -n "${PANELS_OVERRIDE:-}" ]; then read -r -a PANELS <<< "$PANELS_OVERRIDE"; fi
# copy inputs for selected panels outside ALL_PANELS (e.g. poscontrol, whose inputs
# are written into $SRC_DATA by build_poscontrol_panel.R); fail if a dosage file is missing
for p in "${PANELS[@]}"; do
  for suf in dosage pheno group; do
    f="$SRC_DATA/${p}_${suf}.rds"; [ -f "$f" ] && cp -n "$f" "$RUNDIR/data/" || true
  done
  [ -f "$RUNDIR/data/${p}_dosage.rds" ] || { echo "ERROR: no dosage file for panel '$p' in $SRC_DATA" >&2; exit 1; }
done
SEED=2026
MARKER_QC="${PWF_MARKER_QC:-1}"
MAF_MIN="${PWF_MAF_MIN:-0.05}"
CALL_MIN="${PWF_CALLRATE_MIN:-0.8}"

# --- job lists --------------------------------------------------------------
STAGES=(geom real perm addctrl d2spike structspike pred)
: > "$RUNDIR/jobs.tsv"; : > "$RUNDIR/build_jobs.tsv"
for p in "${PANELS[@]}"; do
  printf "%s\tbuild\n" "$p" >> "$RUNDIR/build_jobs.tsv"
  for s in "${STAGES[@]}"; do printf "%s\t%s\n" "$p" "$s" >> "$RUNDIR/jobs.tsv"; done
done
NBUILD=$(wc -l < "$RUNDIR/build_jobs.tsv"); NJOBS=$(wc -l < "$RUNDIR/jobs.tsv")

EXP="RUNDIR=$RUNDIR,POLYDX_SRC=$RUNDIR/polyDx_R,PWF_NPERM=$NPERM,PWF_NREF=$NREF,PWF_SEED=$SEED,PWF_MAXTR=$MAXTR,PWF_REPS=$REPS,PWF_USE_CACHE=1,PWF_MARKER_QC=$MARKER_QC,PWF_MAF_MIN=$MAF_MIN,PWF_CALLRATE_MIN=$CALL_MIN"
RLIB=/path/to/work/Rlibs
echo "RUNDIR=$RUNDIR  MODE=$MODE  panels=${#PANELS[@]}  build=$NBUILD  main=$NJOBS  (NPERM=$NPERM NREF=$NREF MAF=$MAF_MIN CALL=$CALL_MIN QC=$MARKER_QC)"
cd "$RUNDIR"

# (1) build-cache array: 8 cores, 40 min, generous mem (miscanthus 268k markers)
JID_BUILD=$(sbatch --parsable --account=YOUR_ACCOUNT --job-name=panelbuild \
  --cpus-per-task=8 --mem=72G --time=00:40:00 --output="$RUNDIR/build_%A_%a.log" --array=1-"$NBUILD" \
  --export=ALL,"$EXP" \
  --wrap="module load R/4.4.1; export R_LIBS_USER=$RLIB OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 DATA_DIR=$RUNDIR/data LIB_PANEL=$RUNDIR/lib_panel_workflow.R LIB_ORTHO=$RUNDIR/lib_ortho.R POLYDX_SRC=$RUNDIR/polyDx_R; L=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $RUNDIR/build_jobs.tsv); Rscript $RUNDIR/run_panel_workflow.R \$(echo \"\$L\"|cut -f1) build $RUNDIR/out")
echo "build array: $JID_BUILD ($NBUILD tasks)"

# (2) main array starts after build FINISHES (afterany: a failed build just means
#     that panel's main tasks rebuild the cache themselves, so don't block on it)
JID_MAIN=$(sbatch --parsable --dependency=afterany:"$JID_BUILD" \
  --array=1-"$NJOBS" --export=ALL,"$EXP" "$RUNDIR/panel_workflow.sbatch")
echo "main array:  $JID_MAIN ($NJOBS tasks, afterany:$JID_BUILD)"
echo "$JID_BUILD $JID_MAIN" > "$RUNDIR/.jobids"
echo "monitor: squeue -u $USER ;  results: $RUNDIR/out/<panel>__<stage>.csv"

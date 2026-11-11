#!/usr/bin/env Rscript
# =====================================================================
# run_panel_workflow.R  --  one (panel, stage) unit of the panel-conditioned
# diagnostic workflow.  Decomposed by stage so the whole battery runs as a wide
# SLURM array (one task per panel x stage), each task a fraction of a node.
#
# Usage:
#   Rscript run_panel_workflow.R <panel_key> <stage> <out_dir>
#
# Stages:
#   build        cache panel kernels (A, D, O[1..ndeg]) + geometry -> <panel>_KP.rds
#   geom         kernel-geometry diagnostics row
#   real         real-trait variance fits (standard A+D + orthogonal degrees)
#   perm         phenotype-permutation null (representative trait[s])
#   addctrl      panel-conditioned additive-only control  (true D2 = 0)
#   d2spike      panel-conditioned D2 spike-in positive control (s in s_spike)
#   structspike  panel-conditioned structure-spike control (PC1 / family block)
#   pred         paired 5-fold x reps prediction CV (A, A+D, orthogonal)
#
# Env:
#   DATA_DIR        dir holding <panel>_dosage.rds / _pheno.rds / _group.rds
#   PWF_NPERM=30    permutations per trait        PWF_NREF=30  sim reps per cell
#   PWF_SEED=2026   global seed                   PWF_MAXTR=8  max real traits
#   PWF_PERM_TRAITS=3   traits to run the null on PWF_REPS=5   CV repeats
#   PWF_PERM_ORTHO=0    also fit orthogonal in perm null (P_D2_perm/L_D2_perm)
#   PWF_H2          override simulated heritability (else 0.6)
#   CV_CORES / SLURM_CPUS_PER_TASK   parallel workers
#   PWF_FRAC_NITER/BURN/THIN, PWF_PT_NITER/BURN   BGLR controls (see lib)
# Output: <out_dir>/<panel>__<stage>.csv   (concatenated by the aggregator)
#
# The panel-conditioned diagnostic workflow itself now lives in the R package
# polyDx (software/polyDx); this driver loads it via lib_panel_workflow.R (a
# thin bootstrap shim) and delegates every computational stage to
# polyDx::run_workflow_stage(), keeping only CLI plumbing here: argument
# parsing, panel loading, kernel caching, and CSV writing.
# =====================================================================
LIB <- Sys.getenv("LIB_PANEL", "experiments/empirical/multi/lib_panel_workflow.R")
source(LIB)
suppressMessages(library(parallel))

args <- commandArgs(TRUE)
if (length(args) < 3) stop("usage: run_panel_workflow.R <panel_key> <stage> <out_dir>", call. = FALSE)
panel_key <- args[1]; stage <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DATA_DIR    <- Sys.getenv("DATA_DIR", "data/processed/panel_inputs")
NPERM       <- as.integer(Sys.getenv("PWF_NPERM", "30"))
NREF        <- as.integer(Sys.getenv("PWF_NREF", "30"))
SEED        <- as.integer(Sys.getenv("PWF_SEED", "2026"))
MAXTR       <- as.integer(Sys.getenv("PWF_MAXTR", "8"))
PERM_TRAITS <- as.integer(Sys.getenv("PWF_PERM_TRAITS", "3"))
PERM_ORTHO  <- as.integer(Sys.getenv("PWF_PERM_ORTHO", "0")) == 1L
CV_REPS     <- as.integer(Sys.getenv("PWF_REPS", "5"))
CORES       <- as.integer(Sys.getenv("CV_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
H2          <- suppressWarnings(as.numeric(Sys.getenv("PWF_H2", "")))
if (!is.finite(H2)) H2 <- PWF_H2_DEFAULT
S_SPIKE     <- c(0.25, 0.50, 0.75)
LAMBDA      <- c(0, 0.25, 0.50, 0.75, 1.00)

# ---- panel registry: key -> display name, ploidy spec, group file ----------
REG <- panel_registry()
# simulated positive control (inputs written by build_poscontrol_panel.R); it is
# not an empirical panel, so it is registered here rather than in polyDx
REG$poscontrol <- list(name = "Positive control (sim)", ploidy = "4", group = "none")
if (is.null(REG[[panel_key]])) stop(sprintf("unknown panel '%s'", panel_key), call. = FALSE)
meta <- REG[[panel_key]]
PANEL <- meta$name
ploidy_arg <- if (grepl("%s", meta$ploidy)) sprintf(meta$ploidy, DATA_DIR) else meta$ploidy
group_rds  <- if (identical(meta$group, "none")) "none" else file.path(DATA_DIR, meta$group)
dosage_rds <- file.path(DATA_DIR, paste0(panel_key, "_dosage.rds"))
pheno_rds  <- file.path(DATA_DIR, paste0(panel_key, "_pheno.rds"))
kp_cache   <- file.path(out_dir, paste0(panel_key, "_KP.rds"))
out_csv    <- file.path(out_dir, sprintf("%s__%s.csv", panel_key, stage))
qc_csv     <- file.path(out_dir, sprintf("%s__qcaudit.csv", panel_key))

cat(sprintf("[%s | %s] data=%s ploidy=%s cores=%d  H2=%.2f  nperm=%d nref=%d\n",
            PANEL, stage, DATA_DIR, ploidy_arg, CORES, H2, NPERM, NREF))

P <- load_panel(dosage_rds, pheno_rds, ploidy_arg, group_rds, max_traits = MAXTR)
cat(sprintf("  n=%d markers=%d traits=%d group=%s ploidy=%s\n",
            P$n, P$markers, length(P$traits), if (is.null(P$group)) "none" else paste0(length(unique(P$group)), "grp"), P$pcfg$label))
cat(sprintf("  QC: raw_markers=%d pre_qc=%d post_qc=%d dropped=%d MAF>=%.3f callrate>%.2f\n",
            P$qc$raw_markers, P$qc$markers_pre_qc, P$qc$markers_post_qc,
            P$qc$markers_dropped_qc, P$qc$maf_min, P$qc$callrate_min))
if (stage == "build" || !file.exists(qc_csv)) {
  qcout <- cbind(data.frame(panel = PANEL, panel_key = panel_key, ploidy = P$pcfg$label,
                            stringsAsFactors = FALSE), P$qc)
  tmp <- paste0(qc_csv, ".tmp", Sys.getpid())
  write.csv(qcout, tmp, row.names = FALSE)
  file.rename(tmp, qc_csv)
}

# ---- kernel cache: optional reuse across stages (atomic write, race-safe) --
USE_CACHE <- Sys.getenv("PWF_USE_CACHE", "0") == "1"
get_KP <- function(need_X = FALSE, force_write = FALSE) {
  KP <- NULL
  if (USE_CACHE && file.exists(kp_cache)) KP <- tryCatch(readRDS(kp_cache), error = function(e) NULL)
  if (is.null(KP)) {
    KP <- build_panel(P$M, P$pcfg, need_X = FALSE)
    if (USE_CACHE || force_write) { tmp <- paste0(kp_cache, ".tmp", Sys.getpid()); saveRDS(KP, tmp); file.rename(tmp, kp_cache) }
  }
  if (need_X) {                                    # add degree-1/2 contrasts for simulators
    cc <- build_contrasts12(P$M, P$pcfg); KP$X <- cc$X
    stopifnot(cc$valid == KP$valid)
  }
  KP
}

set.seed(SEED)

# =====================================================================
# STAGE: build  (cache kernels only)
# =====================================================================
if (stage == "build") {
  KP <- get_KP(need_X = FALSE, force_write = TRUE)
  ids <- rownames(P$M); check_kernels(KP, ids)
  geo <- kernel_geometry(KP)
  cat(sprintf("  cached kernels: valid=%d AD_cos=%.3f eff_rank=%.2f cond=%.1f\n",
              KP$valid, geo$AD_cosine, geo$eff_rank, geo$condition_number))
  cat("BUILD DONE\n"); quit(save = "no")
}

# =====================================================================
# STAGE: geom  (kernel-geometry diagnostics row + full degree cosine matrix)
# =====================================================================
if (stage == "geom") {
  KP <- get_KP()
  res <- run_workflow_stage("geom", P, KP, panel_name = PANEL, panel_key = panel_key,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  write.csv(res$extra$geom_cosine,
            file.path(out_dir, sprintf("%s__geom_cosine.csv", panel_key)), row.names = FALSE)
  cat("GEOM DONE\n"); quit(save = "no")
}

# =====================================================================
# STAGE: real  (real-trait variance fits)
# =====================================================================
if (stage == "real") {
  if (!length(P$traits)) { cat("no traits; skipping\n"); write.csv(data.frame(), out_csv, row.names = FALSE); quit(save = "no") }
  KP <- get_KP()
  res <- run_workflow_stage("real", P, KP, panel_name = PANEL, panel_key = panel_key,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  cat(sprintf("REAL DONE (%d/%d traits)\n", if (is.null(res$table)) 0 else nrow(res$table), length(P$traits))); quit(save = "no")
}

# =====================================================================
# STAGE: perm  (phenotype-permutation null on representative trait[s])
# =====================================================================
if (stage == "perm") {
  if (!length(P$traits)) { cat("no traits; skipping\n"); write.csv(data.frame(), out_csv, row.names = FALSE); quit(save = "no") }
  KP <- get_KP(need_X = FALSE)
  res <- run_workflow_stage("perm", P, KP, panel_name = PANEL, panel_key = panel_key,
                            nperm = NPERM, perm_traits = PERM_TRAITS, perm_ortho = PERM_ORTHO,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  use_tr <- P$traits[seq_len(min(PERM_TRAITS, length(P$traits)))]
  cat(sprintf("PERM DONE (%d traits x %d perms)\n", length(use_tr), NPERM)); quit(save = "no")
}

# =====================================================================
# STAGE: addctrl  (additive-only control; true D2 = 0)
# =====================================================================
if (stage == "addctrl") {
  KP <- get_KP(need_X = TRUE)
  res <- run_workflow_stage("addctrl", P, KP, panel_name = PANEL, panel_key = panel_key,
                            nref = NREF, h2 = H2, s_spike = S_SPIKE, lambda = LAMBDA,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  cat(sprintf("ADDCTRL DONE (%d/%d reps)\n", if (is.null(res$table)) 0 else nrow(res$table), NREF)); quit(save = "no")
}

# =====================================================================
# STAGE: d2spike  (D2 spike-in positive control)
# =====================================================================
if (stage == "d2spike") {
  KP <- get_KP(need_X = TRUE)
  res <- run_workflow_stage("d2spike", P, KP, panel_name = PANEL, panel_key = panel_key,
                            nref = NREF, h2 = H2, s_spike = S_SPIKE, lambda = LAMBDA,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  cat(sprintf("D2SPIKE DONE (%d/%d cells)\n", if (is.null(res$table)) 0 else nrow(res$table), length(S_SPIKE) * NREF)); quit(save = "no")
}

# =====================================================================
# STAGE: structspike  (structure-driven false-signal control)
# =====================================================================
if (stage == "structspike") {
  KP <- get_KP(need_X = TRUE)
  res <- run_workflow_stage("structspike", P, KP, panel_name = PANEL, panel_key = panel_key,
                            nref = NREF, h2 = H2, s_spike = S_SPIKE, lambda = LAMBDA,
                            seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  n_stypes <- length(c("PC", if (!is.null(P$group)) "family"))
  cat(sprintf("STRUCTSPIKE DONE (%d/%d cells)\n", if (is.null(res$table)) 0 else nrow(res$table),
              n_stypes * length(LAMBDA) * NREF)); quit(save = "no")
}

# =====================================================================
# STAGE: pred  (paired 5-fold x reps prediction CV; A, A+D, orthogonal)
# =====================================================================
if (stage == "pred") {
  if (!length(P$traits)) { cat("no traits; skipping\n"); write.csv(data.frame(), out_csv, row.names = FALSE); quit(save = "no") }
  KP <- get_KP(need_X = FALSE)
  res <- run_workflow_stage("pred", P, KP, panel_name = PANEL, panel_key = panel_key,
                            cv_reps = CV_REPS, seed = SEED, cores = CORES)
  write.csv(res$table, out_csv, row.names = FALSE)
  n_traits <- if (is.null(res$table)) 0 else length(unique(res$table$trait))
  cat(sprintf("PRED DONE (%d traits x %d reps)\n", n_traits, CV_REPS)); quit(save = "no")
}

stop(sprintf("unknown stage '%s'", stage), call. = FALSE)

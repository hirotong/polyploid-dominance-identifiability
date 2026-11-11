#!/usr/bin/env Rscript
# =====================================================================
# True-digenic CALIBRATION-CURVE simulation (known truth).
# Question: as the TRUE digenic fraction s = Var(G_D2)/(Var(G_A)+Var(G_D2))
# increases from 0 to 0.75, does the fitted STANDARD A+D dominance share
#   std_D_share = Vd / (Va + Vd)
# track it (identity line)? And does the ORTHOGONAL A+D2+D3+D4 model's
# fitted D2 fraction track it better?
#
# Reuses the EXACT Figure-1 known-truth pipeline (poly_ortho_framework.R):
#   - one shared 450-founder runMacs pool -> three structures
#     (unrelated / half-sib / full-sib; n~400 each), same ~1100-marker panel,
#     same 200 causal loci polymorphic in all three structures;
#   - same orthonormal dosage contrasts, standard (VanRaden A, digenic D) and
#     orthogonal (A^,D2^,D3^,D4^) kernels, same BGLR RKHS fits (4000/1000).
# The ONLY new piece: mixed additive+digenic traits with a controlled true
# digenic fraction (instead of the pure-degree truths used in Figure 1):
#   G = sqrt(1-s) * G_A* + sqrt(s) * G_D2*,  with G_A*,G_D2* standardized to
#   unit variance within the panel (degree-1 and degree-2 contrasts resp.),
#   then residual noise added to target H2 = 0.6.
#
# Models fit per replicate: (1) A only, (2) standard A+D, (3) orthogonal
# A+D2+D3+D4. Output: one LONG-format CSV row per (replicate x model) with
# all columns required by the experiment plan.
#
# Run modes:
#   - SLURM array: SLURM_ARRAY_TASK_ID in 1..18 selects ONE (structure x
#     fraction) cell; writes out/calib_cell_<NN>.csv.
#   - Local/legacy: no array id -> runs every cell into out/calibration_primary.csv
#     (set DRYRUN=1 for a fast 1-cell, REPS=2 smoke test).
# =====================================================================
suppressMessages({library(AlphaSimR); library(BGLR); library(parallel)})

OUT  <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "experiments/simulations/calibration/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PL   <- 4L
H2   <- 0.6
NITER <- 4000L; BURN <- 1000L
REPS <- as.integer(Sys.getenv("REPS", "100"))
NC   <- as.integer(Sys.getenv("NC", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
FRACTIONS  <- c(0, 0.05, 0.10, 0.25, 0.50, 0.75)   # target true D2 fraction s
STRUCTURES <- c("unrelated", "halfsib", "fullsib")

## ---- orthonormal dosage contrasts under binomial genotype freqs (Fig-1 code) ----
ortho <- function(p, ploidy = PL) {
  d <- 0:ploidy; w <- choose(ploidy, d) * p^d * (1 - p)^(ploidy - d); w <- w / sum(w)
  inn <- function(a, b) sum(w * a * b); Q <- matrix(0, ploidy + 1, ploidy + 1)
  for (k in 1:(ploidy + 1)) { v <- d^(k - 1)
    if (k > 1) for (j in 1:(k - 1)) v <- v - inn(v, Q[, j]) / inn(Q[, j], Q[, j]) * Q[, j]
    nrm <- sqrt(inn(v, v)); Q[, k] <- if (nrm > 0) v / nrm else v }
  Q }  # cols 1..5 = degree 0..4 ; rows d=0..4

## standard kernels: VanRaden additive, digenic (heterozygosity) dominance
Amat <- function(M, pl = PL) { p <- colMeans(M) / pl; Z <- sweep(M, 2, pl * p); d <- sum(pl * p * (1 - p)); if (d <= 0) d <- 1
  K <- tcrossprod(Z) / d; K / mean(diag(K)) + diag(1e-6, nrow(M)) }
Dmat <- function(M, pl = PL) { Hh <- (M * (pl - M)) / choose(pl, 2); Hh <- sweep(Hh, 2, colMeans(Hh))
  K <- tcrossprod(Hh); K / mean(diag(K)) + diag(1e-6, nrow(M)) }

## orthogonal kernels A^,D2^,D3^,D4^ from per-marker degree contrasts
build_ortho_kernels <- function(M, ploidy = PL) {
  n <- nrow(M); poly <- which(apply(M, 2, function(g) { p <- mean(round(g)) / ploidy; p > 0 && p < 1 }))
  X <- list(matrix(0, n, length(poly)), matrix(0, n, length(poly)), matrix(0, n, length(poly)), matrix(0, n, length(poly)))
  for (idx in seq_along(poly)) { j <- poly[idx]; g <- round(M[, j]); p <- mean(g) / ploidy; Q <- ortho(p, ploidy)
    for (deg in 1:4) X[[deg]][, idx] <- Q[g + 1, deg + 1] }
  mkK <- function(Z) { K <- tcrossprod(Z) / ncol(Z); K / mean(diag(K)) + diag(1e-6, nrow(Z)) }
  list(A = mkK(X[[1]]), D2 = mkK(X[[2]]), D3 = mkK(X[[3]]), D4 = mkK(X[[4]])) }

## collinearity diagnostics: A-D Frobenius cosine + 4-kernel effective rank
dc   <- function(K) { n <- nrow(K); J <- diag(n) - 1 / n; J %*% K %*% J }
fcos <- function(A, B) { a <- dc(A); b <- dc(B); sum(a * b) / sqrt(sum(a * a) * sum(b * b)) }
design_diag <- function(Klist) {
  vc <- sapply(Klist, function(K) { v <- as.vector(dc(K)); v / sqrt(sum(v * v)) })
  G <- crossprod(vc); ev <- eigen(G, only.values = TRUE)$values; ev <- pmax(ev, 0)
  list(cond = max(ev) / max(min(ev), 1e-12), eff_rank = sum(ev)^2 / sum(ev^2)) }

## one standardized degree-`deg` genetic value from the causal loci (var 1 in panel)
build_comp <- function(M, causal, deg) {
  gv <- numeric(nrow(M))
  for (j in causal) { g <- round(M[, j]); p <- mean(g) / PL; if (p <= 0 || p >= 1) next
    Q <- ortho(p, PL); gv <- gv + rnorm(1) * Q[g + 1, deg + 1] }
  gv }

## ---- ONE shared coalescent founder pool -> three structures (identical to Fig 1) ----
build_panel <- function() {
  set.seed(123)
  N0 <- 450L
  fp <- runMacs(nInd = N0, nChr = 8, segSites = 600, ploidy = PL, species = "GENERIC")
  SP <- SimParam$new(fp); SP$addSnpChip(nSnpPerChr = 300L)
  founders <- newPop(fp, simParam = SP)
  pop_unrel <- founders[1:400]
  pop_fs    <- randCross(founders[1:40], nCrosses = 50, nProgeny = 8, simParam = SP)
  sires <- 1:50; damsM <- matrix(51:450, nrow = 50, byrow = TRUE)
  plan <- do.call(rbind, lapply(1:50, function(f) cbind(damsM[f, ], rep(sires[f], 8))))
  pop_hs <- makeCross(founders, crossPlan = plan, nProgeny = 1L, simParam = SP)
  POPS <- list(unrelated = pullSnpGeno(pop_unrel, simParam = SP),
               halfsib   = pullSnpGeno(pop_hs,   simParam = SP),
               fullsib   = pullSnpGeno(pop_fs,   simParam = SP))[STRUCTURES]
  polyOf <- function(M) which(apply(M, 2, function(g) { p <- mean(round(g)) / PL; p > 0 && p < 1 }))
  common <- Reduce(intersect, lapply(POPS, polyOf))
  set.seed(7); causal <- sort(sample(common, 200))
  list(POPS = POPS, causal = causal, nmark = ncol(POPS[[1]]), ncommon = length(common)) }

## prebuild + cache kernels and structure-level diagnostics for a structure
prep_struct <- function(M) {
  KA <- Amat(M); KD <- Dmat(M); KO <- build_ortho_kernels(M)
  list(M = M, n = nrow(M), KA = KA, KD = KD, KO = KO,
       ad_cos = fcos(KA, KD),
       eff_rank = design_diag(list(KO$A, KO$D2, KO$D3, KO$D4))$eff_rank) }

vget <- function(fm, nm) fm$ETA[[nm]]$varU

## fit one BGLR model with warning/error capture + runtime; returns list(varU=named, ok, warn, secs)
fit_capture <- function(y, ETA, pfx) {
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fm <- withCallingHandlers(
    tryCatch(BGLR(y = y, ETA = ETA, nIter = NITER, burnIn = BURN, verbose = FALSE, saveAt = pfx),
             error = function(e) { warns <<- c(warns, paste0("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(fm)) return(list(varU = NULL, ok = FALSE, warn = warns, secs = secs))
  vu <- sapply(names(ETA), function(nm) vget(fm, nm))
  list(varU = vu, ok = all(is.finite(vu)), warn = warns, secs = secs) }

NA_row <- function() data.frame(estimated_VA = NA_real_, estimated_VD_standard = NA_real_, standard_D_share = NA_real_,
  orth_A_fraction = NA_real_, orth_D2_fraction = NA_real_, orth_D3_fraction = NA_real_, orth_D4_fraction = NA_real_,
  model_converged = NA, runtime_seconds = NA_real_, warnings = "")

## ---- one replicate: simulate mixed A+D2 trait at fraction s, fit 3 models ----
run_one <- function(P, s, scode, fcode, rep, nmark) {
  set.seed(100000L * scode + 1000L * fcode + rep)
  M <- P$M; n <- P$n
  gA  <- scale(build_comp(M, P$`causal`, 1L))[, 1]   # standardized additive component
  gD2 <- scale(build_comp(M, P$`causal`, 2L))[, 1]   # standardized digenic  component
  G   <- sqrt(1 - s) * gA + sqrt(s) * gD2
  ve  <- var(G) * (1 - H2) / H2; y <- G + rnorm(n, 0, sqrt(ve))
  # realized true variance fractions (marginal, matching the s definition)
  vA_comp <- var(sqrt(1 - s) * gA); vD_comp <- var(sqrt(s) * gD2)
  realA <- vA_comp / (vA_comp + vD_comp); realD2 <- vD_comp / (vA_comp + vD_comp)
  ad_cor <- suppressWarnings(cor(gA, gD2)); realH2 <- var(G) / var(y)
  pfx <- file.path(tempdir(), sprintf("c_%d_%d_%d_", scode, fcode, rep))
  base <- data.frame(replicate_id = rep, ploidy = PL, population_structure = NA_character_,  # set in run_cell
                     n = n, M = nmark, H2 = H2,
                     target_true_A_fraction = 1 - s, target_true_D2_fraction = s,
                     realized_true_A_fraction = realA, realized_true_D2_fraction = realD2,
                     realized_gA_gD2_cor = ad_cor, realized_H2 = realH2,
                     A_D_Frobenius_cosine = P$ad_cos, effective_rank = P$eff_rank,
                     stringsAsFactors = FALSE)
  mk <- function(model, fields) cbind(base, data.frame(model = model, stringsAsFactors = FALSE), fields)
  # (1) A only
  fA <- fit_capture(y, list(A = list(K = P$KA, model = "RKHS")), paste0(pfx, "a_"))
  rowA <- NA_row(); rowA$estimated_VA <- if (fA$ok) fA$varU[["A"]] else NA_real_
  rowA$model_converged <- fA$ok; rowA$runtime_seconds <- fA$secs; rowA$warnings <- paste(fA$warn, collapse = " | ")
  # (2) standard A + D
  fS <- fit_capture(y, list(A = list(K = P$KA, model = "RKHS"), D = list(K = P$KD, model = "RKHS")), paste0(pfx, "s_"))
  rowS <- NA_row()
  if (fS$ok) { vA <- fS$varU[["A"]]; vD <- fS$varU[["D"]]
    rowS$estimated_VA <- vA; rowS$estimated_VD_standard <- vD; rowS$standard_D_share <- vD / (vA + vD) }
  rowS$model_converged <- fS$ok; rowS$runtime_seconds <- fS$secs; rowS$warnings <- paste(fS$warn, collapse = " | ")
  # (3) orthogonal A + D2 + D3 + D4
  fO <- fit_capture(y, list(A = list(K = P$KO$A, model = "RKHS"), D2 = list(K = P$KO$D2, model = "RKHS"),
                            D3 = list(K = P$KO$D3, model = "RKHS"), D4 = list(K = P$KO$D4, model = "RKHS")), paste0(pfx, "o_"))
  rowO <- NA_row()
  if (fO$ok) { v <- c(fO$varU[["A"]], fO$varU[["D2"]], fO$varU[["D3"]], fO$varU[["D4"]]); fr <- v / sum(v)
    rowO$estimated_VA <- fO$varU[["A"]]
    rowO$orth_A_fraction <- fr[1]; rowO$orth_D2_fraction <- fr[2]; rowO$orth_D3_fraction <- fr[3]; rowO$orth_D4_fraction <- fr[4] }
  rowO$model_converged <- fO$ok; rowO$runtime_seconds <- fO$secs; rowO$warnings <- paste(fO$warn, collapse = " | ")
  rbind(mk("A", rowA), mk("A+D", rowS), mk("A+D2+D3+D4", rowO)) }

## run all reps for one (structure, fraction) cell, parallel over reps
run_cell <- function(PREP, st, s, panel) {
  scode <- match(st, STRUCTURES); fcode <- match(s, FRACTIONS); P <- PREP[[st]]
  P$`causal` <- panel$causal
  resl <- mclapply(1:REPS, function(r) {
    df <- run_one(P, s, scode, fcode, r, panel$nmark); df$population_structure <- st; df }, mc.cores = NC)
  bad <- which(!sapply(resl, is.data.frame))
  if (length(bad)) cat(sprintf("  WARN [%s s=%.2f]: %d/%d reps failed; first error: %s\n",
      st, s, length(bad), length(resl), as.character(resl[[bad[1]]])))
  ok <- resl[setdiff(seq_along(resl), bad)]
  do.call(rbind, ok) }

## ===================== driver =====================
cat(sprintf("calibration-curve sim: PL=%d H2=%.2f REPS=%d NC=%d\n", PL, H2, REPS, NC))
cat("building shared founder panel (3 structures)...\n")
panel <- build_panel()
cat(sprintf("  founder pool n=450; %d markers; %d polymorphic in all structures; 200 causal\n",
            panel$nmark, panel$ncommon))
cat("building + caching kernels per structure...\n")
PREP <- lapply(STRUCTURES, function(st) prep_struct(panel$POPS[[st]])); names(PREP) <- STRUCTURES
for (st in STRUCTURES) cat(sprintf("  [%-9s] n=%d  A-D Frobenius cos=%.3f  4-kernel eff_rank=%.3f\n",
                                   st, PREP[[st]]$n, PREP[[st]]$ad_cos, PREP[[st]]$eff_rank))

cells <- expand.grid(structure = STRUCTURES, s = FRACTIONS, stringsAsFactors = FALSE)
AID <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (Sys.getenv("DRYRUN") == "1") { cells <- cells[cells$structure == "fullsib" & cells$s == 0.50, , drop = FALSE]; REPS <- 2L }

if (nzchar(AID)) {                                   # ---- ARRAY MODE: one cell per task ----
  ci <- as.integer(AID)
  st <- cells$structure[ci]; s <- cells$s[ci]
  cat(sprintf("ARRAY cell %d/%d: structure=%s s=%.2f, REPS=%d NC=%d\n", ci, nrow(cells), st, s, REPS, NC))
  res <- run_cell(PREP, st, s, panel)
  f <- file.path(OUT, sprintf("calib_cell_%02d.csv", ci))
  write.csv(res, f, row.names = FALSE)
  cat(sprintf("ARRAY cell %d wrote %d rows -> %s\n", ci, nrow(res), f))
} else {                                             # ---- LOCAL/LEGACY: all cells ----
  cat(sprintf("running %d cells x %d reps x 3 models, NC=%d\n", nrow(cells), REPS, NC))
  resl <- lapply(seq_len(nrow(cells)), function(ci) {
    st <- cells$structure[ci]; s <- cells$s[ci]
    cat(sprintf("  cell %d/%d structure=%s s=%.2f\n", ci, nrow(cells), st, s))
    run_cell(PREP, st, s, panel) })
  res <- do.call(rbind, resl)
  out <- if (Sys.getenv("DRYRUN") == "1") file.path(OUT, "calibration_smoke.csv") else file.path(OUT, "calibration_primary.csv")
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("wrote %d rows -> %s\n", nrow(res), out))
  ## quick summary: standard D share & orthogonal D2 fraction vs true s
  ad <- res[res$model == "A+D", ]; or <- res[res$model == "A+D2+D3+D4", ]
  cat("mean standard_D_share by structure x target s:\n"); print(aggregate(standard_D_share ~ population_structure + target_true_D2_fraction, ad, mean))
  cat("mean orth_D2_fraction by structure x target s:\n"); print(aggregate(orth_D2_fraction ~ population_structure + target_true_D2_fraction, or, mean))
}
cat("CALIBRATION SIM DONE\n")

#!/usr/bin/env Rscript
# =====================================================================
# SECONDARY calibration stress grid (plan ��9). Does the calibration of the
# standard A+D dominance share (and the orthogonal D2 fraction) against the TRUE
# digenic fraction s improve with larger n, more markers, or different LD?
#
# Reuses the validated, memory-safe machinery of the EXPANDED regime map
# (block-wise kernel accumulation -> O(n^2) peak; truncated-EVD low-rank RKHS fit
# -> no M<<n rank-deficiency artifact, MCMC O(n*rank); LD set by chromosome genetic
# length via runMacs2; realized adjacent-marker r2 measured). We source that file
# with REGIME_FUNCS_ONLY=1 to import make_pops(), measure_LD(), apply_dr_vec(),
# ortho(), design_diag(), fcos() and the engineering constants, then add:
#   - build_blockwise_mix(): accumulates standard A,D + (optionally) the four
#     orthogonal kernels OA,OD2,OD3,OD4 AND two standardized polygenic components
#     g_A (degree 1) and g_D2 (degree 2), all block-wise over markers;
#   - a MIXED trait  G = sqrt(1-s)*g_A* + sqrt(s)*g_D2*  at controlled true s;
#   - standard A+D fit -> standard_D_share, and (gated) orthogonal fit -> orth D2.
#
# Grid (full-sib; the plan's main stress structure):
#   N  in {400, 1200, 4800, 19200}
#   M  in {5000, 50000, 200000}
#   LD in {independent (low LD, r2~0.11), high (r2~0.22)}
#   s  in {0, 0.10, 0.25, 0.50}     dr = 0
#   reps 50 (N<=4800) / 30 (N=19200)
#
# Cost control: the orthogonal fit (4 extra n-sized eigendecompositions + a
# 4-kernel build over ALL markers) is gated to N <= CALIB_ORTH_MAXN (default 4800);
# at larger N orth_D2_fraction is NA and eff_rank uses a capped marker subsample
# (MCAP_DIAG), exactly as the expanded regime map does. The cheap standard share
# is computed at every cell. Set CALIB_ORTH_MAXN=19200 to force orth everywhere.
#
# Run modes (mirrors regime_map_expanded.R dispatcher):
#   grain=cell    -> task = one (N,M,LD,s) cell, mclapply over reps, BLAS pinned 1
#   grain=cellrep -> task = one cell-rep, threaded BLAS (for the large-N eigendecomp)
#   CHUNK groups consecutive jobs/task; SLURM_ARRAY_TASK_ID slices the work-list.
#   DRYRUN=1 -> tiny smoke (N=400, M=5000, s in {0,0.5}, 2 reps).
# Output: out/calib_sec_<tier>_<NNNN>.csv ; concatenate afterwards.
# =====================================================================
Sys.setenv(REGIME_FUNCS_ONLY = "1")
source("experiments/simulations/regime_map/regime_map_expanded.R")   # imports funcs + constants (PL, COLCAP, KCAP, MCAP_DIAG, NE, NITER, BURN, NC, make_pops, ...)
suppressMessages({library(AlphaSimR); library(BGLR); library(parallel)})

H2 <- 0.6                                                 # match the primary calibration + plan (expanded used 0.5)
OUTS <- Sys.getenv("CALIB_OUTDIR", "experiments/simulations/calibration/out"); dir.create(OUTS, showWarnings = FALSE, recursive = TRUE)
ORTH_MAXN <- as.integer(Sys.getenv("CALIB_ORTH_MAXN", "4800"))

## ---- grid axes (own; independent of the expanded NS/MS) --------------
NS2  <- c(400L, 1200L, 4800L, 19200L)
MS2  <- c(5000L, 50000L, 200000L, 1000L, 20000L, 100000L)  # APPENDED in run order so existing M keep their which(MS2==M) seed index -> already-computed data stays reproducible (1000 dropped from the figure; 20000,100000 added for the 5-col set {5k,20k,50k,100k,200k})
LDS2 <- c("independent", "high")          # low-LD vs high-LD extremes (genlen_for from sourced file)
SS2  <- c(0, 0.10, 0.25, 0.50)            # target true digenic fraction s
reps_for2 <- function(N) if (N <= 4800L) as.integer(Sys.getenv("REPS_SMALL", "50")) else as.integer(Sys.getenv("REPS_BIG", "30"))
seed_of2 <- function(N, M, LD, s, rep)
  as.integer((5e4 + which(NS2 == N) * 9001L + which(MS2 == M) * 811L +
              which(LDS2 == LD) * 101L + round(s * 100) * 13L + rep) %% 2e9L)

## ---- block-wise mixed-trait kernels ---------------------------------
## Always builds standard A (VanRaden) and D (digenic), plus standardized
## polygenic g_A (degree-1) and g_D2 (degree-2) over ALL markers. When want_orth,
## also accumulates the four orthogonal kernels OA,OD2,OD3,OD4 over ALL markers
## (used for the orthogonal fit and exact eff_rank); otherwise it accumulates a
## capped Od from the first MCAP_DIAG markers for eff_rank only (cheap at high M).
build_blockwise_mix <- function(pop, SP, nChr, k, dr, want_orth) {
  n <- nInd(pop)
  SA <- matrix(0, n, n); dd <- 0; SD <- matrix(0, n, n)
  SO <- replicate(4, matrix(0, n, n), simplify = FALSE)   # orth kernels (full if want_orth, else capped)
  gA <- numeric(n); gD2 <- numeric(n); vcount <- 0L; od_used <- 0L
  for (cc in 1:nChr) {
    G <- pullSnpGeno(pop, snpChip = 1, chr = cc, simParam = SP); storage.mode(G) <- "double"
    if (dr > 0) G <- apply_dr_vec(G, dr)
    nm <- ncol(G)
    for (s0 in seq(1L, nm, by = COLCAP)) {
      idx <- s0:min(s0 + COLCAP - 1L, nm)
      Gb <- G[, idx, drop = FALSE]
      pk <- colMeans(round(Gb)) / PL
      keep <- which(pk > 0 & pk < 1); if (!length(keep)) next
      Gb <- Gb[, keep, drop = FALSE]; pk <- pk[keep]; nk <- length(keep); vcount <- vcount + nk
      # standard additive + digenic kernels (all markers)
      Z <- sweep(Gb, 2, PL * pk); SA <- SA + tcrossprod(Z); dd <- dd + sum(PL * pk * (1 - pk))
      Hh <- (Gb * (PL - Gb)) / choose(PL, 2); Hh <- sweep(Hh, 2, colMeans(Hh)); SD <- SD + tcrossprod(Hh)
      # orthogonal degree-1..4 contrasts (grouped by rounded allele freq)
      gint <- round(Gb) + 1L
      pkey <- pmin(pmax(round(pk, 3), 0.001), 0.999)
      X1 <- matrix(0, n, nk); X2 <- matrix(0, n, nk)
      X3 <- if (want_orth) matrix(0, n, nk) else NULL; X4 <- if (want_orth) matrix(0, n, nk) else NULL
      for (key in unique(pkey)) { gcc <- which(pkey == key); Q <- ortho(key, PL)
        X1[, gcc] <- matrix(Q[, 2][gint[, gcc, drop = FALSE]], n)
        X2[, gcc] <- matrix(Q[, 3][gint[, gcc, drop = FALSE]], n)
        if (want_orth) { X3[, gcc] <- matrix(Q[, 4][gint[, gcc, drop = FALSE]], n)
                         X4[, gcc] <- matrix(Q[, 5][gint[, gcc, drop = FALSE]], n) } }
      gA  <- gA  + X1 %*% rnorm(nk)
      gD2 <- gD2 + X2 %*% rnorm(nk)
      if (want_orth) {
        SO[[1]] <- SO[[1]] + tcrossprod(X1); SO[[2]] <- SO[[2]] + tcrossprod(X2)
        SO[[3]] <- SO[[3]] + tcrossprod(X3); SO[[4]] <- SO[[4]] + tcrossprod(X4)
      } else if (od_used < MCAP_DIAG) {              # capped eff_rank kernels only
        take <- min(nk, MCAP_DIAG - od_used); sc <- 1:take
        Xc <- list(X1[, sc, drop = FALSE], X2[, sc, drop = FALSE], matrix(0, n, take), matrix(0, n, take))
        pks <- pkey[sc]; gis <- gint[, sc, drop = FALSE]
        for (key in unique(pks)) { gcc <- which(pks == key); Q <- ortho(key, PL)
          for (d in 3:4) Xc[[d]][, gcc] <- matrix(Q[, d + 1][gis[, gcc, drop = FALSE]], n) }
        for (d in 1:4) SO[[d]] <- SO[[d]] + tcrossprod(Xc[[d]])
        od_used <- od_used + take; rm(Xc)
      }
      rm(Gb, Z, Hh, X1, X2, X3, X4, gint);
    }
    rm(G); gc(FALSE)
  }
  nrm <- function(S) S / mean(diag(S))
  KA <- SA / ifelse(dd > 0, dd, 1); KA <- KA / mean(diag(KA))
  list(A = KA, D = nrm(SD), Od = lapply(SO, nrm),
       gA = scale(gA)[, 1], gD2 = scale(gD2)[, 1], valid = vcount)
}

## low-rank truncated-EVD RKHS term (drop numerical-zero directions)
lr_term <- function(K) { e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-8 * max(e$values)
  list(V = e$vectors[, keep, drop = FALSE], d = e$values[keep], model = "RKHS") }

## ---- one cell-rep ---------------------------------------------------
fail_row <- function(N, M, LD, s, rep) data.frame(N = N, nMarker = M, LD = LD, adjLD_r2 = NA_real_, dr = 0,
  structure = "fullsib", rep = rep, target_D2 = s, realized_D2 = NA_real_,
  standard_D_share = NA_real_, orth_D2_fraction = NA_real_,
  calibration_bias_standard = NA_real_, calibration_bias_orth = NA_real_,
  A_D_Frobenius_cosine = NA_real_, effective_rank = NA_real_, valid_markers = NA_integer_,
  H2 = H2, failed = 1L)

run_cellrep_calib <- function(N, M, LD, s, rep) {
  base_seed <- seed_of2(N, M, LD, s, rep); set.seed(base_seed)
  mk <- tryCatch(make_pops(N, M, LD, fullsib_only = TRUE), error = function(e) { cat("make_pops ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(mk)) return(fail_row(N, M, LD, s, rep))
  SP <- mk$SP
  adjLD <- tryCatch(measure_LD(mk$founders, SP, mk$nChr, 0), error = function(e) NA_real_)
  pop <- mk$pops[["fullsib"]]; if (is.null(pop) || nInd(pop) < 20L) return(fail_row(N, M, LD, s, rep))
  want_orth <- (Sys.getenv("CALIB_ORTH", "1") == "1") && (N <= ORTH_MAXN)
  set.seed(base_seed * 7L + 1L)
  K <- build_blockwise_mix(pop, SP, mk$nChr, mk$k, 0, want_orth = want_orth)
  G  <- sqrt(1 - s) * K$gA + sqrt(s) * K$gD2
  ve <- var(G) * (1 - H2) / H2; y <- G + rnorm(length(G), 0, sqrt(ve))
  realD2 <- { vA <- var(sqrt(1 - s) * K$gA); vD <- var(sqrt(s) * K$gD2); if (vA + vD > 0) vD / (vA + vD) else 0 }
  pfx <- file.path(tempdir(), sprintf("cs_%d_%d_", base_seed, rep))
  # standard A+D (low-rank)
  fs <- BGLR(y = y, ETA = list(A = lr_term(K$A), D = lr_term(K$D)),
             nIter = NITER, burnIn = BURN, verbose = FALSE, saveAt = paste0(pfx, "s_"))
  std_share <- fs$ETA$D$varU / (fs$ETA$A$varU + fs$ETA$D$varU)
  adcos <- round(fcos(K$A, K$D), 3); effr <- round(design_diag(K$Od)$eff_rank, 3)
  orthD2 <- NA_real_
  if (want_orth) {
    fo <- BGLR(y = y, ETA = list(A = lr_term(K$Od[[1]]), D2 = lr_term(K$Od[[2]]),
                                 D3 = lr_term(K$Od[[3]]), D4 = lr_term(K$Od[[4]])),
               nIter = NITER, burnIn = BURN, verbose = FALSE, saveAt = paste0(pfx, "o_"))
    v <- c(fo$ETA$A$varU, fo$ETA$D2$varU, fo$ETA$D3$varU, fo$ETA$D4$varU); orthD2 <- v[2] / sum(v)
  }
  vmk <- K$valid; rm(K, fs); gc(FALSE)
  data.frame(N = N, nMarker = M, LD = LD, adjLD_r2 = round(adjLD, 4), dr = 0, structure = "fullsib", rep = rep,
    target_D2 = s, realized_D2 = round(realD2, 4),
    standard_D_share = round(std_share, 4), orth_D2_fraction = round(orthD2, 4),
    calibration_bias_standard = round(std_share - s, 4), calibration_bias_orth = round(orthD2 - s, 4),
    A_D_Frobenius_cosine = adcos, effective_rank = effr, valid_markers = vmk,
    H2 = H2, failed = 0L)
}

## ---- dispatcher (mirrors regime_map_expanded.R) ---------------------
tier  <- Sys.getenv("CALIB_TIER", "1")
grain <- Sys.getenv("REGIME_GRAIN", "cell")
Nsel <- NS2; Msel <- MS2; LDsel <- LDS2; Ssel <- SS2
e <- Sys.getenv("CALIB_NSET", ""); if (nzchar(e)) Nsel <- as.integer(strsplit(e, ",")[[1]])
e <- Sys.getenv("CALIB_MSET", ""); if (nzchar(e)) Msel <- as.integer(strsplit(e, ",")[[1]])
e <- Sys.getenv("CALIB_LDSET", ""); if (nzchar(e)) LDsel <- strsplit(e, ",")[[1]]
e <- Sys.getenv("CALIB_SSET", ""); if (nzchar(e)) Ssel <- as.numeric(strsplit(e, ",")[[1]])

cells <- expand.grid(N = Nsel, nMarker = Msel, LD = LDsel, s = Ssel, stringsAsFactors = FALSE)
if (grain == "cellrep") {
  jobs <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) data.frame(cells[i, ], rep = seq_len(reps_for2(cells$N[i])))))
} else jobs <- cbind(cells, rep = NA_integer_)

if (Sys.getenv("DRYRUN") == "1") {
  jobs <- expand.grid(N = 400L, nMarker = 5000L, LD = "independent", s = c(0, 0.5), rep = NA_integer_, stringsAsFactors = FALSE)
  grain <- "cell"; assign("reps_for2", function(N) 2L, envir = environment())
}

## ---- TOPUP mode: run an explicit (N,M,LD,s,rep) work-list -----------
## When TOPUP_JOBS points to a CSV (column s or target_D2 for the true digenic
## fraction), ignore the generated grid and run exactly those cell-reps (grain
## forced to cellrep). seed_of2() is index-based on the fixed NS2/MS2/LDS2/SS2
## vectors, so listed reps reuse the SAME deterministic seeds as the original
## run. SLURM_ARRAY_TASK_ID + CHUNK slice the work-list as usual.
topf <- Sys.getenv("TOPUP_JOBS", "")
if (nzchar(topf)) {
  tj <- read.csv(topf, stringsAsFactors = FALSE)
  sv <- if ("s" %in% names(tj)) tj$s else tj$target_D2
  jobs <- data.frame(N = as.integer(tj$N), nMarker = as.integer(tj$nMarker),
                     LD = as.character(tj$LD), s = as.numeric(sv),
                     rep = as.integer(tj$rep), stringsAsFactors = FALSE)
  grain <- "cellrep"
  cat(sprintf("CALIB TOPUP mode: %d cell-reps from %s\n", nrow(jobs), topf))
}

AID   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
CHUNK <- as.integer(Sys.getenv("CHUNK", as.character(nrow(jobs))))
sel   <- ((AID - 1L) * CHUNK + 1L):min(AID * CHUNK, nrow(jobs)); sel <- sel[sel >= 1L & sel <= nrow(jobs)]
cat(sprintf("CALIB-SECONDARY tier=%s grain=%s H2=%.2f ORTH_MAXN=%d: %d jobs, task %d rows %s (NC=%d)\n",
            tier, grain, H2, ORTH_MAXN, nrow(jobs), AID, paste(range(sel), collapse = "-"), NC))

emit <- list()
for (i in sel) {
  jr <- jobs[i, ]
  if (grain == "cellrep") {
    emit[[length(emit) + 1L]] <- tryCatch(run_cellrep_calib(jr$N, jr$nMarker, jr$LD, jr$s, jr$rep),
                                          error = function(e) { cat("ERR cellrep:", conditionMessage(e), "\n"); NULL })
  } else {
    R <- reps_for2(jr$N)
    rr <- mclapply(1:R, function(rp) tryCatch(run_cellrep_calib(jr$N, jr$nMarker, jr$LD, jr$s, rp),
                   error = function(e) NULL), mc.cores = NC)
    emit[[length(emit) + 1L]] <- do.call(rbind, Filter(Negate(is.null), rr))
  }
}
res <- do.call(rbind, Filter(Negate(is.null), emit))
outf <- file.path(OUTS, sprintf("calib_sec_t%s_%04d.csv", tier, AID))
write.csv(res, outf, row.names = FALSE)
cat(sprintf("wrote %d rows -> %s\n", if (is.null(res)) 0L else nrow(res), outf))
if (!is.null(res)) print(aggregate(standard_D_share ~ N + nMarker + target_D2, res, function(z) round(mean(z, na.rm = TRUE), 3)))
cat("CALIB-SECONDARY DONE\n")

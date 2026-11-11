#!/usr/bin/env Rscript
# =====================================================================
# EXPANDED known-truth regime-map (idea:010) — comprehensive stress test.
# Known-truth = PURELY ADDITIVE (zero dominance by construction). Measures the
# SPURIOUS standard digenic Vd-share V_D/(V_A+V_D), plus the additive-digenic
# Frobenius cosine and four-kernel effective rank, across a LARGE factorial:
#     N      in {150,300,600,1200,2400,4800,9600,19200}
#     nMarker in {1e3,5e3,1e4,2e4,5e4,1e5,2e5}
#     LD     in {independent,low,medium,high}   (realized adjacent-r2 reported)
#     dr     in {0, 0.1}                         (mean-preserving double reduction)
#     reps   50 (N<=4800) / 30 (N>=9600)
#     structures unrelated/half-sib/full-sib (N<=4800); full-sib only (N>=9600)
#
# LD axis (answers "how is LD set?"): a fixed, always-feasible marker density per
# chromosome (KPC markers/chr) with LD set by CHROMOSOME GENETIC LENGTH via
# runMacs2 genLen (shorter chromosome -> tighter linkage). This decouples the LD
# axis from the coalescent's segregating-site ceiling, so high LD works even at
# M=2e5. Realized adjacent-marker r2 is MEASURED on an unrelated founder sample
# and reported (adjLD_r2); labels are descriptive only.
#
# ENGINEERING for the large corner (n^3 eigendecomposition + n*M kernel build):
#   - BLOCK-WISE kernel accumulation: stream marker columns per chromosome and
#     accumulate ZZ' so peak memory is O(n^2) (~3 GB/kernel at n=19200), NOT
#     O(n*M) (~30 GB/contrast at the corner).
#   - GROUPED orthogonal contrasts: cache ortho(p) by rounded allele freq and
#     vectorize the degree-1..4 contrast fill (fast at M=2e5).
#   - eff_rank uses CAPPED (MCAP_DIAG) markers; AD cosine is pure linear algebra;
#     NO orthogonal BGLR fit (only the standard A+D fit needs n x n eigendecomp).
#   - HYBRID parallelism (env REGIME_GRAIN): 'cell' -> mclapply reps, BLAS=1
#     (small N); 'cellrep' -> one rep/task, threaded BLAS (large N). CHUNK groups
#     consecutive jobs per array task; SLURM_ARRAY_TASK_ID slices the work-list.
# Uses manual NOIA orthonormal contrasts (lib_ortho.R), NOT addTraitAD.
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho.R")          # ortho(), design_diag(), fcos(), BGLR
suppressMessages({library(AlphaSimR); library(parallel)})
OUT <- Sys.getenv("REGIME_OUTDIR", "experiments/simulations/regime_map/out/expanded"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PL <- 4L; H2 <- 0.5
# BGLR RKHS MCMC is O(n^2)/iteration (eigenvector rotation), so iterations dominate at large n.
# A 2-component variance RATIO converges well under ~1500 iters; default tuned for speed at scale.
NITER <- as.integer(Sys.getenv("NITER", "2000")); BURN <- as.integer(Sys.getenv("BURN", "500"))
NC <- as.integer(Sys.getenv("NC", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
COLCAP    <- as.integer(Sys.getenv("COLCAP", "4000"))    # max marker columns held at once
KCAP      <- as.integer(Sys.getenv("KCAP", "4000"))      # markers per chromosome (keeps coalescent feasible)
MCAP_DIAG <- as.integer(Sys.getenv("MCAP_DIAG", "8000")) # markers used for eff_rank kernels
NE        <- as.integer(Sys.getenv("NE", "1000"))        # runMacs2 Ne (sets segregating-site supply)
STRUCT_BIG_N <- 9600L                                    # at/above this N: full-sib only (cost control)
STRUCTURES <- c("unrelated", "halfsib", "fullsib")

# ---- factorial axes -------------------------------------------------
NS  <- c(150L, 300L, 600L, 1200L, 2400L, 4800L, 9600L, 19200L)
MS  <- c(1000L, 5000L, 10000L, 20000L, 50000L, 100000L, 200000L)
LDS <- c("independent", "low", "medium", "high")
DRS <- c(0, 0.1)
reps_for <- function(N) {
  if (N <= 4800L) as.integer(Sys.getenv("REPS_SMALL", "50")) else as.integer(Sys.getenv("REPS_BIG", "30"))
}
# LD axis = chromosome genetic length (Morgans/chr). SMALL genLen -> tight linkage
# -> high LD. Values kept <=1.5 so the coalescent stays fast (large genLen blows up
# the ancestral recombination graph). Markers-per-chr is fixed at KCAP (feasible
# site count); LD is set purely by genLen. Realized adjacent r2 is measured/reported.
genlen_for <- function(LD) switch(LD, independent = 1.5, low = 0.4, medium = 0.1, high = 0.025)

# Mean-preserving double-reduction proxy (excess homozygosity), VECTORIZED per block:
# a heterozygote dosage g is, with prob dr, sent to {0,PL} with P(->PL)=g/PL so E[new]=g.
apply_dr_vec <- function(Gb, dr) {
  if (dr <= 0) return(Gb)
  het <- Gb > 0 & Gb < PL
  hit <- het & (matrix(runif(length(Gb)), nrow(Gb)) < dr)
  if (any(hit)) { g <- Gb[hit]; Gb[hit] <- ifelse(runif(length(g)) < g / PL, PL, 0) }
  Gb
}

# ---- founder pool per cell-rep -> structures ------------------------
# Coalescent marker-generation cost ~ M x (founder haplotypes), so for the
# full-sib-only large-N path we shrink the founder pool to just the parents
# (~2 x nFam) instead of ~1.125 N -- a ~4-5x speed-up at the corner. The
# all-structure path needs >= N unrelated founders, so keeps the larger pool.
make_pops <- function(N, M, LD, fullsib_only) {
  nChr <- as.integer(max(1L, ceiling(M / KCAP)))
  k    <- as.integer(ceiling(M / nChr))                 # ~KCAP markers/chr (uniform)
  bp   <- max(5e6, k * 2500)                            # physical length to supply k segSites
  nFam <- max(2L, ceiling(N / 8))
  Nf   <- if (fullsib_only) max(50L, 2L * nFam + 10L) else max(50L, 9L * nFam)
  fd <- runMacs2(nInd = Nf, nChr = nChr, segSites = k + 30L, Ne = NE, bp = bp,
                 genLen = genlen_for(LD), mutRate = 2.5e-8, ploidy = PL)
  SP <- SimParam$new(fd); SP$addSnpChip(nSnpPerChr = k)
  founders <- newPop(fd, simParam = SP); nf <- nInd(founders)
  if (fullsib_only) {
    fs <- randCross(founders, nCrosses = nFam, nProgeny = 8, simParam = SP); if (nInd(fs) > N) fs <- fs[1:N]
    return(list(SP = SP, founders = founders, nChr = nChr, k = k, pops = list(fullsib = fs)))
  }
  un <- founders[1:min(N, nf)]
  fs <- randCross(founders[1:min(nf, max(2L, ceiling(N / 4)))], nCrosses = nFam, nProgeny = 8, simParam = SP)
  if (nInd(fs) > N) fs <- fs[1:N]
  hs <- NULL
  if (nf >= nFam + nFam * 8) {
    sires <- 1:nFam; damsM <- matrix(nFam + 1:(nFam * 8), nrow = nFam, byrow = TRUE)
    plan <- do.call(rbind, lapply(1:nFam, function(f) cbind(damsM[f, ], rep(sires[f], 8))))
    hs <- makeCross(founders, crossPlan = plan, nProgeny = 1L, simParam = SP)
    if (nInd(hs) > N) hs <- hs[1:N]
  }
  list(SP = SP, founders = founders, nChr = nChr, k = k,
       pops = list(unrelated = un, halfsib = hs, fullsib = fs))
}

# Realized population LD: mean adjacent-marker r2 on an UNRELATED founder sample
# (consistent across N; family structure would otherwise inflate r2).
measure_LD <- function(founders, SP, nChr, dr) {
  ids <- 1:min(400L, nInd(founders)); r2 <- numeric(0)
  for (cc in 1:nChr) {
    G <- pullSnpGeno(founders[ids], snpChip = 1, chr = cc, simParam = SP); storage.mode(G) <- "double"
    if (dr > 0) G <- apply_dr_vec(G, dr)
    if (ncol(G) < 2) next
    jj <- seq_len(ncol(G) - 1L); if (length(jj) > 1500L) jj <- sort(sample(jj, 1500L))
    rr <- vapply(jj, function(j) { s <- suppressWarnings(cor(G[, j], G[, j + 1L])); if (is.finite(s)) s^2 else NA_real_ }, numeric(1))
    r2 <- c(r2, rr[is.finite(rr)])
    if (length(r2) > 4000L) break
  }
  if (length(r2)) mean(r2) else NA_real_
}

# ---- block-wise kernels + purely-additive gv ------------------------
# Streams markers per chromosome. A (VanRaden) and D (digenic) use ALL markers;
# the additive truth gv is accumulated from the degree-1 contrast over ALL markers.
# When want_diag, the 4 orthogonal kernels for eff_rank are accumulated from only
# the first MCAP_DIAG markers (a representative subsample -> bounded cost at high M).
build_blockwise <- function(pop, SP, nChr, k, dr, want_diag = FALSE, want_truth = FALSE) {
  n <- nInd(pop)
  SA <- matrix(0, n, n); dd <- 0; SD <- matrix(0, n, n)
  SO <- if (want_diag) replicate(4, matrix(0, n, n), simplify = FALSE) else NULL
  gv <- numeric(n); vcount <- 0L; od_used <- 0L
  gA <- numeric(n); gD2 <- numeric(n); gD3 <- numeric(n); gD4 <- numeric(n)   # want_truth components
  for (cc in 1:nChr) {
    G <- pullSnpGeno(pop, snpChip = 1, chr = cc, simParam = SP); storage.mode(G) <- "double"
    if (dr > 0) G <- apply_dr_vec(G, dr)
    nm <- ncol(G)
    for (s0 in seq(1L, nm, by = COLCAP)) {
      idx <- s0:min(s0 + COLCAP - 1L, nm)
      Gb <- G[, idx, drop = FALSE]
      pk <- colMeans(round(Gb)) / PL
      keep <- which(pk > 0 & pk < 1)
      if (!length(keep)) next
      Gb <- Gb[, keep, drop = FALSE]; pk <- pk[keep]; nk <- length(keep); vcount <- vcount + nk
      # standard additive + digenic over ALL markers
      Z <- sweep(Gb, 2, PL * pk); SA <- SA + tcrossprod(Z); dd <- dd + sum(PL * pk * (1 - pk))
      Hh <- (Gb * (PL - Gb)) / choose(PL, 2); Hh <- sweep(Hh, 2, colMeans(Hh)); SD <- SD + tcrossprod(Hh)
      # degree-1 contrast for the additive truth gv (grouped by rounded allele freq)
      gint <- round(Gb) + 1L
      X1 <- matrix(0, n, nk)
      pkey <- pmin(pmax(round(pk, 3), 0.001), 0.999)   # clamp off 0/1 so ortho() Gram-Schmidt stays finite
      ug <- unique(pkey)
      for (key in ug) { gcc <- which(pkey == key); X1[, gcc] <- matrix(ortho(key, PL)[, 2][gint[, gcc, drop = FALSE]], n) }
      if (want_truth) {
        # locus-level truth: accumulate degree-1..4 genetic-value components over ALL
        # markers (each marker draws an i.i.d. effect per degree). Reuse X1 as degree 1.
        gA <- gA + as.numeric(X1 %*% rnorm(nk))   # as.numeric: keep components plain vectors (not n x 1 matrices)
        for (dg in 2:4) {
          Xdg <- matrix(0, n, nk)
          for (key in ug) { gcc <- which(pkey == key); Xdg[, gcc] <- matrix(ortho(key, PL)[, dg + 1][gint[, gcc, drop = FALSE]], n) }
          if (dg == 2L) gD2 <- gD2 + as.numeric(Xdg %*% rnorm(nk))
          else if (dg == 3L) gD3 <- gD3 + as.numeric(Xdg %*% rnorm(nk))
          else gD4 <- gD4 + as.numeric(Xdg %*% rnorm(nk))
          rm(Xdg)
        }
      } else gv <- gv + X1 %*% rnorm(nk)
      # eff_rank kernels (degrees 1..4) from the first MCAP_DIAG markers only
      if (want_diag && od_used < MCAP_DIAG) {
        take <- min(nk, MCAP_DIAG - od_used); sc <- 1:take
        Xd <- list(X1[, sc, drop = FALSE], matrix(0, n, take), matrix(0, n, take), matrix(0, n, take))
        pks <- pkey[sc]; gis <- gint[, sc, drop = FALSE]
        for (key in unique(pks)) { gcc <- which(pks == key); Q <- ortho(key, PL)
          for (d in 2:4) Xd[[d]][, gcc] <- matrix(Q[, d + 1][gis[, gcc, drop = FALSE]], n) }
        for (d in 1:4) SO[[d]] <- SO[[d]] + tcrossprod(Xd[[d]])
        od_used <- od_used + take; rm(Xd)
      }
      rm(Gb, Z, Hh, X1)
    }
    rm(G); gc(FALSE)
  }
  # NO diagonal jitter: the standard A/D kernels are handed to BGLR as a TRUNCATED
  # eigendecomposition (low-rank V/d) in run_cellrep, so the old +diag(1e-6) full-rank
  # padding is both unnecessary and harmful at M<<n -- it manufactures ~(n-M) junk
  # eigendirections (eigenvalue ~1e-6) spanning the whole space, which the D kernel then
  # uses to absorb the additive signal -> spuriously inflated V_D. Dropping the jitter
  # lets the EVD threshold cleanly separate the ~M genuine marker directions from the
  # numerical-zero ones. (See RESUME_HANDOFF.md / diag_fast.R for the reproduction.)
  nrm <- function(S) S / mean(diag(S))
  KA <- SA / ifelse(dd > 0, dd, 1); KA <- KA / mean(diag(KA))
  # In want_truth mode the legacy additive-only gv is unused (stays 0) -> don't scale (sd=0 -> NaN).
  out <- list(A = KA, D = nrm(SD), gv = if (want_truth) gv else scale(gv)[, 1], valid = vcount)
  if (want_diag) out$Od <- lapply(SO, nrm)
  if (want_truth) { out$gA <- gA; out$gD2 <- gD2; out$gD3 <- gD3; out$gD4 <- gD4 }
  out
}

# ---- one (N,M,LD,dr,rep) cell-rep -> structure rows -----------------
seed_of <- function(N, M, LD, dr, rep)
  as.integer(1e3 + (which(NS == N) * 1303L + which(MS == M) * 131L +
             which(LDS == LD) * 17L + (dr > 0) * 7L + rep) %% 2e9L)

run_cellrep <- function(N, M, LD, dr, rep) {
  base_seed <- seed_of(N, M, LD, dr, rep); set.seed(base_seed)
  fullsib_only <- N >= STRUCT_BIG_N
  mk <- tryCatch(make_pops(N, M, LD, fullsib_only), error = function(e) { cat("make_pops ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(mk)) return(data.frame(N = N, nMarker = M, LD = LD, adjLD_r2 = NA_real_, dr = dr,
                                     structure = "ALL", rep = rep, spurious_Vd = NA_real_,
                                     nonadd_misattrib = NA_real_, eff_rank = NA_real_, AD_cos = NA_real_,
                                     valid_markers = NA_integer_, failed = 1L))
  SP <- mk$SP; struct_set <- if (N >= STRUCT_BIG_N) "fullsib" else STRUCTURES
  adjLD_cell <- tryCatch(measure_LD(mk$founders, SP, mk$nChr, dr), error = function(e) NA_real_)
  rows <- list()
  for (st in struct_set) {
    pop <- mk$pops[[st]]; if (is.null(pop) || nInd(pop) < 20L) next
    set.seed(base_seed * 7L + match(st, STRUCTURES))            # reproducible per-structure truth
    is_fs <- st == "fullsib"
    K <- build_blockwise(pop, SP, mk$nChr, mk$k, dr, want_diag = is_fs)
    gv <- K$gv; ve <- var(gv) * (1 - H2) / H2; y <- gv + rnorm(length(gv), 0, sqrt(ve))
    pfx <- file.path(tempdir(), sprintf("re_%d_%s_%d_", base_seed, st, rep))
    # Fit the standard A+D RKHS model from a TRUNCATED eigendecomposition (low-rank V/d):
    # keep only directions with eigenvalue > 1e-8*max -- the ~M genuine marker dimensions --
    # and drop the numerical-zero directions. This removes the M<<n rank-deficiency artifact
    # (the D kernel absorbing the additive signal through junk eigendirections) and makes the
    # RKHS MCMC O(n*rank) rather than O(n^2). For M>=n (full-rank kernels) ~all n directions
    # are kept, so the fit is identical to the previous full-K fit.
    eA <- eigen(K$A, symmetric = TRUE); eD <- eigen(K$D, symmetric = TRUE)
    kA <- eA$values > 1e-8 * max(eA$values); kD <- eD$values > 1e-8 * max(eD$values)
    fs <- BGLR(y = y, ETA = list(
                 A = list(V = eA$vectors[, kA, drop = FALSE], d = eA$values[kA], model = "RKHS"),
                 D = list(V = eD$vectors[, kD, drop = FALSE], d = eD$values[kD], model = "RKHS")),
               nIter = NITER, burnIn = BURN, verbose = FALSE, saveAt = paste0(pfx, "s_"))
    vshare <- fs$ETA$D$varU / (fs$ETA$A$varU + fs$ETA$D$varU)
    effr <- NA_real_; adcos <- NA_real_
    if (is_fs) { adcos <- round(fcos(K$A, K$D), 3)
      effr <- round(design_diag(K$Od)$eff_rank, 3) }
    rows[[st]] <- data.frame(N = N, nMarker = M, LD = LD, adjLD_r2 = adjLD_cell, dr = dr, structure = st,
      rep = rep, spurious_Vd = round(vshare, 3), nonadd_misattrib = NA_real_,
      eff_rank = effr, AD_cos = adcos, valid_markers = K$valid, failed = 0L)
    rm(K, fs); gc(FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---- build the tier's job list and slice it by SLURM_ARRAY_TASK_ID ---
if (Sys.getenv("REGIME_FUNCS_ONLY") != "1") {
tier  <- Sys.getenv("REGIME_TIER", "1")
grain <- Sys.getenv("REGIME_GRAIN", "cell")                   # 'cell' (mclapply reps) | 'cellrep'
Nsel <- NS; Msel <- MS
nset <- Sys.getenv("REGIME_NSET", ""); if (nzchar(nset)) Nsel <- as.integer(strsplit(nset, ",")[[1]])
mset <- Sys.getenv("REGIME_MSET", ""); if (nzchar(mset)) Msel <- as.integer(strsplit(mset, ",")[[1]])
ldset <- Sys.getenv("REGIME_LDSET", ""); if (nzchar(ldset)) LDS <- strsplit(ldset, ",")[[1]]
drset <- Sys.getenv("REGIME_DRSET", ""); if (nzchar(drset)) DRS <- as.numeric(strsplit(drset, ",")[[1]])

cells <- expand.grid(N = Nsel, nMarker = Msel, LD = LDS, dr = DRS, stringsAsFactors = FALSE)
if (grain == "cellrep") {
  jobs <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
    data.frame(cells[i, ], rep = seq_len(reps_for(cells$N[i])))))
} else jobs <- cbind(cells, rep = NA_integer_)

if (Sys.getenv("DRYRUN") == "1") {
  Nsel <- c(150L,300L); jobs <- expand.grid(N = Nsel, nMarker = 1000L, LD = c("independent","high"),
            dr = 0, rep = NA_integer_, stringsAsFactors = FALSE); grain <- "cell"
  assign("reps_for", function(N) 2L, envir = environment())
}

## ---- TOPUP mode: run an explicit (N,M,LD,dr,rep) work-list ----------
## When TOPUP_JOBS points to a CSV, ignore the generated grid and run exactly
## those cell-reps (grain forced to cellrep). seed_of() is index-based on the
## fixed NS/MS/LDS vectors, so listed reps use the SAME deterministic seeds as
## the original run -- existing reps reproduce, new reps are clean independent
## draws. SLURM_ARRAY_TASK_ID + CHUNK slice the work-list as usual.
topf <- Sys.getenv("TOPUP_JOBS", "")
if (nzchar(topf)) {
  tj <- read.csv(topf, stringsAsFactors = FALSE)
  jobs <- data.frame(N = as.integer(tj$N), nMarker = as.integer(tj$nMarker),
                     LD = as.character(tj$LD), dr = as.numeric(tj$dr),
                     rep = as.integer(tj$rep), stringsAsFactors = FALSE)
  grain <- "cellrep"
  cat(sprintf("TOPUP mode: %d cell-reps from %s\n", nrow(jobs), topf))
}

AID   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
CHUNK <- as.integer(Sys.getenv("CHUNK", as.character(nrow(jobs))))
sel   <- ((AID - 1L) * CHUNK + 1L):min(AID * CHUNK, nrow(jobs)); sel <- sel[sel >= 1L & sel <= nrow(jobs)]
cat(sprintf("EXPANDED regime tier=%s grain=%s: %d jobs total, task %d handles rows %s (NC=%d, COLCAP=%d, KCAP=%d)\n",
    tier, grain, nrow(jobs), AID, paste(range(sel), collapse = "-"), NC, COLCAP, KCAP))

emit <- list()
for (i in sel) {
  jr <- jobs[i, ]
  if (grain == "cellrep") {
    emit[[length(emit) + 1L]] <- tryCatch(run_cellrep(jr$N, jr$nMarker, jr$LD, jr$dr, jr$rep),
                                          error = function(e) { cat("ERR cellrep at [", deparse(conditionCall(e)), "] ::", conditionMessage(e), "\n"); NULL })
  } else {
    R <- reps_for(jr$N)
    rr <- mclapply(1:R, function(rp) tryCatch(run_cellrep(jr$N, jr$nMarker, jr$LD, jr$dr, rp),
                   error = function(e) NULL), mc.cores = NC)
    emit[[length(emit) + 1L]] <- do.call(rbind, Filter(Negate(is.null), rr))
  }
}
res <- do.call(rbind, Filter(Negate(is.null), emit))
outf <- file.path(OUT, sprintf("regime_exp_t%s_%04d.csv", tier, AID))
write.csv(res, outf, row.names = FALSE)
cat(sprintf("wrote %d rows -> %s\n", if (is.null(res)) 0L else nrow(res), outf))
if (!is.null(res)) print(aggregate(spurious_Vd ~ N + structure, res, function(z) round(mean(z, na.rm = TRUE), 3)))
cat("EXPANDED REGIME DONE\n")
}  # end dispatcher (skipped when REGIME_FUNCS_ONLY=1)

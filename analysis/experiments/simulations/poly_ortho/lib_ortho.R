# experiments/simulations/poly_ortho/lib_ortho.R
# Shared helpers for the blueberry + regime-map expansion (idea:010).
# Mirrors the validated functions in real_ortho_roza.R / intervals_null.R / diagnostics.R.
.have_BGLR <- requireNamespace("BGLR", quietly = TRUE)
if (.have_BGLR) suppressPackageStartupMessages(library(BGLR))
.BGLR <- function(...) {
  if (!.have_BGLR) stop("install BGLR to run variance-component fitting", call. = FALSE)
  BGLR::BGLR(...)
}

# Orthonormal dosage contrasts under binomial genotype freqs at allele freq p.
# Cols 1..(ploidy+1) = degree 0..ploidy ; rows d=0..ploidy.
ortho <- function(p, ploidy = 4L) {
  d <- 0:ploidy; w <- choose(ploidy, d) * p^d * (1 - p)^(ploidy - d); w <- w / sum(w)
  inn <- function(a, b) sum(w * a * b); Q <- matrix(0, ploidy + 1, ploidy + 1)
  for (k in 1:(ploidy + 1)) {
    v <- d^(k - 1)
    if (k > 1) for (j in 1:(k - 1)) v <- v - inn(v, Q[, j]) / inn(Q[, j], Q[, j]) * Q[, j]
    nrm <- sqrt(inn(v, v)); Q[, k] <- if (nrm > 0) v / nrm else v
  }
  Q
}

# Build standard (A_VanRaden, D_digenic) + orthogonal (OA,OD2,OD3,OD4) kernels
# from a dosage matrix M (rows = individuals, cols = markers; values 0..ploidy).
build_kernels <- function(M, ploidy = 4L) {
  n <- nrow(M); m <- ncol(M)
  X <- replicate(4, matrix(0, n, m), simplify = FALSE); v <- 0L
  for (j in 1:m) {
    g <- round(M[, j]); p <- mean(g) / ploidy; if (p <= 0 || p >= 1) next
    v <- v + 1L; Q <- ortho(p, ploidy); for (dg in 1:4) X[[dg]][, j] <- Q[g + 1, dg + 1]
  }
  mk <- function(Z) { K <- tcrossprod(Z) / v; dimnames(K) <- list(rownames(M), rownames(M)); K / mean(diag(K)) + diag(1e-6, n) }
  p <- colMeans(M) / ploidy; Z <- sweep(M, 2, ploidy * p); dd <- sum(ploidy * p * (1 - p)); if (dd <= 0) dd <- 1
  KA <- tcrossprod(Z) / dd; KA <- KA / mean(diag(KA)) + diag(1e-6, n); dimnames(KA) <- list(rownames(M), rownames(M))
  Hh <- (M * (ploidy - M)) / choose(ploidy, 2); Hh <- sweep(Hh, 2, colMeans(Hh))
  KD <- tcrossprod(Hh); KD <- KD / mean(diag(KD)) + diag(1e-6, n); dimnames(KD) <- list(rownames(M), rownames(M))
  list(A = KA, D = KD, OA = mk(X[[1]]), OD2 = mk(X[[2]]), OD3 = mk(X[[3]]), OD4 = mk(X[[4]]), valid = v)
}

vget <- function(fm, nm) fm$ETA[[nm]]$varU
ms   <- function(fm, nm) c(m = fm$ETA[[nm]]$varU, s = max(fm$ETA[[nm]]$SD.varU, 1e-9))

fit_std  <- function(y, K, nIter = 8000, burnIn = 2000)
  .BGLR(y = y, ETA = list(A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")),
        nIter = nIter, burnIn = burnIn, verbose = FALSE, saveAt = file.path(tempdir(), "s_"))
fit_orth <- function(y, K, nIter = 8000, burnIn = 2000)
  .BGLR(y = y, ETA = list(A = list(K = K$OA, model = "RKHS"), D2 = list(K = K$OD2, model = "RKHS"),
        D3 = list(K = K$OD3, model = "RKHS"), D4 = list(K = K$OD4, model = "RKHS")),
        nIter = nIter, burnIn = burnIn, verbose = FALSE, saveAt = file.path(tempdir(), "o_"))

# Monte-Carlo credible intervals on variance fractions from per-kernel (mean, sd).
# NOTE: rough independent-normal approximation; prefer frac_ci_samples() (proper
# posterior CrIs from aligned MCMC samples). Kept for the within-blueberry N-curve trend.
mc_frac <- function(mss, D = 2e4) {
  X <- sapply(mss, function(p) pmax(rnorm(D, p["m"], p["s"]), 1e-9)); X <- X / rowSums(X)
  apply(X, 2, function(c) c(mean = mean(c), lo = quantile(c, .025), hi = quantile(c, .975)))
}

# PROPER posterior credible intervals on variance fractions, from BGLR's saved
# per-component varU samples (<pfx>ETA_<name>_varU.dat). Aligns the kernels'
# samples per iteration, drops burn-in, then takes quantiles of the fractions.
frac_ci_samples <- function(pfx, etanames, burnIn, thin) {
  S <- lapply(etanames, function(nm) scan(paste0(pfx, "ETA_", nm, "_varU.dat"), quiet = TRUE))
  L <- min(sapply(S, length)); S <- lapply(S, function(v) v[seq_len(L)])
  drop <- floor(burnIn / thin); keep <- (drop + 1):L
  M <- do.call(cbind, lapply(S, function(v) v[keep]))      # iterations x components
  Fr <- M / rowSums(M)
  out <- apply(Fr, 2, function(c) c(mean = mean(c), lo = quantile(c, .025, names = FALSE), hi = quantile(c, .975, names = FALSE)))
  colnames(out) <- etanames; out
}

# Fit a BGLR RKHS model and return the fit + proper fraction CrIs from samples.
run_model <- function(y, ETA, etanames, nIter = 12000, burnIn = 3000, thin = 5) {
  d <- tempfile("bg"); dir.create(d); pfx <- file.path(d, "x_")
  fm <- .BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, verbose = FALSE, saveAt = pfx)
  ci <- frac_ci_samples(pfx, etanames, burnIn, thin)
  unlink(d, recursive = TRUE)
  list(fm = fm, ci = ci)
}

# Collinearity diagnostics on a list of kernels.
dc   <- function(K) { n <- nrow(K); J <- diag(n) - 1 / n; J %*% K %*% J }
fcos <- function(A, B) { a <- dc(A); b <- dc(B); sum(a * b) / sqrt(sum(a * a) * sum(b * b)) }
# Effective rank / condition number from UNIT-FROBENIUS-NORMALIZED double-centred
# kernel vectors (Gram = cosine matrix -> scale-free collinearity). eff_rank in
# [1, length(Klist)]; ideal (mutually orthogonal) = length(Klist).
design_diag <- function(Klist) {
  vc <- sapply(Klist, function(K) { v <- as.vector(dc(K)); v / sqrt(sum(v * v)) })
  G <- crossprod(vc); ev <- eigen(G, only.values = TRUE)$values; ev <- pmax(ev, 0)
  list(cond = max(ev) / max(min(ev), 1e-12), eff_rank = sum(ev)^2 / sum(ev^2))
}

# 5-fold x reps CV predictive ability for a given ETA.
# set.seed(1000 + r) makes the fold partition deterministic per rep so that
# SEPARATE cv_pa() calls for competing models (A, A+D, orth) share IDENTICAL
# folds at each rep -- i.e. model deltas are PAIRED, matching run_pred_cv.R.
# Without this each model is scored on its own random folds and the delta
# conflates the model effect with fold-partition luck. (The cross-crop panels
# in the paper are produced by the paired run_pred_cv.R; this keeps the
# library helper consistent for any direct re-run of run_ortho_dataset.R.)
cv_pa <- function(y, ETA, folds = 5, reps = 5, nIter = 1500, burnIn = 400) {
  ok <- which(!is.na(y)); pa <- c()
  for (r in 1:reps) {
    set.seed(1000 + r); fold <- sample(rep(1:folds, length.out = length(ok))); pr <- rep(NA_real_, length(y))
    for (f in 1:folds) { te <- ok[fold == f]; yN <- y; yN[te] <- NA
      fm <- .BGLR(y = yN, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE,
                  saveAt = file.path(tempdir(), sprintf("cv_%d_%d_", r, f))); pr[te] <- fm$yHat[te] }
    pa <- c(pa, cor(pr[ok], y[ok]))
  }
  mean(pa)
}

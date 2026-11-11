#!/usr/bin/env Rscript
# =====================================================================
# Generalized arbitrary-ploidy kernel construction (m > 4). Sources lib_ortho.R
# for ortho(), run_model(), frac_ci_samples(), fcos(), design_diag() (all generic).
# Adds:
#   build_kernels_m(M, ploidy)        -- fixed ploidy m: A, D (digenic), and
#                                        orthogonal degree kernels O[[1..m]]
#                                        (O[[1]]=additive, O[[2]]=D2, ... O[[m]]=Dm)
#   build_kernels_varm(M, mploidy, maxdeg) -- per-marker variable ploidy (e.g. sugarcane):
#                                        contrasts built at each marker's ploidy,
#                                        degrees 1..maxdeg (default min ploidy).
# ortho(p, ploidy) already supports any ploidy (Gram-Schmidt under the binomial
# genotype-frequency measure), so this only assembles kernels for degrees 1..m.
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho.R")

build_kernels_m <- function(M, ploidy) {
  n <- nrow(M); mrk <- ncol(M)
  X <- lapply(1:ploidy, function(k) matrix(0, n, mrk)); v <- 0L
  for (j in 1:mrk) {
    g <- round(M[, j]); p <- mean(g) / ploidy; if (p <= 0 || p >= 1) next
    v <- v + 1L; Q <- ortho(p, ploidy)
    for (deg in 1:ploidy) X[[deg]][, j] <- Q[g + 1, deg + 1]
  }
  mk <- function(Z) { K <- tcrossprod(Z) / v; dimnames(K) <- list(rownames(M), rownames(M)); K / mean(diag(K)) + diag(1e-6, n) }
  p <- colMeans(M) / ploidy; Z <- sweep(M, 2, ploidy * p); dd <- sum(ploidy * p * (1 - p)); if (dd <= 0) dd <- 1
  KA <- tcrossprod(Z) / dd; KA <- KA / mean(diag(KA)) + diag(1e-6, n); dimnames(KA) <- list(rownames(M), rownames(M))
  Hh <- (M * (ploidy - M)) / choose(ploidy, 2); Hh <- sweep(Hh, 2, colMeans(Hh))
  KD <- tcrossprod(Hh); KD <- KD / mean(diag(KD)) + diag(1e-6, n); dimnames(KD) <- list(rownames(M), rownames(M))
  list(A = KA, D = KD, O = lapply(1:ploidy, function(k) mk(X[[k]])), ploidy = ploidy, valid = v)
}

build_kernels_varm <- function(M, mploidy, maxdeg = NULL) {
  n <- nrow(M); mrk <- ncol(M)
  stopifnot(length(mploidy) == mrk)
  if (is.null(maxdeg)) maxdeg <- min(mploidy)
  MP <- matrix(mploidy, n, mrk, byrow = TRUE)
  X <- lapply(1:maxdeg, function(k) matrix(0, n, mrk)); v <- 0L
  for (j in 1:mrk) {
    pj <- mploidy[j]; g <- round(M[, j]); p <- mean(g) / pj; if (p <= 0 || p >= 1) next
    v <- v + 1L; Q <- ortho(p, pj)
    for (deg in 1:maxdeg) if (deg <= pj) X[[deg]][, j] <- Q[g + 1, deg + 1]
  }
  mk <- function(Z) { K <- tcrossprod(Z) / v; dimnames(K) <- list(rownames(M), rownames(M)); K / mean(diag(K)) + diag(1e-6, n) }
  p <- colMeans(M) / mploidy; Z <- M - MP * matrix(p, n, mrk, byrow = TRUE); dd <- sum(mploidy * p * (1 - p)); if (dd <= 0) dd <- 1
  KA <- tcrossprod(Z) / dd; KA <- KA / mean(diag(KA)) + diag(1e-6, n); dimnames(KA) <- list(rownames(M), rownames(M))
  Hh <- (M * (MP - M)) / matrix(choose(mploidy, 2), n, mrk, byrow = TRUE); Hh <- sweep(Hh, 2, colMeans(Hh))
  KD <- tcrossprod(Hh); KD <- KD / mean(diag(KD)) + diag(1e-6, n); dimnames(KD) <- list(rownames(M), rownames(M))
  list(A = KA, D = KD, O = lapply(1:maxdeg, function(k) mk(X[[k]])), ploidy = paste0("var<=", maxdeg), valid = v)
}

# Convenience: ETA name vector for an orthogonal model with `ndeg` degrees (incl. additive).
ortho_etanames <- function(ndeg) c("A", paste0("D", 2:ndeg))

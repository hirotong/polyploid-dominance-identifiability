#' Names of the orthogonal degree-wise components
#'
#' @param ndeg number of degrees to include (\code{ndeg >= 1}); degree 1 is
#'   additive, degrees 2..\code{ndeg} are digenic, trigenic, ....
#' @return a character vector \code{c("A", "D2", ..., "D{ndeg}")}.
#' @examples
#' ortho_etanames(4L)
#' @export
ortho_etanames <- function(ndeg) c("A", paste0("D", 2:ndeg))

#' Orthogonal degree-wise relationship kernels (A, D2, D3, D4, ...) from dosage
#'
#' Builds one relationship kernel per orthogonal contrast degree
#' (\code{\link{ortho_contrasts}}): additive (A), digenic (D2), trigenic (D3),
#' quadrigenic (D4) and so on up to \code{ploidy}. Each kernel is a
#' cross-product of the per-marker degree contrasts, scaled to mean diagonal 1
#' with a small ridge for numerical stability.
#'
#' @param M individuals x markers dosage matrix, values in \code{0..ploidy}.
#' @param ploidy even ploidy (default 4).
#' @param impute if \code{TRUE}, mean-impute missing dosage calls instead of
#'   erroring; see \code{\link{validate_dosage_matrix}}.
#' @param min_maf minimum marker minor-allele frequency to retain.
#' @param min_call_rate minimum marker call rate to retain.
#' @return named list of \code{n x n} kernels: \code{A} (additive), \code{D2}
#'   (digenic), ..., \code{D{ploidy}}, each with mean diagonal 1 plus a
#'   \code{1e-6} ridge and row/column names matching \code{rownames(M)}.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KO <- orthogonal_kernels(M)
#' names(KO)
#' @export
orthogonal_kernels <- function(M, ploidy = 4L, impute = FALSE, min_maf = 0, min_call_rate = 0) {
  chk <- validate_dosage_matrix(M, ploidy, impute = impute, min_maf = min_maf, min_call_rate = min_call_rate)
  G <- chk$G; p <- chk$p; ploidy <- chk$ploidy
  n <- nrow(G); m <- ncol(G)
  X <- replicate(ploidy, matrix(0, n, m), simplify = FALSE)
  for (j in seq_len(m)) {
    g <- G[, j]; Q <- ortho_contrasts(p[j], ploidy)
    for (deg in seq_len(ploidy)) X[[deg]][, j] <- Q[g + 1L, deg + 1L]
  }
  mkK <- function(Z) { K <- tcrossprod(Z) / m; K / mean(diag(K)) + diag(1e-6, nrow(Z)) }
  ks <- lapply(X, mkK)
  nm <- ortho_etanames(ploidy)
  stats::setNames(
    lapply(seq_along(ks), function(i) { k <- ks[[i]]; dimnames(k) <- list(rownames(G), rownames(G)); k }),
    nm)
}

#' Standard additive and dominance kernels from dosage
#'
#' VanRaden-style additive kernel plus the digenic dominance (heterozygosity)
#' kernel commonly used as the "standard" A/D pair, for comparison against
#' \code{\link{orthogonal_kernels}}.
#'
#' @param M individuals x markers dosage matrix, values in \code{0..ploidy}.
#' @param ploidy even ploidy (default 4).
#' @param impute if \code{TRUE}, mean-impute missing dosage calls instead of
#'   erroring; see \code{\link{validate_dosage_matrix}}.
#' @param min_maf minimum marker minor-allele frequency to retain.
#' @param min_call_rate minimum marker call rate to retain.
#' @return a list with \code{A} (additive) and \code{D} (digenic dominance)
#'   \code{n x n} kernels, each with mean diagonal 1 plus a \code{1e-6} ridge.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KS <- standard_kernels(M)
#' names(KS)
#' @export
standard_kernels <- function(M, ploidy = 4L, impute = FALSE, min_maf = 0, min_call_rate = 0) {
  chk <- validate_dosage_matrix(M, ploidy, impute = impute, min_maf = min_maf, min_call_rate = min_call_rate)
  M <- chk$M; ploidy <- chk$ploidy
  p <- colMeans(M) / ploidy; Z <- sweep(M, 2, ploidy * p)
  dd <- sum(ploidy * p * (1 - p)); if (dd <= 0) dd <- 1
  A <- tcrossprod(Z) / dd; A <- A / mean(diag(A)) + diag(1e-6, nrow(M))
  H <- (M * (ploidy - M)) / choose(ploidy, 2); H <- sweep(H, 2, colMeans(H))
  D <- tcrossprod(H); D <- D / mean(diag(D)) + diag(1e-6, nrow(M))
  dimnames(A) <- dimnames(D) <- list(rownames(M), rownames(M))
  list(A = A, D = D)
}

#' Panel kernels: standard A/D plus orthogonal degree kernels and contrasts
#'
#' Builds the standard additive (A) and dominance (D) kernels together with
#' the orthogonal degree-wise kernels \code{O$A, O$D2, ..., O$D{ndeg}}, from a
#' fixed- or per-marker variable-ploidy panel (\code{\link{parse_ploidy}}).
#' When \code{need_X = TRUE}, also returns the underlying \code{n x m}
#' orthonormal contrast matrices \code{X[[1..ndeg]]} (degree 1 = additive), so
#' downstream simulators can draw marker effects in the identical contrast
#' space the kernels are built from. Normalisation matches
#' \code{\link{orthogonal_kernels}} and \code{\link{standard_kernels}} exactly:
#' \code{K / mean(diag(K)) + diag(1e-6, n)}.
#'
#' @param M individuals x markers dosage matrix; \code{NA} not allowed (QC it
#'   first with \code{\link{filter_markers}}).
#' @param pcfg a \code{\link{parse_ploidy}} result, fixed or variable mode.
#' @param need_X if \code{TRUE} (default), also return the per-degree
#'   contrast matrices \code{X}; set \code{FALSE} to skip them when only the
#'   kernels are needed.
#' @return a list with \code{A}, \code{D} (\code{n x n} standard kernels),
#'   \code{O} (named list of \code{n x n} orthogonal degree kernels,
#'   \code{\link{ortho_etanames}(pcfg$ndeg)}), \code{valid} (number of
#'   polymorphic markers used), \code{ndeg}, \code{etanames}, and — when
#'   \code{need_X = TRUE} — \code{X} (list of \code{ndeg} \code{n x m}
#'   contrast matrices).
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)))
#' names(KP$O)
#' @export
build_panel <- function(M, pcfg, need_X = TRUE) {
  n <- nrow(M); mrk <- ncol(M)
  mk <- function(Z, v) {
    K <- tcrossprod(Z) / v; dimnames(K) <- list(rownames(M), rownames(M))
    K / mean(diag(K)) + diag(1e-6, n)
  }
  ndeg <- pcfg$ndeg
  X <- lapply(seq_len(ndeg), function(k) matrix(0, n, mrk)); v <- 0L

  if (pcfg$mode == "fixed") {
    pl <- pcfg$ploidy
    for (j in seq_len(mrk)) {
      g <- round(M[, j]); p <- mean(g) / pl; if (p <= 0 || p >= 1) next
      v <- v + 1L; Q <- ortho_contrasts(p, pl)
      for (deg in seq_len(ndeg)) X[[deg]][, j] <- Q[g + 1, deg + 1]
    }
    p <- colMeans(M) / pl; Z <- sweep(M, 2, pl * p); dd <- sum(pl * p * (1 - p)); if (dd <= 0) dd <- 1
    KA <- tcrossprod(Z) / dd; KA <- KA / mean(diag(KA)) + diag(1e-6, n); dimnames(KA) <- list(rownames(M), rownames(M))
    Hh <- (M * (pl - M)) / choose(pl, 2); Hh <- sweep(Hh, 2, colMeans(Hh))
    KD <- tcrossprod(Hh); KD <- KD / mean(diag(KD)) + diag(1e-6, n); dimnames(KD) <- list(rownames(M), rownames(M))
  } else {
    mp <- pcfg$mploidy; MP <- matrix(mp, n, mrk, byrow = TRUE)
    for (j in seq_len(mrk)) {
      pj <- mp[j]; g <- round(M[, j]); p <- mean(g) / pj; if (p <= 0 || p >= 1) next
      v <- v + 1L; Q <- ortho_contrasts(p, pj)
      for (deg in seq_len(ndeg)) if (deg <= pj) X[[deg]][, j] <- Q[g + 1, deg + 1]
    }
    p <- colMeans(M) / mp; Z <- M - MP * matrix(p, n, mrk, byrow = TRUE); dd <- sum(mp * p * (1 - p)); if (dd <= 0) dd <- 1
    KA <- tcrossprod(Z) / dd; KA <- KA / mean(diag(KA)) + diag(1e-6, n); dimnames(KA) <- list(rownames(M), rownames(M))
    Hh <- (M * (MP - M)) / matrix(choose(mp, 2), n, mrk, byrow = TRUE); Hh <- sweep(Hh, 2, colMeans(Hh))
    KD <- tcrossprod(Hh); KD <- KD / mean(diag(KD)) + diag(1e-6, n); dimnames(KD) <- list(rownames(M), rownames(M))
  }
  if (v < 1L) stop("no polymorphic markers after filtering", call. = FALSE)
  O <- lapply(seq_len(ndeg), function(k) mk(X[[k]], v))
  names(O) <- ortho_etanames(ndeg)                          # A, D2, ..., D{ndeg}
  out <- list(A = KA, D = KD, O = O, valid = v, ndeg = ndeg, etanames = ortho_etanames(ndeg))
  if (need_X) out$X <- X                                    # keep contrasts only when simulating
  out
}

#' Degree-1 and degree-2 orthonormal contrast matrices only
#'
#' Lightweight companion to \code{\link{build_panel}} that builds only the
#' additive (degree-1) and digenic (degree-2) orthonormal contrast matrices,
#' skipping the \code{n x n} kernel cross-products. Useful for simulators that
#' only need low-degree marker effects.
#'
#' @param M individuals x markers dosage matrix; \code{NA} not allowed.
#' @param pcfg a \code{\link{parse_ploidy}} result, fixed or variable mode.
#' @return a list with \code{X} (a length-2 list of \code{n x m} contrast
#'   matrices, degree 1 then degree 2) and \code{valid} (number of
#'   polymorphic markers used).
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' cc <- build_contrasts12(M, parse_ploidy("4", ncol(M)))
#' length(cc$X)
#' @export
build_contrasts12 <- function(M, pcfg) {
  n <- nrow(M); mrk <- ncol(M)
  X1 <- matrix(0, n, mrk); X2 <- matrix(0, n, mrk); v <- 0L
  if (pcfg$mode == "fixed") {
    pl <- pcfg$ploidy
    for (j in seq_len(mrk)) {
      g <- round(M[, j]); p <- mean(g) / pl; if (p <= 0 || p >= 1) next
      v <- v + 1L; Q <- ortho_contrasts(p, pl); X1[, j] <- Q[g + 1, 2]; X2[, j] <- Q[g + 1, 3]
    }
  } else {
    mp <- pcfg$mploidy
    for (j in seq_len(mrk)) {
      pj <- mp[j]; g <- round(M[, j]); p <- mean(g) / pj; if (p <= 0 || p >= 1) next
      v <- v + 1L; Q <- ortho_contrasts(p, pj); X1[, j] <- Q[g + 1, 2]; if (pj >= 2) X2[, j] <- Q[g + 1, 3]
    }
  }
  list(X = list(X1, X2), valid = v)
}

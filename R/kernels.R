## Internal: cross-product relationship kernel, renormalised to mean diagonal 1
## with a 1e-6 ridge. Every polyDx kernel is renormalised this way, so any
## positive scalar denominator applied to the cross-product (VanRaden's
## sum(ploidy * p * (1 - p)), a marker count, ...) cancels exactly and is
## deliberately not computed.
kernel_from <- function(Z, ids = rownames(Z), what = "kernel") {
  K <- tcrossprod(Z)
  md <- mean(diag(K))
  # A panel with no between-individual variation in this component (e.g. every
  # tetraploid dosage is 1 or 3, so heterozygosity is constant) makes Z all
  # zeros. Dividing by md == 0 used to return an entirely NaN kernel, which
  # propagated as AD_cosine = NaN and only died later inside BGLR.
  if (!is.finite(md) || md <= 0) {
    stop(sprintf("%s is degenerate: no between-individual variation to normalise (mean diagonal %s)",
                 what, format(md)), call. = FALSE)
  }
  K <- K / md + diag(1e-6, nrow(Z))
  dimnames(K) <- list(ids, ids)
  K
}

## Internal: per-marker orthonormal degree contrasts. Returns `ndeg` n x m
## matrices (degree 1 = additive, 2 = digenic, ...) plus the number of
## polymorphic markers used; monomorphic markers are left as zero columns.
## `mp` is the per-marker ploidy vector (length ncol(M)). Shared by
## build_panel(), build_contrasts12() and orthogonal_kernels().
degree_contrasts <- function(M, mp, ndeg) {
  X <- lapply(seq_len(ndeg), function(k) matrix(0, nrow(M), ncol(M)))
  v <- 0L
  for (j in seq_len(ncol(M))) {
    g <- round(M[, j]); p <- mean(g) / mp[j]
    if (p <= 0 || p >= 1) next                       # monomorphic: leave zeros
    v <- v + 1L; Q <- ortho_contrasts(p, mp[j])
    for (deg in seq_len(ndeg)) X[[deg]][, j] <- Q[g + 1L, deg + 1L]
  }
  list(X = X, valid = v)
}

#' Names of the orthogonal degree-wise components
#'
#' @param ndeg number of degrees to include (\code{ndeg >= 1}); degree 1 is
#'   additive, degrees 2..\code{ndeg} are digenic, trigenic, ....
#' @return a character vector \code{c("A", "D2", ..., "D{ndeg}")}.
#' @examples
#' ortho_etanames(4L)
#' @export
ortho_etanames <- function(ndeg) {
  if (!is.numeric(ndeg) || length(ndeg) != 1L || !is.finite(ndeg) || ndeg < 1) {
    stop("ndeg must be a single integer value >= 1", call. = FALSE)
  }
  if (ndeg < 2) "A" else c("A", paste0("D", 2:ndeg))   # 2:1 would count backwards
}

#' Orthogonal degree-wise relationship kernels (A, D2, D3, D4, ...) from dosage
#'
#' Builds one relationship kernel per orthogonal contrast degree
#' (\code{\link{ortho_contrasts}}): additive (A), digenic (D2), trigenic (D3),
#' quadrigenic (D4) and so on up to \code{ploidy}. Each kernel is a
#' cross-product of the per-marker degree contrasts, scaled to mean diagonal 1
#' with a small ridge for numerical stability.
#'
#' @param M individuals x markers dosage matrix, values in \code{0..ploidy}.
#' @param ploidy the panel's ploidy. Required, with no default: an unstated
#'   ploidy is silently wrong rather than noisily wrong (a diploid matrix
#'   analysed as tetraploid halves every allele frequency, makes the
#'   dominance kernel monotone in dosage instead of peaked at the
#'   heterozygote, and inflates the reported A-D cosine).
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
#' KO <- orthogonal_kernels(M, ploidy = 4L)
#' names(KO)
#' @export
orthogonal_kernels <- function(M, ploidy, impute = FALSE, min_maf = 0, min_call_rate = 0) {
  chk <- validate_dosage_matrix(M, ploidy, impute = impute, min_maf = min_maf, min_call_rate = min_call_rate)
  G <- chk$G; ploidy <- chk$ploidy
  X <- degree_contrasts(G, rep(ploidy, ncol(G)), ploidy)$X
  nm <- ortho_etanames(ploidy)
  stats::setNames(Map(function(Z, k) kernel_from(Z, rownames(G), paste("degree kernel", k)),
                      X, nm), nm)
}

#' Standard additive and dominance kernels from dosage
#'
#' VanRaden-style additive kernel plus the digenic dominance (heterozygosity)
#' kernel commonly used as the "standard" A/D pair, for comparison against
#' \code{\link{orthogonal_kernels}}.
#'
#' @param M individuals x markers dosage matrix, values in \code{0..ploidy}.
#' @param ploidy the panel's ploidy. Required, with no default; see
#'   \code{\link{orthogonal_kernels}}.
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
#' KS <- standard_kernels(M, ploidy = 4L)
#' names(KS)
#' @export
standard_kernels <- function(M, ploidy, impute = FALSE, min_maf = 0, min_call_rate = 0) {
  chk <- validate_dosage_matrix(M, ploidy, impute = impute, min_maf = min_maf, min_call_rate = min_call_rate)
  M <- chk$M; ploidy <- chk$ploidy
  p <- colMeans(M) / ploidy
  H <- (M * (ploidy - M)) / choose(ploidy, 2)
  list(A = kernel_from(sweep(M, 2, ploidy * p), rownames(M), "additive kernel A"),
       D = kernel_from(sweep(H, 2, colMeans(H)), rownames(M), "dominance kernel D"))
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
  n <- nrow(M); mrk <- ncol(M); ndeg <- pcfg$ndeg
  mp <- panel_mploidy(M, pcfg, ndeg)              # fixed ploidy = a constant vector
  MP <- matrix(mp, n, mrk, byrow = TRUE)

  dc <- degree_contrasts(M, mp, ndeg); X <- dc$X; v <- dc$valid
  if (v < 1L) stop("no polymorphic markers after filtering", call. = FALSE)

  p  <- colMeans(M) / mp
  KA <- kernel_from(M - MP * matrix(p, n, mrk, byrow = TRUE), rownames(M), "additive kernel A")
  Hh <- (M * (MP - M)) / matrix(choose(mp, 2), n, mrk, byrow = TRUE)
  KD <- kernel_from(sweep(Hh, 2, colMeans(Hh)), rownames(M), "dominance kernel D")

  O <- Map(function(Z, nm) kernel_from(Z, rownames(M), paste("degree kernel", nm)),
           X, ortho_etanames(ndeg))
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
  degree_contrasts(M, panel_mploidy(M, pcfg, 2L), 2L)
}

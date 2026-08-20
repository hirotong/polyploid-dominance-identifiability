## Internal kernel-geometry primitives. Not exported: experiments/simulations/
## poly_ortho/lib_ortho.R defines dc()/fcos()/design_diag() into the global
## environment, and the Task 8 shim sources both that file and this package;
## identical names would collide, so these are renamed and kept internal.

## Double-center a kernel: J K J with J = I - 1/n.
dcenter <- function(K) {
  n <- nrow(K); J <- diag(n) - 1 / n
  J %*% K %*% J
}

## Frobenius cosine between two double-centred kernels.
fcosine <- function(A, B) {
  a <- dcenter(A); b <- dcenter(B)
  sum(a * b) / sqrt(sum(a * a) * sum(b * b))
}

## Effective rank / condition number from unit-Frobenius-normalised
## double-centred kernel vectors (Gram = cosine matrix -> scale-free
## collinearity). eff_rank in [1, length(Klist)]; ideal (mutually orthogonal)
## = length(Klist).
design_diag_polydx <- function(Klist) {
  vc <- vapply(Klist, function(K) { v <- as.vector(dcenter(K)); v / sqrt(sum(v * v)) },
               numeric(length(Klist[[1]])))
  G <- crossprod(vc); ev <- eigen(G, only.values = TRUE)$values; ev <- pmax(ev, 0)
  list(cond = max(ev) / max(min(ev), 1e-12), eff_rank = sum(ev)^2 / sum(ev^2))
}

#' Kernel-geometry diagnostics (phenotype-free)
#'
#' Standard additive-dominance (A-D) cosine, the full degree-wise cosine
#' matrix over the orthogonal kernels, and the effective rank / condition
#' number of the degree design. These summarise how collinear a panel's
#' kernels are, independent of any fitted model.
#'
#' @param KP a \code{\link{build_panel}} (or \code{\link{orthogonal_kernels}}
#'   plus \code{A}/\code{D}) result: a list with \code{A}, \code{D} (standard
#'   kernels) and \code{O} (named list of orthogonal degree kernels).
#' @return a list with \code{AD_cosine} (Frobenius cosine between \code{A} and
#'   \code{D}), \code{cosine_matrix} (\code{ndeg x ndeg} cosine matrix over
#'   \code{O}), \code{max_abs_degree_cosine} (largest off-diagonal absolute
#'   cosine in \code{cosine_matrix}), \code{eff_rank}, \code{ideal_rank}
#'   (\code{length(O)}), \code{rank_ratio} (\code{eff_rank / ideal_rank}) and
#'   \code{condition_number}.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(30 * 50, 4, 0.4), 30, 50)
#' rownames(M) <- paste0("ind", seq_len(30))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = FALSE)
#' g <- kernel_geometry(KP)
#' g$AD_cosine
#' @export
kernel_geometry <- function(KP) {
  O <- KP$O; deg_names <- names(O)
  C <- matrix(NA_real_, length(O), length(O), dimnames = list(deg_names, deg_names))
  for (i in seq_along(O)) for (j in seq_along(O)) C[i, j] <- fcosine(O[[i]], O[[j]])
  dd <- design_diag_polydx(O)
  off <- C[upper.tri(C)]
  list(AD_cosine = fcosine(KP$A, KP$D),
       cosine_matrix = C,
       max_abs_degree_cosine = if (length(off)) max(abs(off), na.rm = TRUE) else NA_real_,
       eff_rank = dd$eff_rank, ideal_rank = length(O),
       rank_ratio = dd$eff_rank / length(O), condition_number = dd$cond)
}

#' Validate panel kernels before fitting
#'
#' Pre-fit sanity check on a panel's kernels: row/column ID order must match
#' the phenotype's ID order, kernels must contain no missing values and be
#' symmetric. Also warns (does not stop) when a kernel's mean diagonal drifts
#' from 1 or it is not positive semi-definite to tolerance.
#'
#' @param KP a \code{\link{build_panel}} result: a list with \code{A},
#'   \code{D} and \code{O} (named list of orthogonal degree kernels).
#' @param ids character vector of individual IDs in the phenotype's order;
#'   every kernel's row and column names must match this exactly.
#' @param tol eigenvalue tolerance for the positive-semi-definiteness warning.
#' @return invisibly \code{TRUE} if all checks pass; stops with an
#'   informative message on ID-order mismatch, missing values, or asymmetry.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(30 * 50, 4, 0.4), 30, 50)
#' rownames(M) <- paste0("ind", seq_len(30))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = FALSE)
#' check_kernels(KP, rownames(M))
#' @export
check_kernels <- function(KP, ids, tol = 1e-6) {
  Ks <- c(list(A = KP$A, D = KP$D), KP$O)
  for (nm in names(Ks)) {
    K <- Ks[[nm]]
    if (!identical(rownames(K), ids) || !identical(colnames(K), ids)) {
      stop(sprintf("kernel %s row/col order does not match phenotype id order", nm), call. = FALSE)
    }
    if (anyNA(K)) stop(sprintf("kernel %s has missing values", nm), call. = FALSE)
    if (max(abs(K - t(K))) > 1e-8) stop(sprintf("kernel %s not symmetric", nm), call. = FALSE)
    if (abs(mean(diag(K)) - 1) > 0.05) warning(sprintf("kernel %s mean diagonal %.3f not ~1", nm, mean(diag(K))))
    ev <- min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
    if (ev < -tol) warning(sprintf("kernel %s min eigenvalue %.2e < 0 (not PSD to tol)", nm, ev))
  }
  invisible(TRUE)
}

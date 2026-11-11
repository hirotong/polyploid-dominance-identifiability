#' Orthonormal dosage contrasts under binomial genotype frequencies
#'
#' Gram-Schmidt orthonormalisation of the dosage monomials
#' \eqn{1, d, d^2, \ldots, d^{ploidy}} under the binomial genotype weights at
#' allele frequency \code{p}. Column \code{k + 1} is the degree-\code{k}
#' contrast: 0 = mean, 1 = additive, 2 = digenic, 3 = trigenic, 4 = quadrigenic.
#'
#' @param p alternate-allele frequency at the locus, strictly inside (0, 1);
#'   the binomial weights degenerate at 0 and 1, so those are an error.
#' @param ploidy the locus ploidy. Required, with no default: an unstated
#'   ploidy silently rescales the binomial weights and every contrast.
#' @return a \code{(ploidy + 1) x (ploidy + 1)} matrix \code{Q}; rows are dosage
#'   \code{0..ploidy}, columns are degree \code{0..ploidy}, orthonormal under the
#'   binomial weights.
#' @examples
#' Q <- ortho_contrasts(0.3, ploidy = 4L)
#' round(Q[, 2], 3)  # additive contrast
#' @export
ortho_contrasts <- function(p, ploidy) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be a single allele frequency strictly inside (0, 1)", call. = FALSE)
  }
  if (!is.numeric(ploidy) || length(ploidy) != 1L || !is.finite(ploidy) ||
      ploidy < 1 || ploidy != as.integer(ploidy)) {
    stop("ploidy must be a single integer value >= 1", call. = FALSE)
  }
  d <- 0:ploidy
  w <- choose(ploidy, d) * p^d * (1 - p)^(ploidy - d); w <- w / sum(w)
  inn <- function(a, b) sum(w * a * b)
  Q <- matrix(0, ploidy + 1L, ploidy + 1L)
  for (k in 1:(ploidy + 1L)) {
    v <- d^(k - 1L)
    if (k > 1L) for (j in 1:(k - 1L)) v <- v - inn(v, Q[, j]) / inn(Q[, j], Q[, j]) * Q[, j]
    nrm <- sqrt(inn(v, v)); Q[, k] <- if (nrm > 0) v / nrm else v
  }
  Q
}

#' Mean-impute missing dosage calls by marker
#'
#' @param M individuals (rows) x markers (columns) dosage matrix.
#' @return a numeric matrix with \code{NA} replaced by the marker mean.
#' @examples
#' M <- matrix(c(0, 2, 4, NA, 1, 3), nrow = 3)
#' impute_mean_dosage(M)
#' @export
impute_mean_dosage <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "numeric"
  cm <- colMeans(M, na.rm = TRUE)
  bad <- !is.finite(cm)
  if (any(bad)) stop("cannot impute markers with all missing values")
  for (j in which(colSums(is.na(M)) > 0L)) M[is.na(M[, j]), j] <- cm[j]
  M
}

#' Validate and filter a dosage matrix before kernel construction
#'
#' Checks dosage range, ploidy, missingness and polymorphism, optionally
#' mean-imputes, and drops monomorphic markers and markers failing the MAF /
#' call-rate filters. A marker is monomorphic when its rounded dosage takes a
#' single value -- note this drops an all-heterozygous marker, whose allele
#' frequency is 0.5 but whose contrasts are constant. Orthogonal dosage contrasts are defined on integer dosage
#' classes, so the returned \code{G} is the rounded matrix used for contrast
#' lookup; \code{M} keeps the (possibly fractional, imputed) values.
#'
#' @param M individuals x markers dosage matrix in \code{0..ploidy}.
#' @param ploidy single integer ploidy, at least 2.
#' @param impute if \code{TRUE}, mean-impute missing calls instead of erroring.
#' @param min_maf minimum minor-allele frequency to retain a marker, in [0, 0.5].
#' @param min_call_rate minimum pre-imputation call rate to retain a marker, in [0, 1].
#' @return a list with \code{M} (retained markers, imputed), \code{G} (rounded),
#'   \code{p} (allele frequencies), \code{keep} (logical over the input columns)
#'   and \code{ploidy}.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(200, 4, 0.4), nrow = 20)
#' str(validate_dosage_matrix(M, ploidy = 4)$p)
#' @export
validate_dosage_matrix <- function(M, ploidy, impute = FALSE, min_maf = 0, min_call_rate = 0) {
  M <- as.matrix(M)
  storage.mode(M) <- "numeric"
  if (!is.numeric(ploidy) || length(ploidy) != 1L || !is.finite(ploidy) ||
      ploidy < 2 || ploidy != as.integer(ploidy)) {
    stop("ploidy must be a single integer value >= 2")
  }
  ploidy <- as.integer(ploidy)
  if (!is.numeric(min_maf) || length(min_maf) != 1L || !is.finite(min_maf) ||
      min_maf < 0 || min_maf > 0.5) {
    stop("min_maf must be a single finite value in [0, 0.5]")
  }
  if (!is.numeric(min_call_rate) || length(min_call_rate) != 1L ||
      !is.finite(min_call_rate) || min_call_rate < 0 || min_call_rate > 1) {
    stop("min_call_rate must be a single finite value in [0, 1]")
  }
  if (!impute && anyNA(M)) {
    stop("M contains missing values; mean-impute first with impute_mean_dosage(M) or call with impute=TRUE")
  }
  call_rate <- colMeans(!is.na(M))
  if (impute && anyNA(M)) M <- impute_mean_dosage(M)
  if (any(!is.finite(M))) stop("M must contain only finite dosage values after imputation")
  if (any(M < 0 | M > ploidy)) stop("M contains dosage values outside 0..ploidy")
  G <- round(M)
  p <- colMeans(G) / ploidy
  maf <- pmin(p, 1 - p)
  # p in (0, 1) is NOT polymorphism: an all-heterozygous tetraploid marker has
  # p = 0.5 and MAF 0.5 yet zero dosage variance, and its degree contrasts are
  # a constant NON-zero value, so keeping it adds a rank-1 constant matrix to
  # every degree kernel and over-counts `valid`. Test the dosage classes.
  varies <- apply(G, 2, function(z) length(unique(z)) > 1L)
  keep <- is.finite(p) & p > 0 & p < 1 & varies &
    maf >= min_maf & call_rate >= min_call_rate
  if (!any(keep)) stop("no polymorphic markers remain after filtering")
  list(M = M[, keep, drop = FALSE], G = G[, keep, drop = FALSE], p = p[keep],
       keep = keep, ploidy = ploidy)
}

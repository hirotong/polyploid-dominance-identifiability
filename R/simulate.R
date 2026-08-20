## Panel-conditioned simulators for the known-truth diagnostic controls.
## Ported verbatim (draw order preserved) from
## experiments/empirical/multi/lib_panel_workflow.R (sim_comp, gen_additive_only,
## gen_d2_spikein, gen_structure) and, for simulate_pure_degree, ported from
## software/polyortho.R -- removed in this release; see git history for its source.

#' Standardized degree-\emph{d} genetic value in a panel's own contrast space
#'
#' Draws a random marker-effect vector \code{b ~ N(0, 1/valid)} (one effect
#' per column of \code{Xd}), forms the genetic value \code{Xd \%*\% b}, and
#' standardizes it to mean 0, variance 1. \code{Xd} is one of the orthonormal
#' per-degree contrast matrices from \code{\link{build_panel}} (\code{KP$X}),
#' so the resulting genetic value lives in exactly the same contrast space the
#' panel's kernels are built from.
#'
#' @param Xd an n x m orthonormal degree contrast matrix (e.g. \code{KP$X[[1]]}
#'   for additive, \code{KP$X[[2]]} for digenic).
#' @param valid number of polymorphic markers used to build \code{Xd} (e.g.
#'   \code{KP$valid}); sets the marker-effect draw variance to \code{1/valid}.
#' @return a length-\code{nrow(Xd)} numeric vector, standardized to mean 0, sd
#'   1, or all zeros when the drawn genetic value has zero variance.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
#' g <- sim_comp(KP$X[[1]], KP$valid)
#' c(mean(g), sd(g))
#' @export
sim_comp <- function(Xd, valid) {
  b <- stats::rnorm(ncol(Xd), 0, sqrt(1 / valid))   # b ~ N(0, 1/M)
  g <- as.numeric(Xd %*% b)
  s <- stats::sd(g); if (!is.finite(s) || s == 0) return(rep(0, length(g)))
  as.numeric(scale(g))
}

#' Additive-only known-truth phenotype (true digenic fraction 0)
#'
#' Draws an additive genetic value with \code{\link{sim_comp}} and composes
#' \code{y = sqrt(H2) g_A + sqrt(1 - H2) e} with fresh Gaussian noise \code{e}.
#' The true digenic (D2) fraction of this phenotype is exactly 0, making it
#' the workflow's negative control.
#'
#' @param KP a \code{\link{build_panel}} result built with \code{need_X =
#'   TRUE} (needs \code{X[[1]]} and \code{valid}).
#' @param H2 broad-sense heritability in \code{[0, 1]}; defaults to
#'   \code{\link{PWF_H2_DEFAULT}}.
#' @return a length-n numeric phenotype vector.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
#' y <- gen_additive_only(KP, H2 = 0.6)
#' length(y)
#' @export
gen_additive_only <- function(KP, H2 = PWF_H2_DEFAULT) {
  gA <- sim_comp(KP$X[[1]], KP$valid)
  n <- length(gA); e <- stats::rnorm(n)
  as.numeric(sqrt(H2) * as.numeric(scale(gA)) + sqrt(1 - H2) * e)
}

#' Digenic (D2) spike-in known-truth phenotype
#'
#' Draws an additive genetic value \code{g_A} and a digenic genetic value
#' \code{g_D2} (both via \code{\link{sim_comp}}), residualises \code{g_D2} on
#' \code{g_A} so the planted digenic component is orthogonal to the additive
#' one, and composes the genetic value \code{G = sqrt(1 - s) g_A + sqrt(s)
#' rD2}. The phenotype is \code{y = sqrt(H2) G + sqrt(1 - H2) e}. \code{s} is
#' the targeted true digenic fraction of the genetic value.
#'
#' @param KP a \code{\link{build_panel}} result built with \code{need_X =
#'   TRUE} (needs \code{X[[1]]}, \code{X[[2]]} and \code{valid}).
#' @param s targeted true digenic fraction of the genetic value, in
#'   \code{[0, 1]}.
#' @param H2 broad-sense heritability in \code{[0, 1]}; defaults to
#'   \code{\link{PWF_H2_DEFAULT}}.
#' @return a length-n numeric phenotype vector.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
#' y <- gen_d2_spikein(KP, s = 0.5, H2 = 0.6)
#' length(y)
#' @export
gen_d2_spikein <- function(KP, s, H2 = PWF_H2_DEFAULT) {
  gA  <- sim_comp(KP$X[[1]], KP$valid)
  gD2 <- sim_comp(KP$X[[2]], KP$valid)
  rD2 <- stats::residuals(stats::lm(gD2 ~ gA)); rD2 <- as.numeric(scale(rD2))
  G <- sqrt(1 - s) * gA + sqrt(s) * rD2
  n <- length(G); e <- stats::rnorm(n)
  as.numeric(sqrt(H2) * as.numeric(scale(G)) + sqrt(1 - H2) * e)
}

#' Population-structure known-truth phenotype
#'
#' Draws an additive genetic value \code{g_A} (\code{\link{sim_comp}}) and a
#' structure component \code{g_struct}: standardized genomic PC1
#' (\code{structure_type = "PC"}) or a standardized family/group block effect
#' (\code{structure_type = "family"}, one fresh Gaussian effect per group
#' level, which leaks into the dominance kernel). The genetic value is
#' \code{G = sqrt(1 - lambda) g_A + sqrt(lambda) g_struct}, and the phenotype
#' is \code{y = sqrt(H2) G + sqrt(1 - H2) e}.
#'
#' @param KP a \code{\link{build_panel}} result built with \code{need_X =
#'   TRUE} (needs \code{X[[1]]} and \code{valid}).
#' @param PC an n x npc matrix of genomic PCs (e.g. from
#'   \code{\link{panel_pc}}); only column 1 is used, and only when
#'   \code{structure_type = "PC"}.
#' @param group a vector of family/group labels, length n; only used when
#'   \code{structure_type = "family"}.
#' @param lambda structure fraction of the genetic value, in \code{[0, 1]}.
#' @param structure_type \code{"PC"} (genomic-PC structure) or \code{"family"}
#'   (family/group block structure).
#' @param H2 broad-sense heritability in \code{[0, 1]}; defaults to
#'   \code{\link{PWF_H2_DEFAULT}}.
#' @return a length-n numeric phenotype vector.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
#' PC <- panel_pc(KP, 5L)
#' grp <- rep(c("a", "b"), length.out = 40)
#' y <- gen_structure(KP, PC, grp, lambda = 0.5, structure_type = "PC", H2 = 0.6)
#' length(y)
#' @export
gen_structure <- function(KP, PC, group, lambda, structure_type, H2 = PWF_H2_DEFAULT) {
  gA <- sim_comp(KP$X[[1]], KP$valid); n <- length(gA)
  if (structure_type == "PC") {
    gS <- as.numeric(scale(PC[, 1]))
  } else {                                  # family / group block effect (leaks into D kernel)
    lv <- unique(group); eff <- stats::rnorm(length(lv)); names(eff) <- lv
    gS <- as.numeric(scale(eff[group]))
  }
  G <- sqrt(1 - lambda) * gA + sqrt(lambda) * gS
  e <- stats::rnorm(n)
  as.numeric(sqrt(H2) * as.numeric(scale(G)) + sqrt(1 - H2) * e)
}

#' Simulate a pure-degree genotypic value (falsification-check simulator)
#'
#' Assigns one fresh Gaussian effect per causal marker, in marker order, and
#' sums each causal marker's degree-\code{degree} orthogonal contrast scaled
#' by that effect. Unlike \code{\link{sim_comp}}/\code{\link{gen_d2_spikein}},
#' this simulator works directly off the raw dosage matrix \code{M} (via
#' \code{\link{validate_dosage_matrix}}) rather than a pre-built panel, and
#' produces a genotypic value of a single, exactly pure degree — used to check
#' that a degree-\code{d} planted signal is recovered as degree \code{d} and
#' not smeared across other kernels.
#'
#' @param M individuals x markers dosage matrix, values in \code{0..ploidy}.
#' @param causal integer indices (into the columns of \code{M}) of the causal
#'   markers.
#' @param degree contrast degree to simulate: 1 = additive, 2 = digenic, 3 =
#'   trigenic, 4 = quadrigenic, etc.; must be in \code{1..ploidy}.
#' @param ploidy even ploidy (default 4).
#' @param impute if \code{TRUE}, mean-impute missing dosage calls instead of
#'   erroring; see \code{\link{validate_dosage_matrix}}.
#' @param min_maf minimum marker minor-allele frequency to retain.
#' @param min_call_rate minimum marker call rate to retain.
#' @return a length-\code{nrow(M)} numeric genotypic-value vector.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' gv <- simulate_pure_degree(M, causal = 1:10, degree = 2, ploidy = 4)
#' length(gv)
#' @export
simulate_pure_degree <- function(M, causal, degree, ploidy = 4L, impute = FALSE,
                                 min_maf = 0, min_call_rate = 0) {
  if (!is.numeric(degree) || length(degree) != 1L || degree < 1 || degree > ploidy) {
    stop("degree must be in 1..ploidy")
  }
  chk <- validate_dosage_matrix(M, ploidy, impute = impute, min_maf = min_maf, min_call_rate = min_call_rate)
  G <- chk$G; ploidy <- chk$ploidy
  retained <- which(chk$keep)
  causal <- match(as.integer(causal), retained, nomatch = 0L)
  causal <- causal[causal > 0L]
  if (!length(causal)) stop("no causal markers remain after filtering")
  gv <- numeric(nrow(G))
  for (j in causal) {
    g <- G[, j]; p <- mean(g) / ploidy
    gv <- gv + stats::rnorm(1) * ortho_contrasts(p, ploidy)[g + 1L, degree + 1L]
  }
  gv
}

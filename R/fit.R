## Internal: proper posterior credible intervals on variance fractions, from
## BGLR's saved per-component varU samples (<pfx>ETA_<name>_varU.dat). Aligns
## the kernels' samples per iteration, drops burn-in, then takes quantiles of
## the fractions. Ported from experiments/simulations/poly_ortho/lib_ortho.R
## (frac_ci_samples), renamed to avoid colliding with that file's global.
frac_ci_samples_polydx <- function(pfx, etanames, burnIn, thin) {
  S <- lapply(etanames, function(nm) scan(paste0(pfx, "ETA_", nm, "_varU.dat"), quiet = TRUE))
  L <- min(vapply(S, length, integer(1))); S <- lapply(S, function(v) v[seq_len(L)])
  drop <- floor(burnIn / thin)
  if (drop >= L) {
    stop(sprintf("burnIn (%g) discards all %d saved samples; use nIter > burnIn", burnIn, L),
         call. = FALSE)
  }
  keep <- (drop + 1):L
  Mn <- do.call(cbind, lapply(S, function(v) v[keep]))
  Fr <- Mn / rowSums(Mn)
  out <- apply(Fr, 2, function(cc) c(mean = mean(cc),
                                     lo = stats::quantile(cc, .025, names = FALSE),
                                     hi = stats::quantile(cc, .975, names = FALSE)))
  colnames(out) <- etanames; out
}

#' Fit a BGLR RKHS model and return posterior variance-fraction intervals
#'
#' Fits \code{y} on the supplied kernels/fixed effects (\code{ETA}) with
#' \code{\link[BGLR:BGLR]{BGLR::BGLR}}, then reads BGLR's per-iteration
#' variance-component samples off disk, aligns them, drops burn-in and
#' converts each kernel's variance to a posterior *fraction* of the total
#' fitted variance at every retained iteration, returning the mean and a
#' 95\% credible interval per fraction.
#'
#' @param y numeric response vector, length n (\code{NA} allowed for
#'   held-out individuals, per BGLR conventions).
#' @param ETA a BGLR \code{ETA} list, one entry per model component (e.g.
#'   \code{list(A = list(K = KP$A, model = "RKHS"), D = list(K = KP$D, model
#'   = "RKHS"))}).
#' @param etanames character vector naming the components whose variance
#'   fractions are wanted; must match \code{names(ETA)} for the components of
#'   interest, in the order the fraction matrix should report them.
#' @param nIter total number of MCMC iterations.
#' @param burnIn number of initial iterations to discard.
#' @param thin thinning interval BGLR was run with.
#' @return a list with \code{fm} (the fitted \code{BGLR} object) and
#'   \code{ci} (a \code{3 x length(etanames)} matrix with rows \code{mean},
#'   \code{lo}, \code{hi} and columns \code{etanames}, giving the posterior
#'   mean and 95\% credible interval of each component's variance
#'   fraction).
#' @examples
#' if (requireNamespace("BGLR", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30L; m <- 60L
#'   M <- matrix(rbinom(n * m, 4, 0.4), n, m)
#'   rownames(M) <- paste0("ind", seq_len(n))
#'   KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)
#'   y <- as.numeric(scale(KP$X[[1]] %*% rnorm(m))) + rnorm(n, 0, 0.8)
#'   rr <- run_rkhs(y, list(A = list(K = KP$A, model = "RKHS"),
#'                          D = list(K = KP$D, model = "RKHS")),
#'                  c("A", "D"), nIter = 300, burnIn = 100, thin = 5)
#'   rr$ci
#' }
#' @export
run_rkhs <- function(y, ETA, etanames,
                     nIter = pwf_opt("PWF_FRAC_NITER", 12000L),
                     burnIn = pwf_opt("PWF_FRAC_BURN", 3000L),
                     thin = pwf_opt("PWF_FRAC_THIN", 5L)) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("package 'BGLR' is required for variance-component fitting; install.packages(\"BGLR\")", call. = FALSE)
  }
  d <- tempfile("bg"); dir.create(d); pfx <- file.path(d, "x_")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fm <- BGLR::BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,
                   verbose = FALSE, saveAt = pfx)
  list(fm = fm, ci = frac_ci_samples_polydx(pfx, etanames, burnIn, thin))
}

#' Standard A+D and orthogonal-degree dominance-share fit metrics
#'
#' Fits two RKHS models on the panel's kernels (the standard additive/
#' dominance pair \code{KP$A}/\code{KP$D}, and the orthogonal degree kernels
#' \code{KP$O}) via \code{\link{run_rkhs}}, and packages the manuscript's
#' fraction metrics: the standard dominance share \code{s_hat_AD} with its
#' credible interval, and the orthogonal digenic (D2) fraction of total
#' fitted genetic variance with its credible interval and whether D2 is the
#' largest non-additive degree.
#'
#' @param y numeric response vector, length n.
#' @param KP a \code{\link{build_panel}} result (needs \code{A}, \code{D},
#'   \code{O}, \code{ndeg}, \code{etanames}).
#' @param nIter total number of MCMC iterations for both fits.
#' @param burnIn number of initial iterations to discard for both fits.
#' @param thin thinning interval for both fits.
#' @return a one-row \code{data.frame} with columns \code{s_hat_AD},
#'   \code{s_hat_AD_lo}, \code{s_hat_AD_hi} (standard dominance share and its
#'   95\% CrI), \code{s_hat_D2_total} (orthogonal digenic fraction of total
#'   fitted genetic variance), \code{P_D2} (digenic fraction of \code{A + }
#'   non-additive degrees; equals \code{s_hat_D2_total}), \code{L_D2}
#'   (\code{D2} fraction minus the largest higher-degree fraction),
#'   \code{D2_CrI_low}, \code{D2_CrI_high}, \code{D2_CrI_width} and
#'   \code{D2_largest} (logical: is D2 the largest non-additive degree).
#'   When the panel has no degree above 2 (\code{ndeg == 2}, i.e. ploidy 2)
#'   there is no competing degree, so \code{L_D2} and \code{D2_largest} are
#'   \code{NA} rather than trivially favourable.
#' @examples
#' if (requireNamespace("BGLR", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30L; m <- 60L
#'   M <- matrix(rbinom(n * m, 4, 0.4), n, m)
#'   rownames(M) <- paste0("ind", seq_len(n))
#'   KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)
#'   y <- as.numeric(scale(KP$X[[1]] %*% rnorm(m))) + rnorm(n, 0, 0.8)
#'   fr <- fit_fraction(y, KP, nIter = 300, burnIn = 100, thin = 5)
#'   fr$s_hat_AD
#' }
#' @export
fit_fraction <- function(y, KP,
                         nIter = pwf_opt("PWF_FRAC_NITER", 12000L),
                         burnIn = pwf_opt("PWF_FRAC_BURN", 3000L),
                         thin = pwf_opt("PWF_FRAC_THIN", 5L)) {
  en <- KP$etanames                            # c("A","D2",...,"Dndeg")
  rs <- run_rkhs(y, list(A = list(K = KP$A, model = "RKHS"), D = list(K = KP$D, model = "RKHS")),
                 c("A", "D"), nIter = nIter, burnIn = burnIn, thin = thin)
  eta_o <- list(); for (k in seq_len(KP$ndeg)) eta_o[[en[k]]] <- list(K = KP$O[[k]], model = "RKHS")
  ro <- run_rkhs(y, eta_o, en, nIter = nIter, burnIn = burnIn, thin = thin)
  s <- rs$ci; f <- ro$ci; gg <- function(m, r, c) as.numeric(m[r, c])
  D2 <- gg(f, "mean", "D2")
  afrac <- gg(f, "mean", "A")
  higher <- if (KP$ndeg >= 3) en[3:KP$ndeg] else character(0)
  vhigh <- if (length(higher)) vapply(higher, function(c) gg(f, "mean", c), numeric(1)) else numeric(0)
  nonadd <- D2 + sum(vhigh)
  # With no higher degrees (ndeg == 2, i.e. ploidy 2) there is nothing for D2
  # to be larger *than*. Substituting 0 made D2_largest unconditionally TRUE
  # -- BGLR variance fractions are strictly positive -- so the scorecard's
  # recovery gate passed without evidence, and L_D2 collapsed onto P_D2,
  # costing resembles_spike its independent check. NA says "not evaluated".
  max_high <- if (length(vhigh)) max(vhigh) else NA_real_
  data.frame(
    s_hat_AD       = gg(s, "mean", "D"),
    s_hat_AD_lo    = gg(s, "lo", "D"),
    s_hat_AD_hi    = gg(s, "hi", "D"),
    s_hat_D2_total = if (afrac + nonadd > 0) D2 else NA_real_,
    # P_D2 = V_D2 / (V_A + V_D2 + ... + V_Dm): digenic fraction of TOTAL fitted
    # genetic variance (additive INCLUDED in the denominator). Equals
    # s_hat_D2_total since the degree fractions sum to one.
    P_D2           = if (afrac + nonadd > 0) D2 / (afrac + nonadd) else NA_real_,
    L_D2           = D2 - max_high,               # NA when no higher degree exists
    D2_CrI_low     = gg(f, "lo", "D2"),
    D2_CrI_high    = gg(f, "hi", "D2"),
    D2_CrI_width   = gg(f, "hi", "D2") - gg(f, "lo", "D2"),
    D2_largest     = if (is.na(max_high)) NA else isTRUE(D2 > max_high),
    stringsAsFactors = FALSE
  )
}

#' Fast point-estimate dominance share, optionally controlled for structure
#'
#' Fits a single BGLR RKHS/FIXED model and returns the point-estimate
#' dominance share \eqn{V_D / (V_A + V_D)} directly from the fitted variance
#' components (no credible interval; use \code{\link{fit_fraction}} for
#' that). Optionally adds a genomic-PC fixed-effect block or a family
#' relationship kernel to the model to control for population structure or
#' relatedness before reading off the dominance share.
#'
#' @param y numeric response vector, length n.
#' @param KP a \code{\link{build_panel}} result (needs \code{A}, \code{D}).
#' @param PC n x npc matrix of genomic PCs (e.g. from \code{\link{panel_pc}});
#'   only used when \code{what = "PC10"}.
#' @param Fk n x n family/relatedness kernel (e.g. from
#'   \code{\link{family_kernel}}); only used when \code{what = "family"}.
#' @param what one of \code{"raw"} (A + D only), \code{"PC10"} (A + D plus
#'   the first \code{min(10, ncol(PC))} PCs as fixed effects) or
#'   \code{"family"} (A + D plus \code{Fk} as an RKHS block).
#' @param nIter total number of MCMC iterations.
#' @param burnIn number of initial iterations to discard.
#' @return a numeric scalar, the posterior point estimate of
#'   \eqn{V_D / (V_A + V_D)}.
#' @examples
#' if (requireNamespace("BGLR", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30L; m <- 60L
#'   M <- matrix(rbinom(n * m, 4, 0.4), n, m)
#'   rownames(M) <- paste0("ind", seq_len(n))
#'   KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)
#'   y <- as.numeric(scale(KP$X[[1]] %*% rnorm(m))) + rnorm(n, 0, 0.8)
#'   share_pt(y, KP, what = "raw", nIter = 300, burnIn = 100)
#' }
#' @export
share_pt <- function(y, KP, PC = NULL, Fk = NULL, what = "raw",
                     nIter = pwf_opt("PWF_PT_NITER", 3000L),
                     burnIn = pwf_opt("PWF_PT_BURN", 750L)) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("package 'BGLR' is required for variance-component fitting; install.packages(\"BGLR\")", call. = FALSE)
  }
  if (burnIn >= nIter) {
    stop(sprintf("burnIn (%g) leaves no post-burn-in draws of nIter (%g)", burnIn, nIter),
         call. = FALSE)
  }
  if (!what %in% c("raw", "PC10", "family")) {
    stop("what must be one of \"raw\", \"PC10\" or \"family\"", call. = FALSE)
  }
  # Falling back to the raw share when the requested control's input is missing
  # would report S = raw - raw = 0, i.e. "structure-insensitive", for a control
  # that was never actually fitted.
  if (what == "PC10" && is.null(PC)) stop("what = \"PC10\" requires PC", call. = FALSE)
  if (what == "family" && is.null(Fk)) stop("what = \"family\" requires Fk", call. = FALSE)
  ETA <- list()
  if (what == "PC10") ETA$PC <- list(X = PC[, seq_len(min(10L, ncol(PC))), drop = FALSE], model = "FIXED")
  if (what == "family") ETA$Fam <- list(K = Fk, model = "RKHS")
  ETA$A <- list(K = KP$A, model = "RKHS"); ETA$D <- list(K = KP$D, model = "RKHS")
  d <- tempfile("pt"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fm <- BGLR::BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE,
                   saveAt = file.path(d, paste0(what, "_")))
  as.numeric(fm$ETA[["D"]]$varU / (fm$ETA[["A"]]$varU + fm$ETA[["D"]]$varU))
}

#' Partition fitted variance across an arbitrary set of kernels
#'
#' Fits a single BGLR RKHS model with one block per kernel in \code{kernels}
#' and returns each block's posterior variance as a fraction of their sum.
#'
#' @param y numeric response vector, length n.
#' @param kernels a named list of n x n relationship kernels, one RKHS block
#'   per kernel.
#' @param nIter total number of MCMC iterations.
#' @param burnIn number of initial iterations to discard.
#' @return a named numeric vector (names = \code{names(kernels)}) of
#'   variance fractions summing to 1.
#' @examples
#' if (requireNamespace("BGLR", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30L; m <- 60L
#'   M <- matrix(rbinom(n * m, 4, 0.4), n, m)
#'   rownames(M) <- paste0("ind", seq_len(n))
#'   KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)
#'   y <- as.numeric(scale(KP$X[[1]] %*% rnorm(m))) + rnorm(n, 0, 0.8)
#'   partition_variance(y, list(A = KP$A, D = KP$D), nIter = 300, burnIn = 100)
#' }
#' @export
partition_variance <- function(y, kernels, nIter = 6000L, burnIn = 1500L) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("package 'BGLR' is required for variance-component fitting; install.packages(\"BGLR\")", call. = FALSE)
  }
  if (burnIn >= nIter) {
    stop(sprintf("burnIn (%g) leaves no post-burn-in draws of nIter (%g)", burnIn, nIter),
         call. = FALSE)
  }
  ETA <- lapply(kernels, function(K) list(K = K, model = "RKHS"))
  d <- tempfile("po"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  fm <- BGLR::BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE,
                   saveAt = file.path(d, "x_"))
  v <- vapply(names(kernels), function(nm) fm$ETA[[nm]]$varU, numeric(1))
  v / sum(v)
}

#' Top genomic PCs from a panel's additive kernel
#'
#' Eigendecomposes the panel's additive kernel \code{KP$A} and returns the
#' leading \code{npc} eigenvectors as a genomic-PC matrix, for use as a
#' structure-control fixed effect (see \code{\link{share_pt}}).
#'
#' @param KP a \code{\link{build_panel}} result (needs \code{A}).
#' @param npc number of leading PCs to return (capped at \code{ncol(KP$A)}).
#' @return an n x \code{npc} matrix of eigenvectors, row names matching
#'   \code{rownames(KP$A)}.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(30 * 50, 4, 0.4), 30, 50)
#' rownames(M) <- paste0("ind", seq_len(30))
#' KP <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = FALSE)
#' PC <- panel_pc(KP, 5L)
#' dim(PC)
#' @export
panel_pc <- function(KP, npc = 10L) {
  ev <- eigen(KP$A, symmetric = TRUE)
  PC <- ev$vectors[, seq_len(min(npc, ncol(ev$vectors))), drop = FALSE]
  rownames(PC) <- rownames(KP$A); PC
}

#' Family/group incidence kernel
#'
#' Builds an n x n relationship kernel from a family/group label vector: the
#' cross-product of the group-indicator design matrix, scaled to mean
#' diagonal 1 with a small ridge for numerical stability. Individuals in the
#' same group have kernel entries near 1; individuals in different groups
#' have entries near 0.
#'
#' @param group a vector of family/group labels, length n, or \code{NULL};
#'   \code{NA} labels are pooled into a single \code{"NA_grp"} level.
#' @param ids row/column names to assign to the kernel (typically
#'   individual IDs matching \code{group}'s order).
#' @return an n x n kernel matrix with \code{dimnames} set to
#'   \code{list(ids, ids)}, or \code{NULL} when \code{group} is \code{NULL}.
#' @examples
#' grp <- rep(c("f1", "f2"), each = 3)
#' ids <- paste0("ind", seq_along(grp))
#' Fk <- family_kernel(grp, ids)
#' dim(Fk)
#' is.null(family_kernel(NULL, ids))
#' @export
family_kernel <- function(group, ids) {
  if (is.null(group)) return(NULL)
  # model.matrix() drops NA rows, which would silently shrink the kernel below
  # length(ids); label them instead, matching panel_from_matrix().
  group <- as.character(group); group[is.na(group)] <- "NA_grp"
  kernel_from(stats::model.matrix(~ 0 + factor(group)), ids)
}

#' Unnamed quantile
#'
#' Thin wrapper around \code{stats::quantile} that drops the percentage-style
#' name and ignores \code{NA}s, for use in summaries/scorecards where a bare
#' numeric is wanted.
#'
#' @param x numeric vector.
#' @param p quantile probability between 0 and 1.
#' @return an unnamed numeric scalar.
#' @examples
#' q_named(1:100, 0.5)
#' @export
q_named <- function(x, p) as.numeric(stats::quantile(x, p, names = FALSE, na.rm = TRUE))

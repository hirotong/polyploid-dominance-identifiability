# =====================================================================
# lib_reml.R -- restricted maximum likelihood for kernel variance-component
# models with an intercept only, in base R. Used by estimator_sensitivity.R as a
# prior-free counterpart to the BGLR fits.
#
#   y = mu + sum_k g_k + e,  g_k ~ N(0, s2_k K_k),  e ~ N(0, s2_e I)
#
# The intercept is removed by projecting onto contrasts orthogonal to 1, so every
# likelihood below is a restricted likelihood in n - 1 dimensions.
#
#   reml_contrasts(n)                 n x (n-1) orthonormal basis orthogonal to 1
#   reml_profile_share(Y, KA, KD, s)  two-kernel model: profile restricted
#                                     log-likelihood over the share s = s2_D /
#                                     (s2_A + s2_D), for every column of Y at once
#   reml_share_summary(ll, s)         estimate and likelihood-ratio interval
#   reml_multi(y, Klist)              general k-kernel REML by L-BFGS-B on log variances
# =====================================================================

reml_contrasts <- function(n) {
  # Helmert-type orthonormal contrasts: columns orthogonal to the 1 vector
  Q <- qr.Q(qr(cbind(1, diag(n)[, -1, drop = FALSE])))
  Q[, -1, drop = FALSE]
}

# Profile restricted log-likelihood (up to a constant) of a one-kernel model with
# eigenvalues d of the projected kernel and rotated data eta (a matrix: one column
# per phenotype). For each column, maximises over delta = s2_e / s2_g.
.reml_profile_1k <- function(d, eta) {
  m <- length(d)
  apply(eta, 2, function(e2) {
    f <- function(ld) {
      w <- d + exp(ld)
      -0.5 * (m * log(sum(e2 / w) / m) + sum(log(w)) + m)
    }
    grid <- seq(-12, 12, length.out = 49)                 # coarse scan, then refine
    v <- vapply(grid, f, numeric(1)); i <- which.max(v)
    lo <- grid[max(1, i - 1)]; hi <- grid[min(length(grid), i + 1)]
    o <- optimize(f, c(lo, hi), maximum = TRUE)
    max(o$objective, v[i])
  }, simplify = TRUE)
}

# Y: n x q matrix of phenotypes without missing values (column 1 the real trait,
# other columns e.g. permutations). Returns a length(s) x q matrix of profile
# restricted log-likelihoods.
reml_profile_share <- function(Y, KA, KD, s = seq(0, 1, by = 0.02)) {
  Y <- as.matrix(Y); n <- nrow(Y)
  A <- reml_contrasts(n)
  Ys <- crossprod(A, Y)
  KAs <- crossprod(A, KA %*% A); KDs <- crossprod(A, KD %*% A)
  ll <- matrix(NA_real_, length(s), ncol(Y))
  for (i in seq_along(s)) {
    ev <- eigen((1 - s[i]) * KAs + s[i] * KDs, symmetric = TRUE)
    d <- pmax(ev$values, 0)
    eta2 <- crossprod(ev$vectors, Ys)^2
    ll[i, ] <- .reml_profile_1k(d, eta2)
  }
  ll
}

# Estimate (grid maximiser) and 95% likelihood-ratio interval for the share from
# one column of profile log-likelihoods. The chi-square(1) cut-off is approximate
# at the boundaries s = 0 and s = 1.
reml_share_summary <- function(ll, s) {
  i <- which.max(ll); inside <- which(2 * (ll[i] - ll) <= qchisq(0.95, 1))
  c(share = s[i], lo = s[min(inside)], hi = s[max(inside)],
    lr_s0 = 2 * (ll[i] - ll[which.min(abs(s - 0))]),   # LR against additive only
    lr_s1 = 2 * (ll[i] - ll[which.min(abs(s - 1))]))   # LR against D kernel only
}

# General REML for k kernels. Returns variance components and genetic fractions.
reml_multi <- function(y, Klist, starts = 2L, maxit = 200L) {
  n <- length(y); A <- reml_contrasts(n); ys <- as.numeric(crossprod(A, y)); m <- n - 1L
  Ks <- lapply(Klist, function(K) crossprod(A, K %*% A)); k <- length(Ks); vy <- stats::var(y)
  nll <- function(th) {                                   # th = log variances (k genetic, then residual)
    v <- exp(th); V <- diag(v[k + 1], m); for (j in seq_len(k)) V <- V + v[j] * Ks[[j]]
    R <- tryCatch(chol(V), error = function(e) NULL); if (is.null(R)) return(1e12)
    z <- backsolve(R, ys, transpose = TRUE)
    sum(log(diag(R))) + 0.5 * sum(z^2)
  }
  gr <- function(th) {
    v <- exp(th); V <- diag(v[k + 1], m); for (j in seq_len(k)) V <- V + v[j] * Ks[[j]]
    Vi <- chol2inv(chol(V)); a <- as.numeric(Vi %*% ys)
    g <- numeric(k + 1)
    for (j in seq_len(k)) g[j] <- 0.5 * (sum(Vi * Ks[[j]]) - sum(a * (Ks[[j]] %*% a))) * v[j]
    g[k + 1] <- 0.5 * (sum(diag(Vi)) - sum(a^2)) * v[k + 1]
    g
  }
  lower <- rep(log(vy) - 14, k + 1); upper <- rep(log(vy) + 3, k + 1)
  inits <- list(log(c(rep(0.5 * vy / k, k), 0.5 * vy)),
                log(c(0.4 * vy, rep(0.1 * vy / max(1, k - 1), k - 1), 0.5 * vy)))[seq_len(min(starts, 2L))]
  best <- NULL
  for (th0 in inits) {
    o <- tryCatch(optim(th0, nll, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                        control = list(maxit = maxit)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(NULL)
  v <- exp(best$par); vg <- v[seq_len(k)]; names(vg) <- names(Klist)
  list(var = vg, var_e = v[k + 1], frac = vg / sum(vg), nll = best$value, convergence = best$convergence)
}

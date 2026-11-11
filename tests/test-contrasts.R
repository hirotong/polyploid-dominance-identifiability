library(polyDx)

## ortho_contrasts: column 1 is the constant (degree-0) contrast
Q <- ortho_contrasts(0.3, 4L)
stopifnot(is.matrix(Q), identical(dim(Q), c(5L, 5L)))
stopifnot(all(abs(Q[, 1] - Q[1, 1]) < 1e-10))

## orthonormal under the binomial weights w_d at p
p <- 0.3; ploidy <- 4L; d <- 0:ploidy
w <- choose(ploidy, d) * p^d * (1 - p)^(ploidy - d); w <- w / sum(w)
G <- t(Q) %*% diag(w) %*% Q
stopifnot(max(abs(G - diag(ploidy + 1L))) < 1e-8)

## degree-1 contrast is monotone in dosage (it is the additive contrast)
stopifnot(all(diff(Q[, 2]) > 0) || all(diff(Q[, 2]) < 0))

## ploidy 6 also works
Q6 <- ortho_contrasts(0.45, 6L)
stopifnot(identical(dim(Q6), c(7L, 7L)))

## impute_mean_dosage replaces NA with the marker mean and leaves others alone
M <- matrix(c(0, 2, 4, NA, 1, 3), nrow = 3)
Mi <- impute_mean_dosage(M)
stopifnot(!anyNA(Mi))
stopifnot(abs(Mi[1, 2] - mean(c(1, 3))) < 1e-12)
stopifnot(identical(Mi[, 1], as.numeric(M[, 1])))

## an all-NA marker is an error, not a silent NaN
stopifnot(inherits(try(impute_mean_dosage(matrix(NA_real_, 3, 1)), silent = TRUE), "try-error"))

## validate_dosage_matrix: rejects out-of-range dosage and bad ploidy
set.seed(11)
Mok <- matrix(rbinom(60, 4, 0.4), nrow = 10)
stopifnot(inherits(try(validate_dosage_matrix(Mok, ploidy = 1), silent = TRUE), "try-error"))
stopifnot(inherits(try(validate_dosage_matrix(Mok + 5, ploidy = 4), silent = TRUE), "try-error"))

## refuses NA unless impute = TRUE, accepts it when imputing
Mna <- Mok; Mna[1, 1] <- NA
stopifnot(inherits(try(validate_dosage_matrix(Mna, ploidy = 4), silent = TRUE), "try-error"))
v <- validate_dosage_matrix(Mna, ploidy = 4, impute = TRUE)
stopifnot(!anyNA(v$M), all(v$p > 0), all(v$p < 1), identical(v$ploidy, 4L))

## monomorphic markers are dropped; all-monomorphic is an error
Mmono <- cbind(Mok, rep(0, 10))
v2 <- validate_dosage_matrix(Mmono, ploidy = 4)
stopifnot(ncol(v2$M) == ncol(Mok), isFALSE(v2$keep[ncol(Mmono)]))
stopifnot(inherits(try(validate_dosage_matrix(matrix(0, 10, 3), ploidy = 4), silent = TRUE), "try-error"))

## pwf_opt honours the environment and the default's type
old <- Sys.getenv("POLYDX_TEST_OPT", unset = NA)
Sys.setenv(POLYDX_TEST_OPT = "0.25")
stopifnot(abs(pwf_opt("POLYDX_TEST_OPT", 0.5) - 0.25) < 1e-12)
Sys.unsetenv("POLYDX_TEST_OPT")
stopifnot(abs(pwf_opt("POLYDX_TEST_OPT", 0.5) - 0.5) < 1e-12)
stopifnot(identical(pwf_opt("POLYDX_TEST_OPT", 7L), 7L))
if (!is.na(old)) do.call(Sys.setenv, setNames(list(old), "POLYDX_TEST_OPT"))

stopifnot(abs(PWF_H2_DEFAULT - 0.6) < 1e-12)
cat("test-contrasts.R OK\n")

## --- regression: ortho_contrasts validates its (0, 1) domain ---------------
## p at 0 or 1 makes the binomial weights degenerate; Gram-Schmidt used to fail
## with "missing value where TRUE/FALSE needed" instead of a clear message.
for (bad in list(0, 1, -0.1, 1.1, NA_real_, c(0.2, 0.3))) {
  stopifnot(inherits(try(ortho_contrasts(bad, 4L), silent = TRUE), "try-error"))
}
stopifnot(inherits(try(ortho_contrasts(0.3, 0L), silent = TRUE), "try-error"))
stopifnot(is.matrix(ortho_contrasts(0.3, 4L)))

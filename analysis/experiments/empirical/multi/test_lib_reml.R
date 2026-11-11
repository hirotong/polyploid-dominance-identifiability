#!/usr/bin/env Rscript
# Known-truth check of lib_reml.R (base R only; about 1 minute).
#   Rscript experiments/empirical/multi/test_lib_reml.R
# Simulates phenotypes from two well-separated kernels with known share s and checks
# that the profile REML share recovers it on average, that a permuted phenotype gives
# an interval covering the whole range, and that the general optimiser agrees with
# the profile for two kernels.
source("experiments/empirical/multi/lib_reml.R"); set.seed(7)
n <- 500; mk <- function(m) { Z <- scale(matrix(rnorm(n * m), n, m)); K <- tcrossprod(Z) / m; K / mean(diag(K)) }
KA <- mk(300); KD <- mk(300)
sim <- function(s, h2 = 0.6) { L <- function(K) t(chol(K + diag(1e-8, n))) %*% rnorm(n)
  g <- sqrt(1 - s) * L(KA) + sqrt(s) * L(KD); as.numeric(sqrt(h2) * g + sqrt(1 - h2) * rnorm(n)) + 3 }
sg <- seq(0, 1, by = 0.02)
for (s in c(0, 0.3, 0.7)) {
  est <- replicate(8, { y <- sim(s); reml_share_summary(reml_profile_share(cbind(y), KA, KD, sg)[, 1], sg)[["share"]] })
  cat(sprintf("true s=%.1f  mean REML share=%.3f  sd=%.3f\n", s, mean(est), sd(est)))
  stopifnot(abs(mean(est) - s) < 0.1)
}
y <- sim(0.3); ll <- reml_profile_share(cbind(y, sample(y)), KA, KD, sg)
r <- reml_share_summary(ll[, 1], sg); p <- reml_share_summary(ll[, 2], sg)
m2 <- reml_multi(y, list(A = KA, D = KD))
cat(sprintf("profile share=%.2f [%.2f, %.2f]; general optimiser=%.3f; permuted interval=[%.2f, %.2f]\n",
            r[["share"]], r[["lo"]], r[["hi"]], m2$frac[["D"]], p[["lo"]], p[["hi"]]))
stopifnot(abs(r[["share"]] - m2$frac[["D"]]) < 0.03, p[["lo"]] == 0, p[["hi"]] == 1)
cat("lib_reml checks passed\n")

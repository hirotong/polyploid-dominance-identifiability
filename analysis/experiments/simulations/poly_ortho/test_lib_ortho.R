# experiments/simulations/poly_ortho/test_lib_ortho.R
source("experiments/simulations/poly_ortho/lib_ortho.R")
stopifnot(exists("ortho"), exists("build_kernels"), exists("design_diag"))

# (1) contrasts are orthonormal under the binomial weight measure at p=0.3
p <- 0.3; ploidy <- 4L; Q <- ortho(p, ploidy)
d <- 0:ploidy; w <- choose(ploidy, d) * p^d * (1 - p)^(ploidy - d); w <- w / sum(w)
G <- t(Q) %*% diag(w) %*% Q                      # should be identity
stopifnot(max(abs(G - diag(ncol(Q)))) < 1e-8)

# (2) kernels from a tiny random dosage matrix are symmetric and finite
set.seed(1); M <- matrix(sample(0:4, 30 * 50, replace = TRUE), 30, 50)
rownames(M) <- paste0("i", 1:30); K <- build_kernels(M)
for (nm in c("A","D","OA","OD2","OD3","OD4"))
  stopifnot(isSymmetric(unname(K[[nm]]), tol = 1e-8), all(is.finite(K[[nm]])))

# (3) design diagnostics return positive, sane numbers
dd <- design_diag(list(K$OA, K$OD2, K$OD3, K$OD4))
stopifnot(dd$cond > 0, dd$eff_rank >= 1, dd$eff_rank <= 4)
cat("lib_ortho smoke test: PASS\n")

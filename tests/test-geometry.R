library(polyDx)
set.seed(9)
n <- 60; m <- 250
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))
KP <- build_panel(M, parse_ploidy("4", m), need_X = FALSE)

g <- kernel_geometry(KP)
stopifnot(is.finite(g$AD_cosine), g$AD_cosine > 0, g$AD_cosine <= 1 + 1e-8)
stopifnot(identical(dim(g$cosine_matrix), c(4L, 4L)))
stopifnot(max(abs(diag(g$cosine_matrix) - 1)) < 1e-8)     # self-cosine is 1
stopifnot(max(abs(g$cosine_matrix - t(g$cosine_matrix))) < 1e-10)
stopifnot(g$ideal_rank == 4L)
stopifnot(g$eff_rank >= 1 - 1e-8, g$eff_rank <= 4 + 1e-8) # eff_rank in [1, ndeg]
stopifnot(abs(g$rank_ratio - g$eff_rank / 4) < 1e-12)
stopifnot(is.finite(g$condition_number), g$condition_number >= 1)
stopifnot(is.finite(g$max_abs_degree_cosine), g$max_abs_degree_cosine <= 1 + 1e-8)

## ideal_rank/eff_rank scale with the number of components passed in, not just
## with KP itself: a reduced-degree geometry over only the first two O
## components reports ideal_rank == 2 and eff_rank bounded accordingly
g2 <- kernel_geometry(list(A = KP$A, D = KP$D, O = KP$O[1:2]))
stopifnot(g2$ideal_rank == 2L)
stopifnot(g2$eff_rank >= 1 - 1e-8, g2$eff_rank <= 2 + 1e-8)

## --- check_kernels --------------------------------------------------------
ids <- rownames(M)
stopifnot(isTRUE(check_kernels(KP, ids)))
bad <- KP; rownames(bad$A) <- rev(ids)
stopifnot(inherits(try(check_kernels(bad, ids), silent = TRUE), "try-error"))
bad2 <- KP; bad2$D[1, 1] <- NA
stopifnot(inherits(try(check_kernels(bad2, ids), silent = TRUE), "try-error"))
bad3 <- KP; bad3$A[1, 2] <- bad3$A[1, 2] + 1
stopifnot(inherits(try(check_kernels(bad3, ids), silent = TRUE), "try-error"))
cat("test-geometry.R OK\n")

## --- regression: check_kernels must validate the ORTHOGONAL A too ---------
## c(list(A = , D = ), KP$O) produced names A D A D2 D3 D4, and Ks[["A"]]
## resolved to the standard A twice, so corrupting or reordering KP$O$A -- the
## kernel every degree-wise fit uses -- passed silently.
KPg <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = FALSE)
ids <- rownames(M)
stopifnot(isTRUE(check_kernels(KPg, ids)))
bad_na <- KPg; bad_na$O$A[1, 2] <- NA
stopifnot(inherits(try(check_kernels(bad_na, ids), silent = TRUE), "try-error"))
bad_ord <- KPg; n_g <- nrow(M); bad_ord$O$A <- bad_ord$O$A[n_g:1, n_g:1]
stopifnot(inherits(try(check_kernels(bad_ord, ids), silent = TRUE), "try-error"))
bad_sym <- KPg; bad_sym$O$A[1, 2] <- bad_sym$O$A[1, 2] + 1
stopifnot(inherits(try(check_kernels(bad_sym, ids), silent = TRUE), "try-error"))

## dcenter's O(n^2) rank-2 form must equal the explicit J K J
Kd <- KPg$A; nd <- nrow(Kd); J <- diag(nd) - 1 / nd
stopifnot(max(abs(polyDx:::dcenter(Kd) - J %*% Kd %*% J)) < 1e-10)

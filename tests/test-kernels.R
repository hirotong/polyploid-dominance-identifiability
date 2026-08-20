library(polyDx)
set.seed(7)
n <- 80L; m <- 300L
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))

## --- ortho_etanames -------------------------------------------------------
stopifnot(identical(ortho_etanames(4L), c("A", "D2", "D3", "D4")))
stopifnot(identical(ortho_etanames(6L), c("A", "D2", "D3", "D4", "D5", "D6")))

## --- orthogonal_kernels ---------------------------------------------------
KO <- orthogonal_kernels(M)
stopifnot(identical(names(KO), c("A", "D2", "D3", "D4")))
for (K in KO) {
  stopifnot(identical(dim(K), c(n, n)))
  stopifnot(max(abs(K - t(K))) < 1e-8)             # symmetric
  stopifnot(abs(mean(diag(K)) - 1) < 0.02)         # mean diagonal ~ 1 (+1e-6 ridge)
  stopifnot(identical(rownames(K), rownames(M)))
  stopifnot(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
}

## --- standard_kernels -----------------------------------------------------
KS <- standard_kernels(M)
stopifnot(identical(names(KS), c("A", "D")))
stopifnot(max(abs(KS$A - t(KS$A))) < 1e-8, abs(mean(diag(KS$D)) - 1) < 0.02)

## --- build_panel: fixed ploidy -------------------------------------------
pc <- parse_ploidy("4", m)
KP <- build_panel(M, pc, need_X = TRUE)
stopifnot(identical(names(KP$O), c("A", "D2", "D3", "D4")))
stopifnot(KP$ndeg == 4L, KP$valid > 0, length(KP$X) == 4L)
stopifnot(identical(dim(KP$X[[1]]), c(n, m)))
stopifnot(identical(dim(KP$A), c(n, n)), identical(rownames(KP$D), rownames(M)))
stopifnot(identical(KP$etanames, ortho_etanames(4L)))

## build_panel's orthogonal kernels agree with orthogonal_kernels on clean data
stopifnot(max(abs(KP$O$D2 - KO$D2)) < 1e-8)

## need_X = FALSE drops the contrast matrices but keeps the kernels
KPn <- build_panel(M, pc, need_X = FALSE)
stopifnot(is.null(KPn$X), max(abs(KPn$A - KP$A)) < 1e-12)

## --- build_contrasts12 matches build_panel's first two degrees ------------
cc <- build_contrasts12(M, pc)
stopifnot(cc$valid == KP$valid, length(cc$X) == 2L)
stopifnot(max(abs(cc$X[[1]] - KP$X[[1]])) < 1e-12)
stopifnot(max(abs(cc$X[[2]] - KP$X[[2]])) < 1e-12)

## --- build_panel: variable ploidy ----------------------------------------
tmp <- tempfile(fileext = ".rds")
saveRDS(rep(c(4L, 6L), length.out = m), tmp)
Mv <- M
Mv[, seq(2, m, by = 2)] <- matrix(rbinom(n * length(seq(2, m, by = 2)), 6, 0.4), n)
rownames(Mv) <- rownames(M)
pcv <- parse_ploidy(paste0("var:", tmp, ":4"), m)
KPv <- build_panel(Mv, pcv, need_X = TRUE)
stopifnot(KPv$ndeg == 4L, identical(names(KPv$O), c("A", "D2", "D3", "D4")))
stopifnot(identical(dim(KPv$A), c(n, n)))
unlink(tmp)
cat("test-kernels.R OK\n")

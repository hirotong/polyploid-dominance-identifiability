library(polyDx)
set.seed(7)
n <- 80L; m <- 300L
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))

## --- ortho_etanames -------------------------------------------------------
stopifnot(identical(ortho_etanames(4L), c("A", "D2", "D3", "D4")))
stopifnot(identical(ortho_etanames(6L), c("A", "D2", "D3", "D4", "D5", "D6")))

## --- orthogonal_kernels ---------------------------------------------------
KO <- orthogonal_kernels(M, ploidy = 4L)
stopifnot(identical(names(KO), c("A", "D2", "D3", "D4")))
for (K in KO) {
  stopifnot(identical(dim(K), c(n, n)))
  stopifnot(max(abs(K - t(K))) < 1e-8)             # symmetric
  stopifnot(abs(mean(diag(K)) - 1) < 0.02)         # mean diagonal ~ 1 (+1e-6 ridge)
  stopifnot(identical(rownames(K), rownames(M)))
  stopifnot(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
}

## --- standard_kernels -----------------------------------------------------
KS <- standard_kernels(M, ploidy = 4L)
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

## --- regression: ortho_etanames must not count backwards at ndeg == 1 ------
## `2:ndeg` yields c(2, 1) at ndeg == 1, which produced c("A", "D2", "D1").
stopifnot(identical(ortho_etanames(1L), "A"))
stopifnot(identical(ortho_etanames(2L), c("A", "D2")))
for (bad in list(0L, -1L, NA_integer_, c(2L, 3L), "4")) {
  stopifnot(inherits(try(ortho_etanames(bad), silent = TRUE), "try-error"))
}

## --- build_contrasts12 is build_panel's degree-1/2 contrasts ---------------
pc4 <- parse_ploidy("4", ncol(M))
KPx <- build_panel(M, pc4, need_X = TRUE)
cc <- build_contrasts12(M, pc4)
stopifnot(identical(cc$X, KPx$X[1:2]), identical(cc$valid, KPx$valid))

## --- fixed ploidy == variable ploidy with a constant ploidy vector ---------
pvar <- list(mode = "variable", mploidy = rep(4L, ncol(M)), maxdeg = 4L,
             ndeg = 4L, label = "var<=4")
KPv <- build_panel(M, pvar, need_X = TRUE)
stopifnot(isTRUE(all.equal(KPx[c("A", "D", "O", "X", "valid")],
                           KPv[c("A", "D", "O", "X", "valid")], tolerance = 1e-12)))

## --- regression: preconditions fail loudly, not silently -------------------
## NA dosage used to surface as "missing value where TRUE/FALSE needed" from
## deep inside the marker loop; a degree above a marker's own ploidy used to
## leave an all-zero contrast column that plants no signal at all.
Mna <- M; Mna[1, 1] <- NA
stopifnot(inherits(try(build_panel(Mna, pc4, need_X = FALSE), silent = TRUE), "try-error"))
pmix <- list(mode = "variable", mploidy = rep(c(2L, 4L), length.out = ncol(M)),
             maxdeg = 4L, ndeg = 4L, label = "v")
stopifnot(inherits(try(build_panel(M, pmix, need_X = FALSE), silent = TRUE), "try-error"))

## --- regression: ploidy must be stated, never silently assumed to be 4 -----
## A diploid matrix analysed as tetraploid halves every allele frequency and
## inflates AD_cosine; the old ploidy = 4L defaults made that the quiet path.
for (call in list(function() orthogonal_kernels(M), function() standard_kernels(M))) {
  stopifnot(inherits(try(call(), silent = TRUE), "try-error"))
}

## --- regression: a degenerate component must error, not return NaN --------
## Dosages of only 1 or 3 give constant heterozygosity, so D's cross-product
## is all zeros and K / mean(diag(K)) used to yield an all-NaN kernel that
## surfaced as AD_cosine = NaN and only died later inside BGLR.
set.seed(404)
Mhet <- matrix(sample(c(1, 3), 40 * 40, TRUE), 40, 40)
rownames(Mhet) <- paste0("ind", seq_len(40))
qhet <- filter_markers(Mhet, parse_ploidy("4", 40L))
stopifnot(inherits(try(build_panel(qhet$M, qhet$pcfg, need_X = FALSE), silent = TRUE), "try-error"))

## --- regression: an all-heterozygous marker is monomorphic in dosage ------
## p = 0.5 and MAF 0.5, but zero dosage variance: its contrasts are a constant
## NON-zero value, so keeping it added a rank-1 constant to every degree kernel
## and over-counted `valid`.
Mah <- cbind(M, 2)
colnames(Mah) <- NULL
vd <- validate_dosage_matrix(Mah, ploidy = 4L)
stopifnot(sum(vd$keep) == sum(validate_dosage_matrix(M, ploidy = 4L)$keep))
stopifnot(!vd$keep[ncol(Mah)])

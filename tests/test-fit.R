library(polyDx)

## --- BGLR-free parts ------------------------------------------------------
set.seed(13)
n <- 40L; m <- 150L
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))
KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)

PC <- panel_pc(KP, 5L)
stopifnot(identical(dim(PC), c(n, 5L)), identical(rownames(PC), rownames(M)))

grp <- rep(c("f1", "f2", "f3", "f4"), length.out = n)
Fk <- family_kernel(grp, rownames(M))
stopifnot(identical(dim(Fk), c(n, n)), abs(mean(diag(Fk)) - 1) < 0.02)
stopifnot(identical(rownames(Fk), rownames(M)))
stopifnot(Fk[1, 5] > Fk[1, 2])                  # same family > different family
stopifnot(is.null(family_kernel(NULL, rownames(M))))

stopifnot(abs(q_named(1:100, 0.5) - 50.5) < 1e-9)
stopifnot(is.null(names(q_named(1:100, 0.95))))

## a missing BGLR must produce a clear error, not a cryptic one
if (!requireNamespace("BGLR", quietly = TRUE)) {
  e <- try(partition_variance(rnorm(n), list(A = KP$A)), silent = TRUE)
  stopifnot(inherits(e, "try-error"), grepl("BGLR", conditionMessage(attr(e, "condition"))))
  cat("test-fit.R OK (BGLR absent; fitting tests skipped)\n")
} else {

  ## --- BGLR-backed parts (tiny chains: correctness of shape, not of mixing) --
  y <- as.numeric(scale(KP$X[[1]] %*% rnorm(m))) + rnorm(n, 0, 0.8)

  fr <- fit_fraction(y, KP, nIter = 400, burnIn = 100, thin = 5)
  stopifnot(is.data.frame(fr), nrow(fr) == 1L)
  need <- c("s_hat_AD", "s_hat_AD_lo", "s_hat_AD_hi", "s_hat_D2_total", "P_D2",
            "L_D2", "D2_CrI_low", "D2_CrI_high", "D2_CrI_width", "D2_largest")
  stopifnot(all(need %in% names(fr)))
  stopifnot(fr$s_hat_AD >= 0, fr$s_hat_AD <= 1)
  stopifnot(fr$D2_CrI_low <= fr$D2_CrI_high)
  stopifnot(abs(fr$D2_CrI_width - (fr$D2_CrI_high - fr$D2_CrI_low)) < 1e-12)
  stopifnot(fr$s_hat_AD_lo <= fr$s_hat_AD_hi)
  stopifnot(is.logical(fr$D2_largest))
  ## P_D2 equals s_hat_D2_total (the degree fractions sum to one)
  stopifnot(abs(fr$P_D2 - fr$s_hat_D2_total) < 1e-8)

  ## --- share_pt -------------------------------------------------------------
  s_raw <- share_pt(y, KP, what = "raw", nIter = 400, burnIn = 100)
  stopifnot(length(s_raw) == 1L, is.finite(s_raw), s_raw >= 0, s_raw <= 1)
  s_pc <- share_pt(y, KP, PC = PC, what = "PC10", nIter = 400, burnIn = 100)
  stopifnot(is.finite(s_pc), s_pc >= 0, s_pc <= 1)
  s_fam <- share_pt(y, KP, Fk = Fk, what = "family", nIter = 400, burnIn = 100)
  stopifnot(is.finite(s_fam), s_fam >= 0, s_fam <= 1)

  ## --- partition_variance ---------------------------------------------------
  pv <- partition_variance(y, list(A = KP$A, D = KP$D), nIter = 400, burnIn = 100)
  stopifnot(identical(names(pv), c("A", "D")), abs(sum(pv) - 1) < 1e-8)

  ## --- run_rkhs -------------------------------------------------------------
  rr <- run_rkhs(y, list(A = list(K = KP$A, model = "RKHS"),
                         D = list(K = KP$D, model = "RKHS")),
                 c("A", "D"), nIter = 400, burnIn = 100, thin = 5)
  stopifnot(identical(colnames(rr$ci), c("A", "D")))
  stopifnot(identical(rownames(rr$ci), c("mean", "lo", "hi")))
  stopifnot(abs(sum(rr$ci["mean", ]) - 1) < 1e-8)
  stopifnot(all(rr$ci["lo", ] <= rr$ci["mean", ]), all(rr$ci["mean", ] <= rr$ci["hi", ]))
  cat("test-fit.R OK\n")
}

## --- regression: family_kernel with NA labels -----------------------------
## model.matrix() drops NA rows, which used to make the kernel smaller than
## length(ids) and fail with "length of 'dimnames' [1] not equal to array extent".
gna <- c("f1", "f2", NA, "f1", NA)
ids5 <- paste0("ind", seq_along(gna))
Fkna <- family_kernel(gna, ids5)
stopifnot(identical(dim(Fkna), c(5L, 5L)), identical(rownames(Fkna), ids5))
stopifnot(is.null(family_kernel(NULL, ids5)))

if (requireNamespace("BGLR", quietly = TRUE)) {
  set.seed(123)
  nf <- 30L; mf <- 60L
  Mf <- matrix(rbinom(nf * mf, 4, 0.4), nf, mf)
  rownames(Mf) <- paste0("ind", seq_len(nf))
  KPf <- build_panel(Mf, parse_ploidy("4", mf), need_X = TRUE)
  yf <- as.numeric(scale(KPf$X[[1]] %*% rnorm(mf))) + rnorm(nf, 0, 0.8)

  ## regression: burnIn >= nIter made (drop + 1):L count downwards, which
  ## produced out-of-range indices and an all-NaN CI instead of an error
  stopifnot(inherits(try(run_rkhs(yf, list(A = list(K = KPf$A, model = "RKHS"),
                                           D = list(K = KPf$D, model = "RKHS")),
                                  c("A", "D"), nIter = 200, burnIn = 400, thin = 5),
                         silent = TRUE), "try-error"))

  ## regression: a requested structure control with no input used to fall back
  ## to the raw share, reporting S = raw - raw = 0 for a control never fitted
  for (bad in list(list(what = "PC10", PC = NULL, Fk = NULL),
                   list(what = "family", PC = NULL, Fk = NULL),
                   list(what = "nope", PC = NULL, Fk = NULL))) {
    stopifnot(inherits(try(share_pt(yf, KPf, PC = bad$PC, Fk = bad$Fk, what = bad$what,
                                    nIter = 200, burnIn = 50), silent = TRUE), "try-error"))
  }
}

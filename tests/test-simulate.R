library(polyDx)
set.seed(21)
n <- 120L; m <- 400L
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))
KP <- build_panel(M, parse_ploidy("4", m), need_X = TRUE)

## --- sim_comp: standardized, lives in the contrast column space -----------
set.seed(1); g <- sim_comp(KP$X[[1]], KP$valid)
stopifnot(length(g) == n, abs(mean(g)) < 1e-8, abs(sd(g) - 1) < 1e-8)
set.seed(1); stopifnot(identical(g, sim_comp(KP$X[[1]], KP$valid)))   # seed-reproducible

## a zero contrast matrix degrades to zeros rather than NaN
stopifnot(all(sim_comp(matrix(0, n, 5), 5L) == 0))

## --- gen_additive_only ----------------------------------------------------
set.seed(2); y0 <- gen_additive_only(KP, H2 = 0.6)
stopifnot(length(y0) == n, all(is.finite(y0)))
set.seed(2); stopifnot(identical(y0, gen_additive_only(KP, H2 = 0.6)))

## H2 = 1 makes the phenotype the pure additive genetic value
set.seed(3); y1 <- gen_additive_only(KP, H2 = 1)
set.seed(3); gA <- sim_comp(KP$X[[1]], KP$valid)
stopifnot(max(abs(y1 - as.numeric(scale(gA)))) < 1e-8)

## --- gen_d2_spikein: the planted digenic part is orthogonal to additive ----
set.seed(4)
yd <- gen_d2_spikein(KP, s = 0.5, H2 = 1)
stopifnot(length(yd) == n, all(is.finite(yd)))
## at H2 = 1 and s = 1 the value is the additive-residualised digenic component
set.seed(5); ys1 <- gen_d2_spikein(KP, s = 1, H2 = 1)
set.seed(5)
gA <- sim_comp(KP$X[[1]], KP$valid); gD2 <- sim_comp(KP$X[[2]], KP$valid)
rD2 <- as.numeric(scale(residuals(lm(gD2 ~ gA))))
stopifnot(abs(cor(rD2, gA)) < 1e-8)                       # residualised => orthogonal
stopifnot(max(abs(ys1 - as.numeric(scale(rD2)))) < 1e-8)
## at s = 0 it collapses to the additive-only value
set.seed(6); ys0 <- gen_d2_spikein(KP, s = 0, H2 = 1)
set.seed(6); gA0 <- sim_comp(KP$X[[1]], KP$valid); invisible(sim_comp(KP$X[[2]], KP$valid))
stopifnot(max(abs(ys0 - as.numeric(scale(gA0)))) < 1e-8)

## --- gen_structure --------------------------------------------------------
PC <- panel_pc(KP, 10L)
grp <- rep(c("a", "b", "c", "d"), length.out = n)
set.seed(7); yp <- gen_structure(KP, PC, grp, lambda = 0.5, structure_type = "PC", H2 = 0.6)
stopifnot(length(yp) == n, all(is.finite(yp)))
set.seed(8); yf <- gen_structure(KP, PC, grp, lambda = 0.5, structure_type = "family", H2 = 0.6)
stopifnot(length(yf) == n, all(is.finite(yf)))
## lambda = 1, H2 = 1 with type "PC" reproduces standardized PC1
set.seed(9); ypc <- gen_structure(KP, PC, grp, lambda = 1, structure_type = "PC", H2 = 1)
stopifnot(max(abs(ypc - as.numeric(scale(as.numeric(scale(PC[, 1])))))) < 1e-8)
## lambda = 1, family: individuals in the same family share the value exactly
set.seed(10); yfam <- gen_structure(KP, PC, grp, lambda = 1, structure_type = "family", H2 = 1)
stopifnot(length(unique(round(yfam, 8))) == length(unique(grp)))

## --- simulate_pure_degree -------------------------------------------------
set.seed(11)
gv <- simulate_pure_degree(M, causal = sample(m, 100), degree = 3, ploidy = 4)
stopifnot(length(gv) == n, all(is.finite(gv)), sd(gv) > 0)
stopifnot(inherits(try(simulate_pure_degree(M, 1:10, degree = 9, ploidy = 4), silent = TRUE), "try-error"))
cat("test-simulate.R OK\n")

## --- regression: family structure with NA / factor group labels -----------
## eff[group] indexed a factor by level code rather than label, and NA labels
## propagated an NA phenotype; both are now normalised to character first.
KPg <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
PCg <- panel_pc(KPg, 5L)
gna <- rep(c("a", NA, "b"), length.out = nrow(M))
stopifnot(all(is.finite(gen_structure(KPg, PCg, gna, 0.5, "family"))))
gchr <- rep(c("a", "b"), length.out = nrow(M))
set.seed(4); ychr <- gen_structure(KPg, PCg, gchr, 0.5, "family")
set.seed(4); yfac <- gen_structure(KPg, PCg, factor(gchr), 0.5, "family")
stopifnot(isTRUE(all.equal(ychr, yfac)))

## --- regression: sim_comp's zero-variance guard must survive composition ---
## sim_comp returns rep(0, n) for a degenerate draw; the outer scale() then
## divided by a zero sd and handed BGLR an all-NaN phenotype.
KPz <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
KPz$X[[1]][] <- 0; KPz$X[[2]][] <- 0
stopifnot(all(sim_comp(KPz$X[[1]], KPz$valid) == 0))
stopifnot(!anyNA(gen_additive_only(KPz, H2 = 0.6)))
stopifnot(!anyNA(gen_d2_spikein(KPz, s = 0.5, H2 = 0.6)))

## structure_type is validated: a typo used to fall through to the family branch
KPs <- build_panel(M, parse_ploidy("4", ncol(M)), need_X = TRUE)
PCs <- panel_pc(KPs, 5L); grp <- rep(c("a", "b"), length.out = nrow(M))
stopifnot(inherits(try(gen_structure(KPs, PCs, grp, 0.5, "pc"), silent = TRUE), "try-error"))
stopifnot(inherits(try(gen_structure(KPs, PCs, NULL, 0.5, "family"), silent = TRUE), "try-error"))

## a single group level gives zeros, not NaN
one <- rep("only", nrow(M))
stopifnot(!anyNA(gen_structure(KPs, PCs, one, 0.5, "family")))

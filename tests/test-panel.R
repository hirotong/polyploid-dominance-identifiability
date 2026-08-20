library(polyDx)

## --- parse_ploidy: fixed --------------------------------------------------
p4 <- parse_ploidy("4", 100L)
stopifnot(identical(p4$mode, "fixed"), p4$ploidy == 4L, p4$ndeg == 4L, identical(p4$label, "4"))
p6 <- parse_ploidy("6", 100L)
stopifnot(p6$ndeg == 6L)
stopifnot(inherits(try(parse_ploidy("1", 10L), silent = TRUE), "try-error"))

## --- parse_ploidy: variable ----------------------------------------------
tmp <- tempfile(fileext = ".rds")
saveRDS(c(4L, 6L, 4L, 8L), tmp)
pv <- parse_ploidy(paste0("var:", tmp, ":6"), 4L)
stopifnot(identical(pv$mode, "variable"), identical(pv$mploidy, c(4L, 6L, 4L, 8L)))
stopifnot(pv$maxdeg == 4L)                # capped at min(mploidy)
stopifnot(identical(pv$label, "var<=4"))
stopifnot(inherits(try(parse_ploidy(paste0("var:", tmp, ":6"), 3L), silent = TRUE), "try-error"))
stopifnot(inherits(try(parse_ploidy("var:", 4L), silent = TRUE), "try-error"))

## --- subset_ploidy --------------------------------------------------------
pv2 <- subset_ploidy(pv, c(TRUE, FALSE, TRUE, TRUE))
stopifnot(identical(pv2$mploidy, c(4L, 4L, 8L)))
stopifnot(identical(subset_ploidy(p4, c(TRUE, FALSE))$ploidy, 4L))

## --- filter_markers -------------------------------------------------------
set.seed(3)
n <- 60; m <- 40
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
M[, 1] <- 0                                  # monomorphic -> dropped
M[, 2] <- rbinom(n, 4, 0.005)                # very low MAF -> dropped at 0.05
M[1:40, 3] <- NA                             # call rate 1/3 -> dropped at 0.8
fq <- filter_markers(M, parse_ploidy("4", m), maf_min = 0.05, call_min = 0.8, enabled = TRUE)
stopifnot(ncol(fq$M) < m, !anyNA(fq$M))
stopifnot(fq$qc$markers_pre_qc == m, fq$qc$markers_post_qc == ncol(fq$M))
stopifnot(fq$qc$markers_dropped_qc == m - ncol(fq$M))
stopifnot(fq$qc$marker_qc_enabled)
stopifnot(all(fq$M == round(fq$M)))          # rounded back to dosage classes
stopifnot(all(fq$M >= 0 & fq$M <= 4))

## enabled = FALSE keeps every marker but still imputes
fq0 <- filter_markers(M, parse_ploidy("4", m), enabled = FALSE)
stopifnot(ncol(fq0$M) == m, !anyNA(fq0$M))

## maf_min/call_min are validated (I3): out-of-range or malformed values error
## immediately instead of silently disabling the corresponding filter
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), maf_min = -1), silent = TRUE), "try-error"))
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), maf_min = 0.6), silent = TRUE), "try-error"))
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), maf_min = c(0.1, 0.2)), silent = TRUE), "try-error"))
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), maf_min = NA_real_), silent = TRUE), "try-error"))
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), call_min = -0.1), silent = TRUE), "try-error"))
stopifnot(inherits(try(filter_markers(M, parse_ploidy("4", m), call_min = 1.5), silent = TRUE), "try-error"))

## --- panel_from_matrix ----------------------------------------------------
rownames(M) <- paste0("id", seq_len(n))
phe <- data.frame(id = rownames(M), yield = rnorm(n), height = rnorm(n),
                  subpop = sample(letters[1:3], n, TRUE), stringsAsFactors = FALSE)
P <- panel_from_matrix(M, phe, "4", min_n = 40L)
stopifnot(P$n == n, P$markers == ncol(P$M))
stopifnot(setequal(P$traits, c("yield", "height")))   # 'subpop' dropped by the group regex
stopifnot(identical(rownames(P$M), rownames(P$phe)))
stopifnot(P$qc$raw_markers == m, P$qc$markers == ncol(P$M))

## a constant trait is dropped
phe2 <- phe; phe2$flat <- 1
P2 <- panel_from_matrix(M, phe2, "4", min_n = 40L)
stopifnot(!("flat" %in% P2$traits))

## max_traits caps the list
P3 <- panel_from_matrix(M, phe, "4", max_traits = 1L, min_n = 40L)
stopifnot(length(P3$traits) == 1L)

## group vector is carried through when it has >= 2 levels
grp <- sample(c("A", "B"), n, TRUE)
P4 <- panel_from_matrix(M, phe, "4", group = grp, min_n = 40L)
stopifnot(length(P4$group) == n, length(unique(P4$group)) == 2L)

## --- load_panel round-trips through RDS ----------------------------------
d <- tempfile(); dir.create(d)
saveRDS(M, file.path(d, "x_dosage.rds"))
saveRDS(phe, file.path(d, "x_pheno.rds"))
saveRDS(data.frame(id = rownames(M), grp = grp, stringsAsFactors = FALSE),
        file.path(d, "x_group.rds"))
PL <- load_panel(file.path(d, "x_dosage.rds"), file.path(d, "x_pheno.rds"), "4",
                 file.path(d, "x_group.rds"), min_n = 40L)
stopifnot(PL$n == n, setequal(PL$traits, c("yield", "height")), length(PL$group) == n)
stopifnot(identical(dim(PL$M), dim(P4$M)))

## dosage outside 0..ploidy is rejected
Mbad <- M; Mbad[1, 1] <- 9
saveRDS(Mbad, file.path(d, "bad_dosage.rds"))
stopifnot(inherits(try(load_panel(file.path(d, "bad_dosage.rds"),
                                  file.path(d, "x_pheno.rds"), "4", min_n = 40L),
                       silent = TRUE), "try-error"))
unlink(d, recursive = TRUE)
cat("test-panel.R OK\n")

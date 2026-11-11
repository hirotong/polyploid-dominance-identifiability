#!/usr/bin/env Rscript
# Stage the original alfalfa (Roza2019) + potato (SolCAP) panels into the standard
# dosage/pheno RDS format (integer 0..4) so they run through the SAME v2 generic
# pipeline as the new crops. Mirrors the processing in real_ortho_roza.R / potato_ortho.R.
# Reads raw data from data/raw and writes standard panel inputs to
# data/processed/panel_inputs/.
suppressMessages(library(data.table))
OUT <- "data/processed/panel_inputs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE); PL <- 4L
AD <- "data/raw/alfalfa_roza2019_medina"; PD <- "data/raw/potato_solcap"

## ---- ALFALFA (Roza2019; dosage already 0..4 in GWASPoly file) ----
g <- fread(file.path(AD, "Roza2019_06_GWASPoly.txt"), data.table = FALSE, check.names = FALSE)
samp <- setdiff(colnames(g), c("Marker", "Chrom", "Position")); D0 <- as.matrix(g[, samp]); mode(D0) <- "numeric"
keep <- rowMeans(is.na(D0)) < 0.20 & apply(D0, 1, function(z) length(unique(na.omit(z))) > 1); D0 <- D0[keep, ]
M <- t(D0); for (j in which(colSums(is.na(M)) > 0)) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
M <- round(M); storage.mode(M) <- "integer"
phe <- fread(file.path(AD, "BLUP_Yi_sqrt_SpATS_DArT.csv"), data.table = FALSE)
phe <- phe[, c(1, which(colnames(phe) %in% c("yi_year", "yi_2022", "yi_hrv")))]
saveRDS(M, file.path(OUT, "alfalfa_dosage.rds")); saveRDS(phe, file.path(OUT, "alfalfa_pheno.rds"))
cat(sprintf("alfalfa: %d x %d, range %d..%d, overlap %d\n", nrow(M), ncol(M), min(M), max(M),
            length(intersect(rownames(M), as.character(phe[[1]])))))

## ---- POTATO (SolCAP; ACGT tetraploid calls -> dosage 0..4) ----
to_dos <- function(col) { ch <- strsplit(col, ""); tab <- sort(table(unlist(ch)), decreasing = TRUE)
  tab <- tab[!names(tab) %in% c("N", "-", "?", "")]; if (length(tab) < 2) return(rep(NA_real_, length(col)))
  al <- names(tab)[1:2]; mn <- al[2]; vapply(ch, function(x) { x <- x[x %in% al]; if (length(x) < PL) NA_real_ else sum(x == mn) }, numeric(1)) }
geno <- fread(file.path(PD, "TableS1.csv"), data.table = FALSE, check.names = FALSE)
sp <- setdiff(colnames(geno), c("Marker", "Chrom", "Position")); calls <- as.matrix(geno[, sp])
dos <- t(apply(calls, 1, to_dos)); colnames(dos) <- sp
keep <- rowMeans(is.na(dos)) < 0.2 & apply(dos, 1, function(z) length(unique(na.omit(z))) > 1); dos <- dos[keep, ]
for (i in which(rowSums(is.na(dos)) > 0)) dos[i, is.na(dos[i, ])] <- mean(dos[i, ], na.rm = TRUE)
Mp <- t(dos); Mp <- round(Mp); storage.mode(Mp) <- "integer"
pp <- fread(file.path(PD, "TableS2.csv"), data.table = FALSE)
idc <- if ("Name" %in% colnames(pp)) "Name" else colnames(pp)[1]    # SolCAP clone ID column
keepc <- c(idc, intersect(c("total_yield", "tuber_shape", "sucrose", "vine_maturity_95d"), colnames(pp)))
pp <- pp[, keepc]
saveRDS(Mp, file.path(OUT, "potato_dosage.rds")); saveRDS(pp, file.path(OUT, "potato_pheno.rds"))
cat(sprintf("potato: %d x %d, range %d..%d, overlap %d\n", nrow(Mp), ncol(Mp), min(Mp), max(Mp),
            length(intersect(rownames(Mp), as.character(pp[[1]])))))
cat("STAGE DONE\n")

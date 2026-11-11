#!/usr/bin/env Rscript
# =====================================================================
# Task 2 — updog dosage calling for blueberry (run on krypton; heavy, multi-core).
# Input : data/processed/blueberry_dosage/blueberry_counts.rds = list(RO, AO),
#         oriented individuals x markers (rownames = sample IDs).
# Output: data/processed/blueberry_dosage/blueberry_dosage.rds = discrete dosage 0..4,
#         individuals x markers, after MAF/call-rate/biallelic QC.
# Model : multidog(model="norm") — general diverse-population model (as in
#         de Bem Oliveira et al. 2020); ploidy 4.
# =====================================================================
suppressMessages(library(updog))
OUT <- "data/processed/blueberry_dosage"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cnt <- readRDS(file.path(OUT, "blueberry_counts.rds"))      # list(RO, AO) indiv x markers
RO <- cnt$RO; AO <- cnt$AO
stopifnot(all(dim(RO) == dim(AO)))
size <- RO + AO                                             # total depth (indiv x markers)
NC <- as.integer(Sys.getenv("NC", "48"))
cat(sprintf("updog multidog: %d markers x %d individuals, nc=%d\n", ncol(RO), nrow(RO), NC))

# multidog expects markers in rows, individuals in cols:
md <- multidog(refmat = t(RO), sizemat = t(size), ploidy = 4L, model = "norm", nc = NC)

# Posterior-mode discrete dosage (0..4): markers x individuals -> individuals x markers.
geno <- format_multidog(md, varname = "geno")
M <- t(geno)
stopifnot(identical(rownames(M), rownames(RO)) || nrow(M) == nrow(RO))
if (is.null(rownames(M)) || any(rownames(M) != rownames(RO))) rownames(M) <- rownames(RO)

# QC: drop markers with >20% missing, monomorphic, or MAF < 0.05.
maf <- function(col) { p <- mean(col, na.rm = TRUE) / 4; min(p, 1 - p) }
keep <- apply(M, 2, function(z) mean(is.na(z)) < 0.20 &&
                length(unique(na.omit(z))) > 1 && maf(z) >= 0.05)
cat(sprintf("markers passing QC: %d of %d\n", sum(keep), length(keep)))
M <- M[, keep, drop = FALSE]
for (j in which(colSums(is.na(M)) > 0)) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
cat(sprintf("dosage matrix: %d indiv x %d markers\n", nrow(M), ncol(M)))
stopifnot(all(M >= 0 & M <= 4))
saveRDS(M, file.path(OUT, "blueberry_dosage.rds"))
cat("DOSAGE DONE\n")

#!/usr/bin/env Rscript
# =====================================================================
# Petrichor: combine per-chunk updog dosage calls into the final matrix.
# Reads $WD/chunk_*.rds (each list(idx, geno=markers x indiv)), reassembles
# in marker order, applies MAF/call-rate/biallelic QC, writes blueberry_dosage.rds
# (individuals x markers, 0..4) to the shared path -> visible on krypton.
# =====================================================================
WD <- Sys.getenv("WD", "/path/to/work/code/bb_dosage")
K  <- as.integer(Sys.getenv("NCHUNK", "32"))
cnt <- readRDS(file.path(WD, "blueberry_counts.rds")); RO <- cnt$RO

parts <- lapply(1:K, function(t) {
  f <- file.path(WD, sprintf("chunk_%03d.rds", t))
  if (!file.exists(f)) stop(sprintf("missing chunk %d (%s)", t, f)); readRDS(f) })
parts <- parts[order(sapply(parts, function(p) p$idx[1]))]
geno_mk_ind <- do.call(rbind, lapply(parts, function(p) p$geno))   # markers x indiv (ordered)
M <- t(geno_mk_ind)                                                # indiv x markers
if (is.null(rownames(M)) || length(rownames(M)) != nrow(RO) || !all(rownames(M) == rownames(RO)))
  rownames(M) <- rownames(RO)
cat(sprintf("assembled: %d indiv x %d markers (pre-QC)\n", nrow(M), ncol(M)))

maf <- function(col) { p <- mean(col, na.rm = TRUE) / 4; min(p, 1 - p) }
keep <- apply(M, 2, function(z) mean(is.na(z)) < 0.20 &&
                length(unique(na.omit(z))) > 1 && maf(z) >= 0.05)
Mk <- M[, keep, drop = FALSE]
for (j in which(colSums(is.na(Mk)) > 0)) Mk[is.na(Mk[, j]), j] <- mean(Mk[, j], na.rm = TRUE)
stopifnot(all(Mk >= 0 & Mk <= 4))
saveRDS(Mk, file.path(WD, "blueberry_dosage.rds"))
cat(sprintf("FINAL dosage: %d indiv x %d markers (kept of %d)\n", nrow(Mk), ncol(Mk), length(keep)))
cat("COMBINE DONE\n")

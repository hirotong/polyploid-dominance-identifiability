#!/usr/bin/env Rscript
# =====================================================================
# Build standard-format dosage + phenotype RDS for blueberry panel #2
# (Ferrao et al., Dryad doi:10.5061/dryad.kd4jq6h, autotetraploid Vaccinium).
#
# The tetraploid SNP file is ALREADY discrete dosage 0..4 (pre-filtered),
# NOT ref/alt read counts. Per README:
#   Tetraploid_SNPdata.csv : 80591 markers x 1559 indiv; cols 1-3 =
#       snp.name, chr.name, pos.name; remaining cols = genotype dosages.
#   Tetraploid_phenotypic_data.csv : Sample + 8 traits.
# So we round + orient (indiv x markers) + apply QC; no updog needed.
#
# QC: biallelic (range 0..4), call-rate > 0.8, MAF >= 0.05; mean-impute
# any residual missing (none expected — data are complete).
#
# Output (standard format):
#   blueberry2_dosage.rds : matrix individuals x markers, integer 0..4,
#                           rownames = sample IDs.
#   blueberry2_pheno.rds  : data.frame, FIRST col = Sample ID, rest numeric traits.
# =====================================================================
suppressMessages(library(data.table))

DIR <- "data/raw/blueberry_ferrao_kd4jq6h/extracted/data_set"
OUT <- "data/processed/panel_inputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- genotypes ----
cat("reading SNP matrix (large CSV)...\n")
snp <- fread(file.path(DIR, "Tetraploid_SNPdata.csv"), header = TRUE,
             showProgress = FALSE)
cat(sprintf("raw SNP table: %d rows (markers) x %d cols\n", nrow(snp), ncol(snp)))
meta_cols <- c("snp.name", "chr.name", "pos.name")
stopifnot(all(meta_cols %in% colnames(snp)))
marker_ids <- snp[["snp.name"]]
geno_cols  <- setdiff(colnames(snp), meta_cols)            # sample IDs
G <- as.matrix(snp[, ..geno_cols])                         # markers x individuals
storage.mode(G) <- "double"
rownames(G) <- marker_ids
rm(snp); gc()

# Orient individuals x markers (markers >> individuals -> transpose).
M <- t(G)                                                  # individuals x markers
rm(G); gc()
cat(sprintf("oriented dosage: %d individuals x %d markers\n", nrow(M), ncol(M)))

# Round to integer dosage and clamp to 0..4 (values already discrete).
M <- round(M)
M[M < 0] <- NA; M[M > 4] <- NA
cat(sprintf("value range pre-QC: [%g, %g]; missing cells: %d\n",
            min(M, na.rm = TRUE), max(M, na.rm = TRUE), sum(is.na(M))))

## ---- QC: call-rate > 0.8, polymorphic, MAF >= 0.05 ----
maf <- function(col) { p <- mean(col, na.rm = TRUE) / 4; min(p, 1 - p) }
callrate <- 1 - colMeans(is.na(M))
mafs     <- apply(M, 2, maf)
polymorphic <- apply(M, 2, function(z) length(unique(na.omit(z))) > 1)
keep <- callrate > 0.8 & polymorphic & mafs >= 0.05 & is.finite(mafs)
cat(sprintf("QC: %d of %d markers pass (callrate>0.8 & MAF>=0.05 & polymorphic)\n",
            sum(keep), length(keep)))
M <- M[, keep, drop = FALSE]

# Mean-impute any residual missing (round to nearest integer dosage).
nmiss <- sum(is.na(M))
if (nmiss > 0) {
  cat(sprintf("mean-imputing %d residual missing cells\n", nmiss))
  for (j in which(colSums(is.na(M)) > 0)) {
    mu <- mean(M[, j], na.rm = TRUE); M[is.na(M[, j]), j] <- round(mu)
  }
}
storage.mode(M) <- "integer"
stopifnot(all(M >= 0 & M <= 4), !any(is.na(M)))
cat(sprintf("final dosage matrix: %d indiv x %d markers; range [%d, %d]\n",
            nrow(M), ncol(M), min(M), max(M)))

## ---- phenotypes ----
phe <- fread(file.path(DIR, "Tetraploid_phenotypic_data.csv"), header = TRUE,
             showProgress = FALSE)
phe <- as.data.frame(phe)
colnames(phe)[1] <- "Sample"
phe$Sample <- as.character(phe$Sample)
# Coerce trait columns to numeric.
for (cc in setdiff(colnames(phe), "Sample")) phe[[cc]] <- as.numeric(phe[[cc]])
cat(sprintf("phenotype: %d rows; traits: %s\n", nrow(phe),
            paste(setdiff(colnames(phe), "Sample"), collapse = ", ")))

## ---- ID overlap check ----
ov <- intersect(rownames(M), phe$Sample)
cat(sprintf("ID overlap (dosage rownames vs pheno Sample): %d\n", length(ov)))
stopifnot(length(ov) >= 50)

## ---- write standard-format RDS ----
saveRDS(M, file.path(OUT, "blueberry2_dosage.rds"))
saveRDS(phe, file.path(OUT, "blueberry2_pheno.rds"))
cat("first 3 sample IDs:", paste(head(rownames(M), 3), collapse = ", "), "\n")
cat("first 3 marker IDs:", paste(head(colnames(M), 3), collapse = ", "), "\n")
cat("BUILD DONE\n")

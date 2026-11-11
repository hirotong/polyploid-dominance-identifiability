# experiments/data_preparation/blueberry_dosage/01_parse_blueberry.R
# Parse de Bem Oliveira et al. 2020 autotetraploid blueberry data (Dryad
# doi:10.5061/dryad.8pk0p2nk9). Produces RO/AO read-count matrices oriented as
# individuals x markers, plus a phenotype data.frame (Sample + 3 eBLUE traits).
#
# Object structure found in the .RData files (Step 2 inspection):
#   GenotypeInfo_RO_AO.RData:
#     ROfiltered : matrix 87628 markers (rows) x 1804 individuals (cols)
#     AOfiltered : matrix 87628 markers (rows) x 1804 individuals (cols)
#       rownames = marker IDs (e.g. VaccDscaff1_65041)
#       colnames = sample IDs  (e.g. UFF_192502_P01_WA01)
#   LSMeans_weight_firmness_yield.RData (three 1804x2 data.frames):
#     eBLUEs_firm   : Sample, eblues_firmness
#     eBLUEs_weight : Sample, eblues_weight
#     eBLUEs_yield  : Sample, eblues_yield

DIR <- "data/raw/blueberry_debemoliveira2020"; OUT <- "data/processed/blueberry_dosage"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
eg <- new.env(); load(file.path(DIR, "GenotypeInfo_RO_AO.RData"), envir = eg)
ep <- new.env(); load(file.path(DIR, "LSMeans_weight_firmness_yield.RData"), envir = ep)

## >>> EDIT to match Step-2 names <<<
RO <- get("ROfiltered", eg)   # ref read counts  (markers x individuals)
AO <- get("AOfiltered", eg)   # alt read counts  (markers x individuals)

# Phenotype is split across three single-trait data.frames; merge on Sample.
phe_firm   <- get("eBLUEs_firm",   ep)   # Sample, eblues_firmness
phe_weight <- get("eBLUEs_weight", ep)   # Sample, eblues_weight
phe_yield  <- get("eBLUEs_yield",  ep)   # Sample, eblues_yield
phe <- Reduce(function(a, b) merge(a, b, by = "Sample", all = TRUE),
              list(phe_weight, phe_firm, phe_yield))
# Rename to clean trait names; keep Sample as the first (ID) column.
colnames(phe) <- c("Sample", "weight", "firmness", "yield")
## <<< END EDIT >>>

stopifnot(all(dim(RO) == dim(AO)))
stopifnot(identical(dimnames(RO), dimnames(AO)))

# Orient as individuals x markers. Here markers (87628) >> individuals (1804),
# so the matrices have more rows than cols and must be transposed.
if (ncol(RO) < nrow(RO)) { RO <- t(RO); AO <- t(AO) }   # orient: individuals x markers
cat(sprintf("genotypes: %d individuals x %d markers\n", nrow(RO), ncol(RO)))

# After transpose, rownames(RO) are the sample IDs and colnames(RO) are markers.
stopifnot(!is.null(rownames(RO)), !is.null(colnames(RO)))
cat("first sample IDs:", paste(head(rownames(RO), 3), collapse = ", "), "\n")
cat("first marker IDs:", paste(head(colnames(RO), 3), collapse = ", "), "\n")

# Sanity-check ID overlap with phenotypes.
ov <- length(intersect(rownames(RO), as.character(phe$Sample)))
cat(sprintf("ID overlap (geno rownames vs pheno Sample): %d\n", ov))
if (ov == 0) stop("No ID overlap between genotype rownames and phenotype Sample column")

saveRDS(list(RO = RO, AO = AO), file.path(OUT, "blueberry_counts.rds"))
saveRDS(phe, file.path(OUT, "blueberry_pheno.rds"))
cat("traits:", paste(setdiff(colnames(phe), "Sample"), collapse = ", "), "\n")
cat("PARSE DONE\n")

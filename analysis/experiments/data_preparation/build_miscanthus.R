#!/usr/bin/env Rscript
# =====================================================================
# build_miscanthus.R -- Miscanthus sacchariflorus autotetraploid panel.
#
# SOURCE (current)
#   Illinois Data Bank IDB-4517068 / DOI 10.13012/B2IDB-4829913_V2
#   Njuguna et al. 2023, The Plant Genome 16(4):e20401,
#   "Impact of genotype-calling methodologies on genome-wide association and
#    genomic prediction in polyploids".
#   Genotypes : updog_discrete_Miscanthus_tetraploid_genotype.RData
#               (object `updog_mode_4x_geno`; 226 accessions x 276,159 SNPs,
#                discrete updog calls made AT PLOIDY 4, support {0,1,2,3,4})
#   Phenotypes: Miscanthus_phenotypes.csv -- author-supplied LS-means for
#               culm length, biomass yield, heading date (already transformed).
#   Same dataset and the same `updog discrete` variant as the BB-UF-call*
#   blueberry panel (Updog_discrete_Blueberry_genotype.RData), so the two
#   genotype-calling panels stay comparable.
#
# WHY THIS SCRIPT EXISTS (the bug it replaces)
#   The panel previously came from IDB-6439486 (Njuguna et al., the Miscanthus
#   GWAS/GP dataset) as a finished .rds with no build script -- the only panel
#   in the study without one, which is why the following went unreviewed.
#   That source is MIXED PLOIDY: it calls each accession at that accession's
#   own ploidy and puts everything on a common 0-2 scale.  Doubling it to 0-4
#   therefore yields genuine {0,1,2,3,4} calls for the tetraploid accessions
#   but doubled-diploid {0,2,4} rows for the diploid majority.  Measured on the
#   shipped matrix: 348 of 601 accessions (58%) had stride-2 support, i.e. were
#   diploids analysed at ploidy 4.  Those rows carry no degree-3/4 information
#   -- on a 3-point support the single-locus function space is 3-dimensional,
#   so D3 and D4 are exact linear combinations of {intercept, A, D2}.
#   The panel's reported four-kernel effective rank (1.199, the lowest of the
#   11 panels) was measuring that ploidy mixture, not autotetraploid geometry.
#
#   Note the pooled support of the shipped matrix was a complete {0,1,2,3,4}:
#   a whole-matrix check passes.  Only a PER-INDIVIDUAL support check catches
#   it, which is what assert_polyploid_support() now does.
#
# Usage:  Rscript experiments/data_preparation/build_miscanthus.R
# Env:    MISCANTHUS_RAW  raw dir (default data/raw/miscanthus_gc)
#         PANEL_OUT       output dir (default data/processed/panel_inputs)
#         MISCANTHUS_META optional GWAS metadata .xlsx for the DAPC group file
#                         and an independent ploidy cross-check
# Writes: <PANEL_OUT>/miscanthus_dosage.rds
#         <PANEL_OUT>/miscanthus_pheno.rds
#         <PANEL_OUT>/miscanthus_group.rds   (if metadata available)
# =====================================================================
source("experiments/data_preparation/lib_dosage_checks.R")

PLOIDY <- 4L
RAW    <- Sys.getenv("MISCANTHUS_RAW", "data/raw/miscanthus_gc")
OUT    <- Sys.getenv("PANEL_OUT", "data/processed/panel_inputs")
META   <- Sys.getenv("MISCANTHUS_META", "data/raw/miscanthus/MSA_phenotype_Metadata_GWAS.xlsx")
GENO   <- file.path(RAW, "updog_discrete_Miscanthus_tetraploid_genotype.RData")
PHENO  <- file.path(RAW, "Miscanthus_phenotypes.csv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
for (f in c(GENO, PHENO)) if (!file.exists(f)) stop("missing raw file: ", f, call. = FALSE)

## ---- genotypes ---------------------------------------------------------
e <- new.env()
load(GENO, envir = e)
objs <- ls(e)
cat("objects in RData:", paste(objs, collapse = ", "), "\n")
pick <- NULL; best <- -1
for (o in objs) {
  v <- get(o, envir = e)
  if ((is.matrix(v) || is.data.frame(v)) && length(v) > best) { best <- length(v); pick <- o }
}
if (is.null(pick)) stop("no matrix-like object in ", GENO, call. = FALSE)
G <- get(pick, envir = e); rm(e)
if (is.data.frame(G)) G <- as.matrix(G)
storage.mode(G) <- "double"
cat(sprintf("using object: %s  dim: %d x %d\n", pick, nrow(G), ncol(G)))

## ---- phenotypes --------------------------------------------------------
ph <- read.csv(PHENO, check.names = FALSE, stringsAsFactors = FALSE)
colnames(ph)[1] <- "Sample"
ph$Sample <- as.character(ph$Sample)
cat(sprintf("phenotypes: %d accessions x %d traits (%s)\n",
            nrow(ph), ncol(ph) - 1L, paste(colnames(ph)[-1], collapse = ", ")))

## ---- orient to individuals x markers -----------------------------------
ov_row <- length(intersect(rownames(G), ph$Sample))
ov_col <- length(intersect(colnames(G), ph$Sample))
cat(sprintf("sample-ID overlap: rows=%d cols=%d\n", ov_row, ov_col))
if (ov_col > ov_row) { G <- t(G); cat("transposed -> individuals x markers\n") }
if (is.null(rownames(G))) stop("genotype matrix has no sample IDs", call. = FALSE)

## ---- value scale -------------------------------------------------------
# updog discrete calls are already integer dosages 0..ploidy.  Guard against a
# source that ships allele proportions (0..1) or a lower-ploidy scale (0..2) --
# the 0..2 case is exactly how the previous panel went wrong, so if the range
# implies rescaling we refuse rather than silently multiply.
rng <- range(G, na.rm = TRUE)
cat(sprintf("value range: %g to %g\n", rng[1], rng[2]))
if (rng[2] <= PLOIDY / 2 + 1e-9) {
  stop(sprintf(paste0("genotype values max out at %g, below ploidy %d. These are not ploidy-%d ",
                      "calls; rescaling them here is the exact mistake this rebuild corrects."),
               rng[2], PLOIDY, PLOIDY), call. = FALSE)
}
G <- round(G)
storage.mode(G) <- "integer"

## =====================================================================
## HARD ASSERTION -- the check that was missing.
## Whole-matrix support, plus per-individual support so a mixed-ploidy panel
## cannot pass by pooling real 4x rows with doubled-diploid rows.
## =====================================================================
chk <- assert_polyploid_support(G, PLOIDY, label = "miscanthus")
cat(sprintf("SUPPORT OK: {%s}, rank %d/%d; all %d individuals stride 1 (genuine ploidy %d)\n",
            paste(chk$support, collapse = ","), chk$support_rank, chk$ideal_rank,
            nrow(G), PLOIDY))
cat("dosage class frequencies:\n"); print(round(chk$class_freq, 5))
odd <- rowMeans(matrix(G %% 2L == 1L, nrow(G)), na.rm = TRUE)
cat(sprintf("simplex/triplex fraction per individual: min %.4f median %.4f max %.4f\n",
            min(odd), median(odd), max(odd)))

## ---- independent ploidy cross-check + DAPC groups ----------------------
# The companion GWAS dataset (IDB-6439486) carries a per-accession Ploidy call
# and a DAPC genetic-group assignment.  Neither is used to build the matrix;
# Ploidy is a fully independent confirmation that these accessions are 4x, and
# DAPC supplies the group file the panel workflow's family arm expects.
grp <- NULL
if (file.exists(META) && requireNamespace("readxl", quietly = TRUE)) {
  md <- suppressMessages(readxl::read_excel(META, sheet = "LS means of traits for GWAS", skip = 1))
  md <- as.data.frame(md)
  hit <- match(rownames(G), md$Accession_ID)
  pl  <- md$Ploidy[hit]
  cat(sprintf("independent ploidy cross-check (%s): %d/%d accessions found\n",
              basename(META), sum(!is.na(hit)), nrow(G)))
  print(table(ifelse(is.na(pl), "not-in-GWAS-metadata", pl)))
  if (any(!is.na(pl) & pl != "4x")) {
    stop(sprintf("%d accessions are labelled non-4x in the GWAS metadata: %s",
                 sum(!is.na(pl) & pl != "4x"),
                 paste(rownames(G)[!is.na(pl) & pl != "4x"], collapse = ", ")), call. = FALSE)
  }
  dapc <- md$DAPC[hit]
  if (sum(!is.na(dapc)) > 0 && length(unique(dapc[!is.na(dapc)])) >= 2) {
    grp <- data.frame(ID = rownames(G), group = as.character(dapc), stringsAsFactors = FALSE)
  }
} else {
  cat("NOTE: GWAS metadata unavailable -> no ploidy cross-check, no group file\n")
}

## ---- align genotypes and phenotypes ------------------------------------
# Marker MAF / call-rate QC is deliberately NOT done here: the panel workflow
# and the prediction CV both apply their own (lib_marker_qc.R), and doing it
# twice would make raw_markers vs markers_pre_qc unreadable in the audit.
common <- intersect(rownames(G), ph$Sample)
cat(sprintf("ID overlap geno(%d) n pheno(%d) = %d\n", nrow(G), nrow(ph), length(common)))
if (length(common) < 50L) stop("fewer than 50 accessions with both genotype and phenotype", call. = FALSE)
G  <- G[common, , drop = FALSE]
ph <- ph[match(common, ph$Sample), , drop = FALSE]
rownames(ph) <- NULL
if (!is.null(grp)) grp <- grp[match(common, grp$ID), , drop = FALSE]

## ---- re-assert after subsetting ---------------------------------------
chk2 <- assert_polyploid_support(G, PLOIDY, label = "miscanthus-final")
stopifnot(identical(chk2$support, 0:PLOIDY))

## ---- write -------------------------------------------------------------
saveRDS(G,  file.path(OUT, "miscanthus_dosage.rds"))
saveRDS(ph, file.path(OUT, "miscanthus_pheno.rds"))
if (!is.null(grp)) {
  rownames(grp) <- NULL
  saveRDS(grp, file.path(OUT, "miscanthus_group.rds"))
  cat(sprintf("WROTE miscanthus_group.rds (%d accessions, %d DAPC groups, %d unassigned)\n",
              nrow(grp), length(unique(na.omit(grp$group))), sum(is.na(grp$group))))
  print(table(grp$group, useNA = "ifany"))
}
cat(sprintf("WROTE miscanthus_dosage.rds (%d x %d, ploidy %d) and miscanthus_pheno.rds (%d x %d)\n",
            nrow(G), ncol(G), PLOIDY, nrow(ph), ncol(ph)))
cat("traits:", paste(colnames(ph)[-1], collapse = ", "), "\n")
print(summary(ph[, -1, drop = FALSE]))
cat("BUILD MISCANTHUS DONE\n")

#!/usr/bin/env Rscript
# Rose (tetraploid) preprocessing: SNP-chip dosage (0..4) + in-vitro rooting traits.
# Source: data.uni-hannover.de  Wamhoff et al. 2024
# Genotype CSV layout: markers in ROWS, samples in COLUMNS.
#   col1 = customer_id_PQ (marker id), col2 = direction/Chr, col3 = LG/position,
#   cols 4..end = 106 rose genotype sample codes (e.g. S_AK, AD, BB), values dosage 0..4 (+NA)
# Phenotype: aggregate replicate-level AR rooting + RP anatomy to genotype means.
suppressMessages({library(data.table)})

base   <- "data/raw/rose"
outdir <- "data/processed/panel_inputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- genotype ----
g <- fread(file.path(base, "geno_gwas_ms3_foruse.csv"), data.table = FALSE, na.strings = c("NA",""))
marker_id <- g[[1]]
samp_cols <- colnames(g)[4:ncol(g)]          # 106 sample codes
G <- as.matrix(g[, 4:ncol(g)])               # markers x samples
mode(G) <- "numeric"
rownames(G) <- marker_id
# transpose -> individuals x markers
G <- t(G)                                     # samples x markers
cat("raw dosage matrix:", nrow(G), "indiv x", ncol(G), "markers\n")

# QC: dosage range, call-rate>0.8, MAF>=0.05 ; mean-impute residual missing
stopifnot(all(G[!is.na(G)] %in% 0:4))
callrate <- colMeans(!is.na(G))
G <- G[, callrate > 0.8, drop = FALSE]
af   <- colMeans(G, na.rm = TRUE) / 4         # allele freq (dosage/ploidy)
maf  <- pmin(af, 1 - af)
G <- G[, maf >= 0.05, drop = FALSE]
# mean-impute residual missing, round to nearest dosage 0..4
for (j in which(colSums(is.na(G)) > 0)) {
  m <- mean(G[, j], na.rm = TRUE); G[is.na(G[, j]), j] <- m
}
G <- round(G); storage.mode(G) <- "integer"
G[G < 0] <- 0L; G[G > 4] <- 4L
cat("post-QC dosage matrix:", nrow(G), "indiv x", ncol(G), "markers\n")

## ---- phenotype ----
# XLSX pre-exported to CSV (ar_raw.csv, rp_raw.csv) since readxl is absent in r_env.
ar <- fread(file.path(base, "ar_raw.csv"), data.table = FALSE, na.strings = c("NA",""))
for (cc in c("7d_root","14d_root","21d_root","21d_rootNr"))
  ar[[cc]] <- suppressWarnings(as.numeric(ar[[cc]]))
ar_agg <- aggregate(cbind(`7d_root`, `14d_root`, `21d_root`, `21d_rootNr`) ~ ID,
                    data = ar, FUN = function(x) mean(x, na.rm = TRUE))
colnames(ar_agg) <- c("ID","root_rate_7d","root_rate_14d","root_rate_21d","root_number_21d")

rp <- fread(file.path(base, "rp_raw.csv"), data.table = FALSE, na.strings = c("NA",""))
for (cc in c("RP_total","CrossSection_Area"))
  rp[[cc]] <- suppressWarnings(as.numeric(rp[[cc]]))
rp_agg <- aggregate(cbind(RP_total, CrossSection_Area) ~ ID,
                    data = rp, FUN = function(x) mean(x, na.rm = TRUE))

ph <- merge(ar_agg, rp_agg, by = "ID", all = TRUE)
colnames(ph)[1] <- "Sample"

## ---- intersect ----
common <- intersect(rownames(G), ph$Sample)
cat("ID overlap geno∩pheno:", length(common), "\n")
G  <- G[common, , drop = FALSE]
ph <- ph[match(common, ph$Sample), , drop = FALSE]
rownames(ph) <- NULL

saveRDS(G,  file.path(outdir, "rose_dosage.rds"))
saveRDS(ph, file.path(outdir, "rose_pheno.rds"))
cat("WROTE rose_dosage.rds (", nrow(G), "x", ncol(G), ") and rose_pheno.rds (",
    nrow(ph), "x", ncol(ph), ")\n")
cat("traits:", paste(colnames(ph)[-1], collapse=", "), "\n")
print(summary(ph[,-1]))

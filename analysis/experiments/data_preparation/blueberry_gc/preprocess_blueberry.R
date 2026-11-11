#!/usr/bin/env Rscript
# Blueberry (tetraploid, Vaccinium corymbosum) preprocessing.
# Source: Illinois Data Bank IDB-4517068 (Njuguna et al., genotype-calling methods).
# Genotype: Updog_discrete_Blueberry_genotype.RData = discrete dosage calls 0..4 (updog).
# Phenotype: Blueberry_phenotypes.csv (LS-means): Fruit weight, Firmness, Scar diameter.
suppressMessages({library(data.table)})

base   <- "data/raw/blueberry_gc"
outdir <- "data/processed/panel_inputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- load genotype RData (object name unknown -> inspect) ----
e <- new.env()
load(file.path(base, "Updog_discrete_Blueberry_genotype.RData"), envir = e)
objs <- ls(e)
cat("objects in RData:", paste(objs, collapse=", "), "\n")
# pick the largest matrix/data.frame-like object
pick <- NULL; best <- -1
for (o in objs) {
  v <- get(o, envir = e)
  sz <- tryCatch(length(v), error = function(x) 0)
  if (is.matrix(v) || is.data.frame(v)) if (sz > best) { best <- sz; pick <- o }
}
if (is.null(pick)) pick <- objs[1]
G <- get(pick, envir = e)
cat("using object:", pick, " class:", paste(class(G), collapse="/"), " dim:", paste(dim(G), collapse="x"), "\n")
if (is.data.frame(G)) G <- as.matrix(G)
mode(G) <- "numeric"
cat("rownames head:", paste(head(rownames(G)), collapse=", "), "\n")
cat("colnames head:", paste(head(colnames(G)), collapse=", "), "\n")

## ---- orient to individuals x markers ----
# Sample IDs look like 'UFF_192502_P01_WA01'. Detect which dim carries them.
ph <- fread(file.path(base, "Blueberry_phenotypes.csv"), data.table = FALSE)
ph_ids <- ph[[1]]
ov_row <- length(intersect(rownames(G), ph_ids))
ov_col <- length(intersect(colnames(G), ph_ids))
cat("sample-ID overlap: rows=", ov_row, " cols=", ov_col, "\n")
if (ov_col > ov_row) { G <- t(G); cat("transposed -> individuals x markers\n") }

## ---- value scale: ensure dosage 0..4 ----
rng <- range(G, na.rm = TRUE)
cat("value range:", rng[1], "to", rng[2], "\n")
# if calls are 0..1 ref-allele proportion, scale; updog discrete should already be 0..4 integers
if (rng[2] <= 1.5) { G <- G * 4; cat("scaled proportion -> dosage\n") }
G <- round(G)

## ---- QC ----
callrate <- colMeans(!is.na(G))
G <- G[, callrate > 0.8, drop = FALSE]
af  <- colMeans(G, na.rm = TRUE) / 4
maf <- pmin(af, 1 - af)
G <- G[, maf >= 0.05, drop = FALSE]
for (j in which(colSums(is.na(G)) > 0)) {
  m <- mean(G[, j], na.rm = TRUE); G[is.na(G[, j]), j] <- m
}
G <- round(G); storage.mode(G) <- "integer"
G[G < 0] <- 0L; G[G > 4] <- 4L
cat("post-QC dosage matrix:", nrow(G), "indiv x", ncol(G), "markers\n")

## ---- intersect ----
colnames(ph)[1] <- "Sample"
common <- intersect(rownames(G), ph$Sample)
cat("ID overlap geno∩pheno:", length(common), "\n")
G  <- G[common, , drop = FALSE]
ph <- ph[match(common, ph$Sample), , drop = FALSE]
rownames(ph) <- NULL

saveRDS(G,  file.path(outdir, "blueberry_gc_dosage.rds"))
saveRDS(ph, file.path(outdir, "blueberry_gc_pheno.rds"))
cat("WROTE blueberry_gc_dosage.rds (", nrow(G), "x", ncol(G), ") and blueberry_gc_pheno.rds (",
    nrow(ph), "x", ncol(ph), ")\n")
cat("traits:", paste(colnames(ph)[-1], collapse=", "), "\n")
print(summary(ph[,-1]))

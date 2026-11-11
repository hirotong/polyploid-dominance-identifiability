#!/usr/bin/env Rscript
# =====================================================================
# Orthogonal A+D2+..+Dm decomposition under PHENOTYPE PERMUTATION (null).
# For one primary trait per panel: fit the orthogonal model on the REAL y
# (perm 0) and on NP permuted y; record the per-degree variance FRACTION each
# time. This shows what the orthogonal decomposition reports when there is NO
# genotype-phenotype signal, and whether the real split is distinguishable
# from that null.  (Complements the standard-model _null.csv which only tracks
# the A+D Vd-share under permutation.)
#
# Usage: run_ortho_perm.R <dosage_rds> <pheno_rds> <out_prefix> <trait|auto> <ploidy> [nperm=30]
#   <ploidy>: 4 | 6 | "auto" | "varploidy:<mploidy_rds>"
# Output: {prefix}_ortho_perm.csv  (cols: perm, trait, degree, frac)
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho_mp.R")   # -> build_kernels*, ortho_etanames, BGLR
qc_lib <- Sys.getenv("MARKER_QC_LIB", "experiments/empirical/multi/lib_marker_qc.R")
if (!file.exists(qc_lib) && file.exists("lib_marker_qc.R")) qc_lib <- "lib_marker_qc.R"
source(qc_lib)
a <- commandArgs(TRUE); stopifnot(length(a) >= 5)
dosage <- a[1]; pheno <- a[2]; prefix <- a[3]; trait_arg <- a[4]; ploidy_arg <- a[5]
NP <- if (length(a) >= 6) as.integer(a[6]) else 30L
OUTDIR <- dirname(prefix); if (nzchar(OUTDIR)) dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

M <- readRDS(dosage); phe <- readRDS(pheno); if (is.matrix(phe)) phe <- as.data.frame(phe)
idc <- colnames(phe)[1]; rownames(phe) <- as.character(phe[[idc]])
common <- intersect(rownames(M), rownames(phe)); stopifnot(length(common) >= 50)
M <- M[common, , drop = FALSE]; phe <- phe[common, ]
TR <- setdiff(colnames(phe), idc)
drop <- grepl("group|gen[ae]tic|subpop|^pop$|cluster|^id$|clone|geno", TR, ignore.case = TRUE)
if (any(drop)) cat("dropping non-phenotype columns:", paste(TR[drop], collapse = ", "), "\n")
TR <- TR[!drop]
TR <- TR[sapply(TR, function(t) suppressWarnings(sd(as.numeric(phe[[t]]), na.rm = TRUE)) > 0)]
trait <- if (trait_arg == "auto") TR[1] else trait_arg
stopifnot(trait %in% colnames(phe))

# ---- kernels (same construction as the attribution runners) ----
if (grepl("^varploidy:", ploidy_arg)) {
  mp <- readRDS(sub("^varploidy:", "", ploidy_arg)); stopifnot(length(mp) == ncol(M))
  fq <- apply_marker_qc(M, ploidy = NA_integer_, marker_ploidy = mp)
  M <- fq$M; mp <- fq$marker_ploidy
  maxdeg <- min(6L, min(mp)); K <- build_kernels_varm(M, mp, maxdeg); ndeg <- maxdeg
} else {
  PL <- if (ploidy_arg == "auto") max(round(M), na.rm = TRUE) else as.integer(ploidy_arg)
  fq <- apply_marker_qc(M, ploidy = PL)
  M <- fq$M
  stopifnot(all(M >= 0 & M <= PL, na.rm = TRUE))
  if (PL == 4L) { Kb <- build_kernels(M); K <- list(O = list(Kb$OA, Kb$OD2, Kb$OD3, Kb$OD4)); ndeg <- 4L }
  else { Kx <- build_kernels_m(M, PL); K <- list(O = Kx$O); ndeg <- PL }
}
en <- ortho_etanames(ndeg); names(K$O) <- en
write_marker_qc_audit(prefix, basename(prefix), fq$qc)
ETA <- function() { e <- list(); for (k in seq_len(ndeg)) e[[en[k]]] <- list(K = K$O[[k]], model = "RKHS"); e }
cat(sprintf("[%s] perm-null: %d indiv, %d degrees, trait=%s, NP=%d\n",
            basename(prefix), length(common), ndeg, trait, NP))
cat(sprintf("QC: pre=%d post=%d dropped=%d MAF>=%.3f callrate>%.2f\n",
            fq$qc$markers_pre_qc, fq$qc$markers_post_qc, fq$qc$markers_dropped_qc,
            fq$qc$maf_min, fq$qc$callrate_min))

fit_fr <- function(y) {                                  # point-estimate degree fractions
  fm <- BGLR(y = y, ETA = ETA(), nIter = 4000, burnIn = 1000, verbose = FALSE,
             saveAt = file.path(tempdir(), "perm_"))
  v <- sapply(en, function(nm) fm$ETA[[nm]]$varU); v / sum(v)
}

y0 <- as.numeric(phe[[trait]])
out <- data.frame()
fr0 <- fit_fr(y0)                                        # perm 0 = real data
out <- rbind(out, data.frame(perm = 0, trait = trait, degree = en, frac = round(fr0, 4)))
cat(sprintf("  real:  %s\n", paste(sprintf("%s=%.3f", en, fr0), collapse = " ")))
set.seed(1)
for (i in 1:NP) {
  yp <- sample(y0); fr <- fit_fr(yp)
  out <- rbind(out, data.frame(perm = i, trait = trait, degree = en, frac = round(fr, 4)))
}
write.csv(out, paste0(prefix, "_ortho_perm.csv"), row.names = FALSE)

# console summary: real vs permutation mean +/- sd per degree
wide <- reshape(out, idvar = "perm", timevar = "degree", drop = "trait", direction = "wide")
pm <- wide[wide$perm > 0, -1, drop = FALSE]
cat("  perm-null mean+/-sd per degree:\n")
for (k in seq_along(en)) cat(sprintf("    %-4s real=%.3f  null=%.3f +/- %.3f\n",
    en[k], fr0[k], mean(pm[[k]]), sd(pm[[k]])))
cat("PERM DONE\n")

#!/usr/bin/env Rscript
# =====================================================================
# Generic autotetraploid (ploidy=4) variance-partition + diagnostics runner for the
# cross-crop dataset extension. Reuses experiments/simulations/poly_ortho/lib_ortho.R.
# Produces {prefix}_ortho.csv, {prefix}_diagnostics.csv, {prefix}_null.csv
# (+ {prefix}_PA.csv if do_cv=1).
# Usage: Rscript run_ortho_dataset.R <dosage_rds> <pheno_rds> <out_prefix> [traits|auto] [do_cv:0/1]
#   dosage_rds: individuals x markers matrix, values 0..4, rownames = sample IDs
#   pheno_rds : data.frame, FIRST column = sample ID, remaining = numeric traits
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho.R")
qc_lib <- Sys.getenv("MARKER_QC_LIB", "experiments/empirical/multi/lib_marker_qc.R")
if (!file.exists(qc_lib) && file.exists("lib_marker_qc.R")) qc_lib <- "lib_marker_qc.R"
source(qc_lib)
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3)
dosage_rds <- args[1]; pheno_rds <- args[2]; prefix <- args[3]
traits_arg <- if (length(args) >= 4) args[4] else "auto"
do_cv <- if (length(args) >= 5) as.integer(args[5]) else 0L
set.seed(123)
OUTDIR <- dirname(prefix); if (nzchar(OUTDIR)) dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

M <- readRDS(dosage_rds)
phe <- readRDS(pheno_rds); if (is.matrix(phe)) phe <- as.data.frame(phe)
idcol <- colnames(phe)[1]; rownames(phe) <- as.character(phe[[idcol]])
common <- intersect(rownames(M), rownames(phe))
stopifnot(length(common) >= 50, all(M >= 0 & M <= 4, na.rm = TRUE))
M <- M[common, , drop = FALSE]; phe <- phe[common, ]
fq <- apply_marker_qc(M, ploidy = 4L)
M <- fq$M
write_marker_qc_audit(prefix, basename(prefix), fq$qc)
TR <- if (traits_arg == "auto") setdiff(colnames(phe), idcol) else strsplit(traits_arg, ",")[[1]]
drop <- grepl("group|gen[ae]tic|subpop|^pop$|cluster|^id$|clone|geno", TR, ignore.case = TRUE)
if (any(drop)) cat("dropping non-phenotype columns:", paste(TR[drop], collapse = ", "), "\n")
TR <- TR[!drop]
TR <- TR[sapply(TR, function(t) suppressWarnings(sd(as.numeric(phe[[t]]), na.rm = TRUE)) > 0)]
cat(sprintf("[%s] %d indiv, %d markers, traits: %s\n", basename(prefix), length(common), ncol(M), paste(TR, collapse = ", ")))
cat(sprintf("QC: pre=%d post=%d dropped=%d MAF>=%.3f callrate>%.2f\n",
            fq$qc$markers_pre_qc, fq$qc$markers_post_qc, fq$qc$markers_dropped_qc,
            fq$qc$maf_min, fq$qc$callrate_min))

K <- build_kernels(M)
g <- function(m, r, c) as.numeric(m[r, c])     # safe scalar extract
res <- data.frame()
for (tr in TR) {
  y <- as.numeric(phe[[tr]]); if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) next
  rs <- run_model(y, list(A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")), c("A", "D"))
  ro <- run_model(y, list(A = list(K = K$OA, model = "RKHS"), D2 = list(K = K$OD2, model = "RKHS"),
                          D3 = list(K = K$OD3, model = "RKHS"), D4 = list(K = K$OD4, model = "RKHS")), c("A", "D2", "D3", "D4"))
  s <- rs$ci; f <- ro$ci                        # rows mean/lo/hi ; cols kernels (PROPER posterior fraction CrIs)
  res <- rbind(res, data.frame(dataset = basename(prefix), trait = tr,
    std_Vd = round(g(s,"mean","D"),3), std_lo = round(g(s,"lo","D"),3), std_hi = round(g(s,"hi","D"),3),
    A = round(g(f,"mean","A"),3),
    D2 = round(g(f,"mean","D2"),3), D2_lo = round(g(f,"lo","D2"),3), D2_hi = round(g(f,"hi","D2"),3),
    D3 = round(g(f,"mean","D3"),3), D3_lo = round(g(f,"lo","D3"),3), D3_hi = round(g(f,"hi","D3"),3),
    D4 = round(g(f,"mean","D4"),3), D4_lo = round(g(f,"lo","D4"),3), D4_hi = round(g(f,"hi","D4"),3),
    nonadd = round(g(f,"mean","D2") + g(f,"mean","D3") + g(f,"mean","D4"), 3)))
  cat(sprintf("  %-16s stdVd=%.2f[%.2f,%.2f] | D2=%.2f D3=%.2f D4=%.2f\n",
              tr, g(s,"mean","D"), g(s,"lo","D"), g(s,"hi","D"), g(f,"mean","D2"), g(f,"mean","D3"), g(f,"mean","D4")))
}
write.csv(res, paste0(prefix, "_ortho.csv"), row.names = FALSE)

dd <- design_diag(list(K$OA, K$OD2, K$OD3, K$OD4))
write.csv(data.frame(dataset = basename(prefix), AD = fcos(K$A, K$D), OA_OD2 = fcos(K$OA, K$OD2),
  OD2_OD3 = fcos(K$OD2, K$OD3), OD3_OD4 = fcos(K$OD3, K$OD4), cond = dd$cond, eff_rank = dd$eff_rank),
  paste0(prefix, "_diagnostics.csv"), row.names = FALSE)
cat(sprintf("diagnostics: AD_cos=%.2f eff_rank=%.2f cond=%.1f\n", fcos(K$A, K$D), dd$eff_rank, dd$cond))

if (length(TR) >= 1) {            # phenotype-permutation null on the first trait (30 perms)
  y0 <- as.numeric(phe[[TR[1]]]); set.seed(1); P <- 30; sn <- numeric(P)
  for (i in 1:P) { yp <- sample(y0); fsx <- fit_std(yp, K, nIter = 3000, burnIn = 750); sn[i] <- vget(fsx, "D") / (vget(fsx, "A") + vget(fsx, "D")) }
  write.csv(data.frame(perm = 1:P, std_Vd = round(sn, 3)), paste0(prefix, "_null.csv"), row.names = FALSE)
  cat(sprintf("perm-null std Vd-share (%d perms): %.2f (sd %.2f)\n", P, mean(sn), sd(sn)))
}

if (do_cv == 1L) {
  pa <- data.frame()
  for (tr in TR) { y <- as.numeric(phe[[tr]]); if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) next
    a <- cv_pa(y, list(A = list(K = K$A, model = "RKHS")))
    s2 <- cv_pa(y, list(A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")))
    o <- cv_pa(y, list(A = list(K = K$OA, model = "RKHS"), D2 = list(K = K$OD2, model = "RKHS"), D3 = list(K = K$OD3, model = "RKHS"), D4 = list(K = K$OD4, model = "RKHS")))
    pa <- rbind(pa, data.frame(trait = tr, PA_additive = round(a, 3), PA_std_AD = round(s2, 3), PA_orthogonal = round(o, 3)))
  }
  write.csv(pa, paste0(prefix, "_PA.csv"), row.names = FALSE)
}
cat("DATASET DONE\n")

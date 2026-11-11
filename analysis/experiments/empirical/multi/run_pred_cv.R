#!/usr/bin/env Rscript
# =====================================================================
# Cross-crop test of the alfalfa PREDICTION conclusion (paper Sec 3.3):
#   "A+D helps phenotype prediction, but orthogonal A+D2+...+Dm does NOT beat A+D."
# Only ever run on Roza2019 alfalfa; the cross-crop expansion (exp:014-016) produced
# variance decompositions only (do_cv=0). This driver runs the missing prediction test.
#
# For each trait, paired K-fold x R-rep CV of THREE models (shared fold partitions):
#   (1) A only              additive baseline
#   (2) A + D               standard digenic dominance ("baseline to beat")
#   (3) A + D2 + ... + Dm    orthogonal full higher-order model
# Metric: Pearson r(predicted, observed) pooled across folds, per rep.
# Paired deltas (orth - AD, AD - A) reported with SD across reps.
# Handles ploidy 4, fixed m>4 (sweetpotato 6x), and per-marker variable ploidy.
#
# PARALLELISM: the unit of mclapply work is a single (trait, model, rep) CV pass,
# so cores are NOT capped by the trait count (matters for many-trait panels like
# miscanthus). Per-unit set.seed(1000+r) makes the fold partition deterministic and
# shared across the three models -> identical numerics to a per-trait sequential run.
#
# Usage: Rscript run_pred_cv.R <dosage_rds> <pheno_rds> <out_prefix> <ploidy> [traits|auto] [reps]
#   <ploidy>: 4 | 6 | auto | varploidy:<mploidy_rds>
# Output: <out_prefix>_PA.csv ; Env: CV_CORES (default 4). Pin BLAS to 1 thread.
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho_mp.R")   # also sources lib_ortho.R
qc_lib <- Sys.getenv("MARKER_QC_LIB", "experiments/empirical/multi/lib_marker_qc.R")
if (!file.exists(qc_lib) && file.exists("lib_marker_qc.R")) qc_lib <- "lib_marker_qc.R"
source(qc_lib)
suppressMessages(library(parallel))

args <- commandArgs(TRUE); stopifnot(length(args) >= 4)
dosage_rds <- args[1]; pheno_rds <- args[2]; prefix <- args[3]; ploidy_arg <- args[4]
traits_arg <- if (length(args) >= 5) args[5] else "auto"
REPS  <- if (length(args) >= 6) as.integer(args[6]) else 5L
CORES <- as.integer(Sys.getenv("CV_CORES", "4"))
DS <- basename(prefix)
OUTDIR <- dirname(prefix); if (nzchar(OUTDIR)) dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

M <- readRDS(dosage_rds); phe <- readRDS(pheno_rds); if (is.matrix(phe)) phe <- as.data.frame(phe)
idcol <- colnames(phe)[1]; rownames(phe) <- as.character(phe[[idcol]])
common <- intersect(rownames(M), rownames(phe)); stopifnot(length(common) >= 50)
M <- M[common, , drop = FALSE]; phe <- phe[common, , drop = FALSE]
TR <- if (traits_arg == "auto") setdiff(colnames(phe), idcol) else strsplit(traits_arg, ",")[[1]]
# drop non-phenotype structure/group columns that auto-selection can pick up
drop <- grepl("group|gen[ae]tic|subpop|^pop$|cluster", TR, ignore.case = TRUE)
if (any(drop)) cat("dropping non-phenotype columns:", paste(TR[drop], collapse = ", "), "\n")
TR <- TR[!drop]
TR <- TR[sapply(TR, function(t) suppressWarnings(sd(as.numeric(phe[[t]]), na.rm = TRUE)) > 0)]
nmap <- sapply(TR, function(t) sum(!is.na(as.numeric(phe[[t]]))))

# ---- build standard (A,D) + orthogonal degree kernels for the right ploidy regime ----
if (grepl("^varploidy:", ploidy_arg)) {
  mp <- readRDS(sub("^varploidy:", "", ploidy_arg)); stopifnot(length(mp) == ncol(M))
  fq <- apply_marker_qc(M, ploidy = NA_integer_, marker_ploidy = mp)
  M <- fq$M; mp <- fq$marker_ploidy
  maxdeg <- min(6L, min(mp)); K <- build_kernels_varm(M, mp, maxdeg); ndeg <- maxdeg; plab <- paste0("var<=", maxdeg)
} else {
  PL <- if (ploidy_arg == "auto") max(round(M), na.rm = TRUE) else as.integer(ploidy_arg)
  fq <- apply_marker_qc(M, ploidy = PL)
  M <- fq$M
  stopifnot(all(M >= 0 & M <= PL, na.rm = TRUE))
  if (PL == 4L) { Kb <- build_kernels(M); K <- list(A = Kb$A, D = Kb$D, O = list(Kb$OA, Kb$OD2, Kb$OD3, Kb$OD4)); ndeg <- 4L }
  else { K <- build_kernels_m(M, PL); ndeg <- PL }
  plab <- as.character(PL)
}
en <- ortho_etanames(ndeg)                         # c("A","D2",...,"Dndeg")
write_marker_qc_audit(prefix, DS, fq$qc)
rm(M); gc(verbose = FALSE)
mk_eta <- function(model) {                        # K-based ETA, identical to alfalfa pipeline
  if (model == "A")  return(list(A = list(K = K$A, model = "RKHS")))
  if (model == "AD") return(list(A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")))
  e <- list(); for (k in seq_len(ndeg)) e[[en[k]]] <- list(K = K$O[[k]], model = "RKHS"); e  # orth
}
cat(sprintf("[%s] ploidy=%s indiv=%d ndeg=%d reps=%d cores=%d traits(%d)\n",
            DS, plab, length(common), ndeg, REPS, CORES, length(TR)))
cat(sprintf("QC: pre=%d post=%d dropped=%d MAF>=%.3f callrate>%.2f\n",
            fq$qc$markers_pre_qc, fq$qc$markers_post_qc, fq$qc$markers_dropped_qc,
            fq$qc$maf_min, fq$qc$callrate_min))

# one (trait, model, rep) -> single pooled-fold PA value
run_unit <- function(ti, model, r, folds = 5, nIter = 1500, burnIn = 400) {
  tr <- TR[ti]; y <- as.numeric(phe[[tr]]); ok <- which(!is.na(y))
  set.seed(1000 + r); fold <- sample(rep(1:folds, length.out = length(ok)))
  ETA <- mk_eta(model); d <- tempfile("cvt"); dir.create(d); pr <- rep(NA_real_, length(y))
  for (f in 1:folds) {
    te <- ok[fold == f]; yN <- y; yN[te] <- NA
    fm <- BGLR(y = yN, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE,
               saveAt = file.path(d, sprintf("%s_%d_%d_", model, r, f))); pr[te] <- fm$yHat[te]
  }
  unlink(d, recursive = TRUE)
  data.frame(trait = tr, model = model, rep = r, pa = suppressWarnings(cor(pr[ok], y[ok])))
}

units <- expand.grid(ti = seq_along(TR), model = c("A", "AD", "orth"), r = 1:REPS, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(units)),
                function(i) tryCatch(run_unit(units$ti[i], units$model[i], units$r[i]),
                                     error = function(e) { cat("ERR", TR[units$ti[i]], units$model[i], conditionMessage(e), "\n"); NULL }),
                mc.cores = max(1, min(CORES, nrow(units))))
res <- do.call(rbind, Filter(is.data.frame, res))

# aggregate to per-trait mean/sd + paired deltas
pa <- data.frame()
for (tr in TR) {
  sub <- res[res$trait == tr, ]
  gv <- function(md) { v <- sub[sub$model == md, ]; v$pa[order(v$rep)] }
  A <- gv("A"); AD <- gv("AD"); O <- gv("orth"); L <- min(length(A), length(AD), length(O))
  if (L < 1) next
  A <- A[1:L]; AD <- AD[1:L]; O <- O[1:L]
  pa <- rbind(pa, data.frame(dataset = DS, ploidy = plab, trait = tr, n = nmap[[tr]], ndeg = ndeg,
    PA_A = round(mean(A), 4),  PA_A_sd = round(sd(A), 4),
    PA_AD = round(mean(AD), 4), PA_AD_sd = round(sd(AD), 4),
    PA_orth = round(mean(O), 4), PA_orth_sd = round(sd(O), 4),
    d_AD_vs_A = round(mean(AD - A), 4),    d_AD_vs_A_sd = round(sd(AD - A), 4),
    d_orth_vs_AD = round(mean(O - AD), 4), d_orth_vs_AD_sd = round(sd(O - AD), 4)))
}
write.csv(pa, paste0(prefix, "_PA.csv"), row.names = FALSE)
for (i in seq_len(nrow(pa))) cat(sprintf("  %-28s A=%.3f  A+D=%.3f  orth=%.3f  | dAD=%+.3f  dOrth=%+.3f\n",
    substr(pa$trait[i], 1, 28), pa$PA_A[i], pa$PA_AD[i], pa$PA_orth[i], pa$d_AD_vs_A[i], pa$d_orth_vs_AD[i]))
cat("PRED CV DONE\n")

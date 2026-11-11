#!/usr/bin/env Rscript
# =====================================================================
# Structure-controlled models: does the standard dominance variance survive
# correction for population (genomic PCs) and family/group structure?
# Reuses experiments/simulations/poly_ortho/lib_ortho.R (build_kernels, run_model -> proper
# sample-based fraction CrIs, fcos). PC control runs on every dataset; family/
# group control runs only when a group RDS is supplied.
# Usage: Rscript run_structure_dataset.R <dosage_rds> <pheno_rds> <out_prefix> <traits|auto> [group_rds|none]
#   group_rds: optional RDS, data.frame col1=sample ID, col2=group label (family/cohort/...)
# Outputs: {prefix}_structure_std.csv, _structure_orth.csv, _structure_cosines.csv, _structure_null.csv
# =====================================================================
source("experiments/simulations/poly_ortho/lib_ortho.R")
qc_lib <- Sys.getenv("MARKER_QC_LIB", "experiments/empirical/multi/lib_marker_qc.R")
if (!file.exists(qc_lib) && file.exists("lib_marker_qc.R")) qc_lib <- "lib_marker_qc.R"
source(qc_lib)
args <- commandArgs(TRUE); stopifnot(length(args) >= 4)
dosage_rds <- args[1]; pheno_rds <- args[2]; prefix <- args[3]; traits_arg <- args[4]
group_rds <- if (length(args) >= 5) args[5] else "none"
set.seed(123); OUTDIR <- dirname(prefix); if (nzchar(OUTDIR)) dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

M <- readRDS(dosage_rds); phe <- readRDS(pheno_rds); if (is.matrix(phe)) phe <- as.data.frame(phe)
idcol <- colnames(phe)[1]; rownames(phe) <- as.character(phe[[idcol]])
common <- intersect(rownames(M), rownames(phe)); stopifnot(length(common) >= 50, all(M >= 0 & M <= 4, na.rm = TRUE))
M <- M[common, , drop = FALSE]; phe <- phe[common, ]
fq <- apply_marker_qc(M, ploidy = 4L)
M <- fq$M
write_marker_qc_audit(prefix, basename(prefix), fq$qc)
TR <- if (traits_arg == "auto") setdiff(colnames(phe), idcol) else strsplit(traits_arg, ",")[[1]]
drop <- grepl("group|gen[ae]tic|subpop|^pop$|cluster|^id$|clone|geno", TR, ignore.case = TRUE)
if (any(drop)) cat("dropping non-phenotype columns:", paste(TR[drop], collapse = ", "), "\n")
TR <- TR[!drop]
TR <- TR[sapply(TR, function(t) suppressWarnings(sd(as.numeric(phe[[t]]), na.rm = TRUE)) > 0)]
cat(sprintf("[%s] %d indiv, %d markers, %d traits\n", basename(prefix), length(common), ncol(M), length(TR)))
cat(sprintf("QC: pre=%d post=%d dropped=%d MAF>=%.3f callrate>%.2f\n",
            fq$qc$markers_pre_qc, fq$qc$markers_post_qc, fq$qc$markers_dropped_qc,
            fq$qc$maf_min, fq$qc$callrate_min))

K <- build_kernels(M)
ev <- eigen(K$A, symmetric = TRUE); PC <- ev$vectors[, 1:min(10, ncol(ev$vectors)), drop = FALSE]; rownames(PC) <- rownames(M)

# optional family/group incidence kernel F
hasF <- FALSE; Fk <- NULL; cosAF <- NA; cosDF <- NA
if (group_rds != "none" && file.exists(group_rds)) {
  gdf <- readRDS(group_rds); if (is.matrix(gdf)) gdf <- as.data.frame(gdf); rownames(gdf) <- as.character(gdf[[1]])
  g <- gdf[common, 2]
  if (sum(!is.na(g)) > 0 && length(unique(g[!is.na(g)])) >= 2) {
    g[is.na(g)] <- "NA_grp"; Z <- model.matrix(~ 0 + factor(g)); Fk <- tcrossprod(Z); Fk <- Fk / mean(diag(Fk)) + diag(1e-6, nrow(Fk))
    dimnames(Fk) <- list(common, common); hasF <- TRUE; cosAF <- fcos(K$A, Fk); cosDF <- fcos(K$D, Fk)
    cat(sprintf("family kernel: %d groups, A-F cos=%.2f D-F cos=%.2f\n", length(unique(g)), cosAF, cosDF))
  }
}

eta_pc <- function(npc) if (npc > 0) list(PC = list(X = PC[, 1:npc, drop = FALSE], model = "FIXED")) else list()
std_share <- function(y, npc, useF) {
  ETA <- eta_pc(npc); if (useF) ETA$Fam <- list(K = Fk, model = "RKHS")
  ETA$A <- list(K = K$A, model = "RKHS"); ETA$D <- list(K = K$D, model = "RKHS")
  ci <- run_model(y, ETA, c("A", "D"))$ci; c(mean = as.numeric(ci["mean","D"]), lo = as.numeric(ci["lo","D"]), hi = as.numeric(ci["hi","D"]))
}
orth_frac <- function(y, npc, useF) {
  ETA <- eta_pc(npc); if (useF) ETA$Fam <- list(K = Fk, model = "RKHS")
  ETA$A <- list(K = K$OA, model = "RKHS"); ETA$D2 <- list(K = K$OD2, model = "RKHS"); ETA$D3 <- list(K = K$OD3, model = "RKHS"); ETA$D4 <- list(K = K$OD4, model = "RKHS")
  run_model(y, ETA, c("A", "D2", "D3", "D4"))$ci
}

std <- data.frame(); orth <- data.frame()
for (tr in TR) { y <- as.numeric(phe[[tr]]); if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) next
  ctrls <- list(c(0,0), c(5,0), c(10,0)); if (hasF) ctrls <- c(ctrls, list(c(0,1)))
  for (cc in ctrls) { npc <- cc[1]; useF <- cc[2] == 1
    lab <- if (useF) "F" else paste0("PC", npc)
    s <- std_share(y, npc, useF)
    std <- rbind(std, data.frame(dataset = basename(prefix), trait = tr, control = lab,
      std_Vd = round(s["mean"],3), std_lo = round(s["lo"],3), std_hi = round(s["hi"],3)))
    cat(sprintf("  %-14s %-4s stdVd=%.2f [%.2f,%.2f]\n", tr, lab, s["mean"], s["lo"], s["hi"]))
  }
  for (cc in (if (hasF) list(c(0,0), c(10,0), c(0,1)) else list(c(0,0), c(10,0)))) { npc <- cc[1]; useF <- cc[2] == 1
    lab <- if (useF) "F" else paste0("PC", npc); f <- orth_frac(y, npc, useF)
    orth <- rbind(orth, data.frame(dataset = basename(prefix), trait = tr, control = lab,
      A = round(as.numeric(f["mean","A"]),3), D2 = round(as.numeric(f["mean","D2"]),3),
      D3 = round(as.numeric(f["mean","D3"]),3), D4 = round(as.numeric(f["mean","D4"]),3),
      D2_lo = round(as.numeric(f["lo","D2"]),3), D2_hi = round(as.numeric(f["hi","D2"]),3)))
  }
}
write.csv(std, paste0(prefix, "_structure_std.csv"), row.names = FALSE)
write.csv(orth, paste0(prefix, "_structure_orth.csv"), row.names = FALSE)

# kernel cosines including PC kernel + family kernel
PCk <- tcrossprod(PC); PCk <- PCk / mean(diag(PCk)) + diag(1e-6, nrow(PCk))
write.csv(data.frame(dataset = basename(prefix),
  AD = fcos(K$A, K$D), A_PC = fcos(K$A, PCk), D_PC = fcos(K$D, PCk),
  D2D3 = fcos(K$OD2, K$OD3), D3D4 = fcos(K$OD3, K$OD4), A_F = cosAF, D_F = cosDF),
  paste0(prefix, "_structure_cosines.csv"), row.names = FALSE)

# permutation null under NO control vs PC10 control (does the ~0.5 default survive?)
y0 <- as.numeric(phe[[TR[1]]]); set.seed(1); P <- 30; n0 <- numeric(P); n10 <- numeric(P)
for (i in 1:P) { yp <- sample(y0)
  f0 <- BGLR(y = yp, ETA = list(A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")), nIter = 3000, burnIn = 750, verbose = FALSE, saveAt = file.path(tempdir(), "n0_"))
  n0[i] <- f0$ETA$D$varU / (f0$ETA$A$varU + f0$ETA$D$varU)
  f1 <- BGLR(y = yp, ETA = list(PC = list(X = PC[, 1:10, drop = FALSE], model = "FIXED"), A = list(K = K$A, model = "RKHS"), D = list(K = K$D, model = "RKHS")), nIter = 3000, burnIn = 750, verbose = FALSE, saveAt = file.path(tempdir(), "n10_"))
  n10[i] <- f1$ETA$D$varU / (f1$ETA$A$varU + f1$ETA$D$varU)
}
write.csv(data.frame(perm = 1:P, null_PC0 = round(n0, 3), null_PC10 = round(n10, 3)), paste0(prefix, "_structure_null.csv"), row.names = FALSE)
cat(sprintf("perm-null D-share: PC0=%.2f PC10=%.2f\n", mean(n0), mean(n10)))
cat("STRUCTURE DONE\n")

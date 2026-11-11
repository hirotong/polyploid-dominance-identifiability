#!/usr/bin/env Rscript
# =====================================================================
# build_poscontrol_panel.R
# -----------------------------------------------------------------------------
# Write the simulated positive control as a panel in the format read by
# load_panel() (<DATA_DIR>/poscontrol_dosage.rds, poscontrol_pheno.rds), so it
# can go through EVERY stage of the panel-conditioned workflow
# (run_panel_workflow.R: geom, real, perm, addctrl, d2spike, structspike, pred)
# under the same rule as the empirical panels.
#
# The simulation is the one in positive_control.R, with the same seed and the
# same order of random draws: independent Hardy-Weinberg loci (linkage
# equilibrium), n = 1500, m = 2500, 300 causal loci with effects on the
# digenic contrast only, H2 = 0.6.
#
# Usage: Rscript build_poscontrol_panel.R <LIB_ortho.R> <DATA_DIR>
# =====================================================================
ARG <- commandArgs(TRUE)
if (length(ARG) < 2) stop("usage: build_poscontrol_panel.R <LIB_ortho.R> <DATA_DIR>", call. = FALSE)
LIB <- ARG[1]; DATA_DIR <- ARG[2]
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
source(LIB)                         # ortho()
set.seed(2026)
PL <- 4L; n <- 1500L; m <- 2500L; ncausal <- 300L; H2 <- 0.6

## ---- independent HWE loci (exact linkage equilibrium) ----
p <- runif(m, 0.15, 0.85)
M <- matrix(rbinom(n * m, size = PL, prob = rep(p, each = n)), nrow = n)
rownames(M) <- paste0("S", seq_len(n)); colnames(M) <- paste0("m", seq_len(m))
keep <- apply(M, 2, function(g) { pp <- mean(g) / PL; pp > 0 && pp < 1 }); M <- M[, keep, drop = FALSE]

## ---- TRUE pure-digenic genotypic value (effects only on the degree-2 contrast phi2) ----
causal <- sample(ncol(M), min(ncausal, ncol(M)))
gv <- numeric(n)
for (j in causal) { g <- M[, j]; pj <- mean(g) / PL; Q <- ortho(pj, PL); gv <- gv + rnorm(1) * Q[g + 1, 3] }  # col 3 = phi2 (digenic)
gv <- scale(gv)[, 1]; ve <- var(gv) * (1 - H2) / H2; y <- gv + rnorm(n, 0, sqrt(ve))

phe <- data.frame(id = rownames(M), d2_sim = y, stringsAsFactors = FALSE)
saveRDS(M,   file.path(DATA_DIR, "poscontrol_dosage.rds"))
saveRDS(phe, file.path(DATA_DIR, "poscontrol_pheno.rds"))
cat(sprintf("poscontrol panel: n=%d markers=%d causal=%d H2=%.2f  sum(M)=%d  mean(y)=%.6f  var(y)=%.6f\n",
            nrow(M), ncol(M), length(causal), H2, sum(M), mean(y), var(y)))
cat(sprintf("written: %s/poscontrol_{dosage,pheno}.rds\n", DATA_DIR))

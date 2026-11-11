#!/usr/bin/env Rscript
# =====================================================================
# estimator_sensitivity.R -- do the fitted variance fractions depend on the
# prior and the estimator?  One panel per call, on the CURRENT panel inputs,
# marker QC, kernels and traits of the panel-conditioned workflow (the same
# load_panel() / build_panel() path as run_panel_workflow.R).
#
# For each of the first PWF_PERM_TRAITS analysed traits (the traits used for the
# permutation null) and for both models (standard A + D; orthogonal A + D2 + ...):
#   BGLR  default priors (R2 = 0.5, df0 = 5), three seeds  -> Monte Carlo spread
#   BGLR  R2 in {0.2, 0.8};  df0 in {1, 30}                -> prior sensitivity
#   REML  standard model: profile restricted likelihood over the share s_D, with
#         a 95% likelihood-ratio interval; the same on EST_NPERM permuted phenotypes
#   REML  orthogonal model: variance fractions by direct maximisation
# BGLR chains are the manuscript's primary setting (12,000 / 3,000 / thin 5) and
# fractions are posterior means of per-iteration fractions.
#
# Usage:  Rscript estimator_sensitivity.R <panel_key> <out_dir>
# Env:    DATA_DIR, LIB_PANEL, LIB_ORTHO, POLYDX_SRC, LIB_REML (as in the workflow),
#         PWF_SEED=2026, PWF_PERM_TRAITS=3, EST_NPERM=30, EST_REML_ORTHO=1,
#         CV_CORES / SLURM_CPUS_PER_TASK
# Output: <out_dir>/<panel_key>__estsens.csv        (long format)
#         <out_dir>/<panel_key>__remlprofile.csv    (profile log-likelihood of s_D)
# =====================================================================
LIB <- Sys.getenv("LIB_PANEL", "experiments/empirical/multi/lib_panel_workflow.R")
source(LIB)
source(Sys.getenv("LIB_REML", "experiments/empirical/multi/lib_reml.R"))
suppressMessages(library(parallel))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: estimator_sensitivity.R <panel_key> <out_dir>", call. = FALSE)
panel_key <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DATA_DIR <- Sys.getenv("DATA_DIR", "data/processed/panel_inputs")
SEED     <- as.integer(Sys.getenv("PWF_SEED", "2026"))
NTR      <- as.integer(Sys.getenv("PWF_PERM_TRAITS", "3"))
NPERM    <- as.integer(Sys.getenv("EST_NPERM", "30"))
DO_RORTH <- Sys.getenv("EST_REML_ORTHO", "1") == "1"
CORES    <- as.integer(Sys.getenv("CV_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
NITER <- 12000L; BURN <- 3000L; THIN <- 5L
if (Sys.getenv("EST_SMOKE", "0") == "1") { NITER <- 600L; BURN <- 200L; NPERM <- 3L }

REG <- panel_registry()
REG$poscontrol <- list(name = "Positive control (sim)", ploidy = "4", group = "none")
meta <- REG[[panel_key]]; if (is.null(meta)) stop(sprintf("unknown panel '%s'", panel_key), call. = FALSE)
PANEL <- meta$name
ploidy_arg <- if (grepl("%s", meta$ploidy)) sprintf(meta$ploidy, DATA_DIR) else meta$ploidy
P <- load_panel(file.path(DATA_DIR, paste0(panel_key, "_dosage.rds")),
                file.path(DATA_DIR, paste0(panel_key, "_pheno.rds")), ploidy_arg, "none", max_traits = 8L)
KP <- build_panel(P$M, P$pcfg, need_X = FALSE)
traits <- P$traits[seq_len(min(NTR, length(P$traits)))]
en <- KP$etanames
cat(sprintf("[%s] n=%d markers=%d ndeg=%d traits=%s cores=%d niter=%d nperm=%d\n",
            PANEL, P$n, P$markers, KP$ndeg, paste(traits, collapse = ","), CORES, NITER, NPERM))

frac_fn <- if (exists("frac_ci_samples_polydx")) frac_ci_samples_polydx else
  utils::getFromNamespace("frac_ci_samples_polydx", "polyDx")

# one BGLR fit -> posterior mean and 95% interval of each component's fraction
bglr_frac <- function(y, Klist, R2, df0) {
  d <- tempfile("bg"); dir.create(d); pfx <- file.path(d, "x_"); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  ETA <- lapply(Klist, function(K) list(K = K, model = "RKHS", df0 = df0))
  BGLR::BGLR(y = y, ETA = ETA, nIter = NITER, burnIn = BURN, thin = THIN, R2 = R2, verbose = FALSE, saveAt = pfx)
  frac_fn(pfx, names(Klist), BURN, THIN)
}
K_std  <- list(A = KP$A, D = KP$D)
K_orth <- stats::setNames(lapply(seq_len(KP$ndeg), function(k) KP$O[[k]]), en)

SETTINGS <- data.frame(
  setting = c("default_seed1", "default_seed2", "default_seed3", "R2=0.2", "R2=0.8", "df0=1", "df0=30"),
  R2  = c(0.5, 0.5, 0.5, 0.2, 0.8, 0.5, 0.5),
  df0 = c(5, 5, 5, 5, 5, 1, 30), stringsAsFactors = FALSE)

row <- function(trait, model, estimator, setting, component, value, lo = NA_real_, hi = NA_real_, extra = NA_real_)
  data.frame(panel = PANEL, panel_key = panel_key, n = P$n, trait = trait, model = model, estimator = estimator,
             setting = setting, component = component, value = value, lo = lo, hi = hi, extra = extra,
             stringsAsFactors = FALSE)

jobs <- list()
for (ti in seq_along(traits)) {
  for (si in seq_len(nrow(SETTINGS))) jobs[[length(jobs) + 1]] <- list(kind = "bglr", ti = ti, si = si)
  jobs[[length(jobs) + 1]] <- list(kind = "reml_std", ti = ti)
  if (DO_RORTH) jobs[[length(jobs) + 1]] <- list(kind = "reml_orth", ti = ti)
}

S_GRID <- seq(0, 1, by = 0.02)
run_job <- function(j) {
  tr <- traits[j$ti]; y <- numeric_trait(P$phe[[tr]]); ok <- is.finite(y)
  if (j$kind == "bglr") {
    st <- SETTINGS[j$si, ]; set.seed(SEED + 1000L * j$ti + j$si)
    cs <- bglr_frac(y, K_std, st$R2, st$df0); co <- bglr_frac(y, K_orth, st$R2, st$df0)
    out <- row(tr, "standard", "BGLR", st$setting, "s_D", cs["mean", "D"], cs["lo", "D"], cs["hi", "D"])
    for (cn in en) out <- rbind(out, row(tr, "orthogonal", "BGLR", st$setting, cn, co["mean", cn], co["lo", cn], co["hi", cn]))
    return(list(rows = out))
  }
  yo <- y[ok]; sub <- function(K) K[ok, ok, drop = FALSE]
  if (j$kind == "reml_std") {
    set.seed(SEED + 7000L + j$ti)
    Y <- cbind(yo, vapply(seq_len(NPERM), function(i) sample(yo), numeric(length(yo))))
    ll <- reml_profile_share(Y, sub(KP$A), sub(KP$D), S_GRID)
    r <- reml_share_summary(ll[, 1], S_GRID)
    pn <- apply(ll[, -1, drop = FALSE], 2, function(v) S_GRID[which.max(v)])
    out <- rbind(
      row(tr, "standard", "REML", "profile", "s_D", r["share"], r["lo"], r["hi"]),
      row(tr, "standard", "REML", "profile", "LR_vs_s0", r["lr_s0"]),
      row(tr, "standard", "REML", "profile", "LR_vs_s1", r["lr_s1"]),
      row(tr, "standard", "REML", sprintf("perm%d", NPERM), "s_D_perm_mean", mean(pn),
          as.numeric(quantile(pn, 0.025)), as.numeric(quantile(pn, 0.975)), extra = sd(pn)),
      row(tr, "standard", "REML", sprintf("perm%d", NPERM), "s_D_perm_at_boundary", mean(pn <= 0 | pn >= 1)))
    prof <- data.frame(panel = PANEL, panel_key = panel_key, trait = tr, s = S_GRID,
                       loglik_minus_max = ll[, 1] - max(ll[, 1]), stringsAsFactors = FALSE)
    return(list(rows = out, prof = prof))
  }
  if (j$kind == "reml_orth") {
    m <- reml_multi(yo, lapply(K_orth, sub))
    if (is.null(m)) return(list(rows = row(tr, "orthogonal", "REML", "optim", "failed", NA_real_)))
    out <- do.call(rbind, lapply(en, function(cn) row(tr, "orthogonal", "REML", "optim", cn, m$frac[[cn]])))
    out <- rbind(out, row(tr, "orthogonal", "REML", "optim", "h2", sum(m$var) / (sum(m$var) + m$var_e),
                          extra = m$convergence))
    return(list(rows = out))
  }
}

res <- mclapply(jobs, function(j) tryCatch(run_job(j), error = function(e) {
  cat(sprintf("JOB ERROR kind=%s trait=%s: %s\n", j$kind, traits[j$ti], conditionMessage(e))); NULL }),
  mc.cores = CORES, mc.preschedule = FALSE)
nfail <- sum(vapply(res, is.null, logical(1)))
res <- Filter(Negate(is.null), res)
rows <- do.call(rbind, lapply(res, `[[`, "rows"))
prof <- do.call(rbind, Filter(Negate(is.null), lapply(res, `[[`, "prof")))
write.csv(rows, file.path(out_dir, sprintf("%s__estsens.csv", panel_key)), row.names = FALSE)
if (!is.null(prof)) write.csv(prof, file.path(out_dir, sprintf("%s__remlprofile.csv", panel_key)), row.names = FALSE)
cat(sprintf("ESTSENS DONE (%d jobs, %d failed, %d rows)\n", length(jobs), nfail, nrow(rows)))
if (nfail > 0) quit(save = "no", status = 1)

#!/usr/bin/env Rscript
# =====================================================================
# make_workflow_scorecard.R
# -----------------------------------------------------------------------------
# Aggregate the per-(panel, stage) CSVs produced by run_panel_workflow.R into
# the canonical workflow outputs + the panel-level diagnostic scorecard + figures
# + a written summary.  Runs locally (ggplot2 only; no BGLR).
#
# Usage:
#   Rscript experiments/empirical/multi/make_workflow_scorecard.R [RAW_DIR] [OUT_DIR]
#     RAW_DIR  dir of <panel>__<stage>.csv
#              (default experiments/empirical/multi/out/panel_workflow/raw)
#     OUT_DIR  dir for canonical CSVs + figures
#              (default experiments/empirical/multi/out/panel_workflow)
#
# Canonical CSVs written to OUT_DIR:
#   kernel_diagnostics.csv  real_trait_variance_diagnostics.csv
#   permutation_null.csv  panel_null_summary.csv
#   additive_only_control.csv  d2_spikein_control.csv  d2_spikein_summary.csv
#   structure_spike_control.csv  structure_spike_summary.csv
#   prediction_diagnostics.csv  panel_prediction_summary.csv
#   workflow_scorecard.csv
# Figures to OUT_DIR/figures:
#   fig_workflow_scorecard.{pdf,png}  fig_panel_conditioned_controls.{pdf,png}
#   fig_structure_spike_controls.pdf  fig_d2_recoverability.pdf
# =====================================================================
args <- commandArgs(TRUE)
OUT_DEFAULT <- "experiments/empirical/multi/out/panel_workflow"
RAW <- if (length(args) >= 1) args[1] else file.path(OUT_DEFAULT, "raw")
OUT <- if (length(args) >= 2) args[2] else OUT_DEFAULT
FIG <- file.path(OUT, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---- decision thresholds (transparent; tunable) ----------------------------
THR <- list(
  struct_large = 0.10,   # |S_struct_real| considered structure-sensitive
  recover_rate = 0.80,   # D2 spike recovery rate to call the panel "recoverable"
  ci_wide      = 0.25,   # D2 fraction CrI width "too wide" above this
  null_margin  = 0.05,   # min delta_null margin for "exceeds null"
  pred_gain    = 0.01    # material mean CV gain for "prediction-useful"
)

# panel display order (BB abbreviation; tetraploids, then higher/variable ploidy)
PANEL_ORDER <- c("Alfalfa", "Potato-SolCAP", "Potato-SLU", "BB-UF-GP", "BB-UF-fruit",
                 "BB-UF-call*", "Panicum", "Rose", "Miscanthus", "Sweetpotato", "Sugarcane",
                 "Positive control (sim)")

rd_glob <- function(stage) {
  fs <- list.files(RAW, pattern = sprintf("__%s\\.csv$", stage), full.names = TRUE)
  fs <- fs[file.size(fs) > 1]
  if (!length(fs)) return(NULL)
  dfs <- lapply(fs, function(f) tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                                         error = function(e) NULL))
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
  if (!length(dfs)) return(NULL)
  cols <- Reduce(intersect, lapply(dfs, names))
  do.call(rbind, lapply(dfs, function(d) d[, cols, drop = FALSE]))
}
ord <- function(df) { if (!is.null(df) && "panel" %in% names(df)) df$panel <- factor(df$panel, levels = PANEL_ORDER); df }
W <- function(df, name) { if (!is.null(df)) write.csv(df, file.path(OUT, name), row.names = FALSE); invisible(df) }
med <- function(x) median(x, na.rm = TRUE)
q95 <- function(x) as.numeric(quantile(x, 0.95, names = FALSE, na.rm = TRUE))

# =====================================================================
# 1. load raw stage tables
# =====================================================================
geom   <- rd_glob("geom")
real   <- rd_glob("real")
perm   <- rd_glob("perm")
addc   <- rd_glob("addctrl")
d2s    <- rd_glob("d2spike")
strs   <- rd_glob("structspike")
pred   <- rd_glob("pred")
# P_D2 = digenic fraction of TOTAL fitted genetic variance (additive INCLUDED), matching
# the appendix / tab:crosscrop_cri definition. This equals the stored s_hat_D2_total from
# the same posterior fit; override the legacy additive-excluded P_D2 column so all
# downstream localization logic (real, additive-null and spike-in) uses one definition.
for (nm in c("real", "addc", "d2s")) {
  d <- get(nm)
  if (!is.null(d) && all(c("P_D2", "s_hat_D2_total") %in% names(d))) { d$P_D2 <- d$s_hat_D2_total; assign(nm, d) }
}
# Single source of truth for the real-trait orthogonal D2 localization metrics:
# reuse the cross-crop ortho-run fits (out/<panel>_ortho.csv) that feed
# tab:crosscrop_cri, so tab:workflow and tab:crosscrop_cri report IDENTICAL
# P_D2 / L_D2 / D2-CrI-width for the tetraploid panels (the two pipelines are
# independent BGLR fits of the same model, so their posterior means differ by
# Monte-Carlo error). 6-degree panels (sugarcane, sweetpotato) are not in
# tab:crosscrop_cri and their *_ortho.csv lacks per-degree D2, so they keep the
# panel-workflow values.
ORTHO_DIR <- dirname(OUT)
ortho_map <- c(Alfalfa = "alfalfa", "Potato-SolCAP" = "potato", "Potato-SLU" = "potato2",
               "BB-UF-GP" = "blueberry", "BB-UF-fruit" = "blueberry2", "BB-UF-call*" = "blueberry_gc",
               Panicum = "panicum", Rose = "rose", Miscanthus = "miscanthus")
ortho_loc <- list()
# The cross-crop tables (tab:crosscrop, tab:crosscrop_cri) and fig:crosscrop are all
# generated from the curated tables/cross_crop_summary.csv, which applies the
# per-panel representative-trait handling (e.g. alfalfa uses its yi_year slice, not
# the 3-trait mean). Use it as the SINGLE source so the workflow scorecard reports
# identical shared quantities. perm SD / percentile come from the matching *_null.csv.
ccs_path <- Sys.getenv("CROSS_CROP_SUMMARY", "tables/cross_crop_summary.csv")
ccs <- if (file.exists(ccs_path)) read.csv(ccs_path, check.names = FALSE, stringsAsFactors = FALSE) else NULL
for (pn in names(ortho_map)) {
  pre <- ortho_map[pn]
  row <- if (!is.null(ccs)) ccs[ccs$name == pre, , drop = FALSE] else NULL
  if (is.null(row) || !nrow(row)) next
  share <- row$std_Vd[1]
  rec <- list(P = row$ortho_D2[1],
              L = row$ortho_D2[1] - max(row$ortho_D3[1], row$ortho_D4[1]),
              W = row$ortho_D2_hi[1] - row$ortho_D2_lo[1],
              share = share,
              perm_mean = row$perm_null[1],
              cos = row$AD_cos[1], rank = row$eff_rank[1])
  # Unrounded localization values. cross_crop_summary.csv stores three decimals, and
  # differences of rounded means (L, W) or a second rounding to two decimals can be
  # off by one in the last digit. Recompute from the per-trait ortho file with the same
  # trait selection as make_figures.R (alfalfa: the representative yi_year trait).
  of <- file.path(ORTHO_DIR, paste0(pre, "_ortho.csv"))
  if (file.exists(of)) {
    o <- read.csv(of, check.names = FALSE, stringsAsFactors = FALSE)
    if (pre == "alfalfa" && "yi_year" %in% o$trait) o <- o[o$trait == "yi_year", , drop = FALSE]
    need <- c("std_Vd", "D2", "D3", "D4", "D2_lo", "D2_hi")
    if (nrow(o) && all(need %in% names(o))) {
      m <- function(cc) mean(o[[cc]], na.rm = TRUE)
      rec$P <- m("D2"); rec$L <- m("D2") - max(m("D3"), m("D4"))
      rec$W <- m("D2_hi") - m("D2_lo"); rec$share <- m("std_Vd")
    }
  }
  # perm SD + percentile from the representative-trait permutation distribution
  fn <- file.path(ORTHO_DIR, paste0(pre, "_null.csv"))
  if (file.exists(fn)) {
    nul <- read.csv(fn)
    if ("std_Vd" %in% names(nul)) {
      rec$perm_sd <- sd(nul$std_Vd, na.rm = TRUE)
      rec$perm_pctile <- mean(nul$std_Vd <= share, na.rm = TRUE)
    }
  }
  ortho_loc[[pn]] <- rec
}
panels <- if (!is.null(geom)) intersect(PANEL_ORDER, unique(geom$panel)) else PANEL_ORDER
cat(sprintf("loaded raw: geom=%s real=%s perm=%s addctrl=%s d2spike=%s structspike=%s pred=%s\n",
  !is.null(geom), !is.null(real), !is.null(perm), !is.null(addc), !is.null(d2s), !is.null(strs), !is.null(pred)))
cat("panels present:", paste(panels, collapse = ", "), "\n")

# =====================================================================
# 2. canonical CSVs (direct passthrough)
# =====================================================================
if (!is.null(geom)) {
  kd <- geom[, c("panel", "ploidy", "n", "markers", "AD_cosine", "effective_rank",
                 "condition_number", "max_abs_degree_cosine", "rank_ratio", "valid_markers")]
  W(ord(kd), "kernel_diagnostics.csv")
}
if (!is.null(real)) {
  rv_cols <- c("panel", "trait", "ploidy", "n", "markers", "s_hat_AD", "s_hat_D2_total",
               "P_D2", "L_D2", "D2_CrI_low", "D2_CrI_high", "D2_CrI_width", "AD_cosine", "effective_rank",
               "D2_largest", "s_hat_AD_raw_pt", "s_hat_AD_PC10", "s_hat_AD_family",
               "S_PC_real", "S_family_real", "S_struct_real")
  rv_cols <- intersect(rv_cols, names(real))
  W(ord(real[, rv_cols]), "real_trait_variance_diagnostics.csv")
}
if (!is.null(perm)) W(ord(perm), "permutation_null.csv")
if (!is.null(addc)) W(ord(addc), "additive_only_control.csv")
if (!is.null(d2s))  W(ord(d2s),  "d2_spikein_control.csv")
if (!is.null(strs)) W(ord(strs), "structure_spike_control.csv")
if (!is.null(pred)) W(ord(pred), "prediction_diagnostics.csv")

# =====================================================================
# 3. permutation-null summary (per panel/trait + panel-level)
# =====================================================================
panel_null_summary <- NULL
if (!is.null(perm) && !is.null(real)) {
  # The permutations use the point estimator share_pt() (ratio of posterior-mean variances,
  # short chain). The like-for-like observed value is the same estimator on the real
  # phenotype, s_hat_AD_raw_pt, so that is the minuend of delta_null. The posterior-mean
  # fraction from the long chain (s_hat_AD) is kept alongside as *_fraction.
  real_s <- real[, intersect(c("panel", "trait", "s_hat_AD", "s_hat_AD_raw_pt"), names(real))]
  has_pt <- "s_hat_AD_raw_pt" %in% names(real_s)
  rows <- list()
  for (pn in unique(perm$panel)) for (tr in unique(perm$trait[perm$panel == pn])) {
    sp <- perm$s_hat_AD_perm[perm$panel == pn & perm$trait == tr]
    rs <- real_s$s_hat_AD[real_s$panel == pn & real_s$trait == tr]
    sfrac <- if (length(rs)) rs[1] else NA_real_
    rp <- if (has_pt) real_s$s_hat_AD_raw_pt[real_s$panel == pn & real_s$trait == tr] else numeric(0)
    sreal <- if (length(rp) && is.finite(rp[1])) rp[1] else sfrac
    pmean <- mean(sp, na.rm = TRUE); psd <- sd(sp, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(panel = pn, trait_or_panel_summary = tr,
      s_hat_AD_real = sreal, perm_mean = pmean, perm_sd = psd,
      perm_q05 = as.numeric(quantile(sp, .05, names = FALSE, na.rm = TRUE)),
      perm_q50 = as.numeric(quantile(sp, .50, names = FALSE, na.rm = TRUE)),
      perm_q95 = as.numeric(quantile(sp, .95, names = FALSE, na.rm = TRUE)),
      delta_null = sreal - pmean, z_null = if (is.finite(psd) && psd > 0) (sreal - pmean) / psd else NA_real_,
      s_hat_AD_real_fraction = sfrac, delta_null_fraction = sfrac - pmean,
      stringsAsFactors = FALSE)
  }
  ptab <- do.call(rbind, rows)
  # panel-level rows
  plev <- list()
  for (pn in unique(ptab$panel)) {
    sub <- ptab[ptab$panel == pn, ]
    pq95 <- mean(sub$perm_q95, na.rm = TRUE)
    plev[[length(plev) + 1]] <- data.frame(panel = pn, trait_or_panel_summary = "PANEL_SUMMARY",
      s_hat_AD_real = mean(sub$s_hat_AD_real, na.rm = TRUE), perm_mean = mean(sub$perm_mean, na.rm = TRUE),
      perm_sd = mean(sub$perm_sd, na.rm = TRUE), perm_q05 = mean(sub$perm_q05, na.rm = TRUE),
      perm_q50 = mean(sub$perm_q50, na.rm = TRUE), perm_q95 = pq95,
      delta_null = mean(sub$delta_null, na.rm = TRUE), z_null = mean(sub$z_null, na.rm = TRUE),
      s_hat_AD_real_fraction = mean(sub$s_hat_AD_real_fraction, na.rm = TRUE),
      delta_null_fraction = mean(sub$delta_null_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  panel_null_summary <- rbind(ptab, do.call(rbind, plev))
  W(ord(panel_null_summary), "panel_null_summary.csv")
}

# =====================================================================
# 4. D2 spike-in summary (per panel x s)
# =====================================================================
d2_spikein_summary <- NULL
if (!is.null(d2s)) {
  rows <- list()
  for (pn in unique(d2s$panel)) for (s in sort(unique(d2s$true_s_D2))) {
    sub <- d2s[d2s$panel == pn & d2s$true_s_D2 == s, ]
    pnull <- if (!is.null(panel_null_summary)) panel_null_summary$perm_mean[
      panel_null_summary$panel == pn & panel_null_summary$trait_or_panel_summary == "PANEL_SUMMARY"] else NA_real_
    rows[[length(rows) + 1]] <- data.frame(panel = pn, true_s_D2 = s, n_rep = nrow(sub),
      D2_recovery_rate = mean(sub$D2_largest, na.rm = TRUE),
      median_P_D2 = med(sub$P_D2), median_L_D2 = med(sub$L_D2),
      median_D2_CrI_width = med(sub$D2_CrI_width), median_s_hat_AD = med(sub$s_hat_AD),
      median_s_hat_D2_total = med(sub$s_hat_D2_total),
      median_delta_vs_perm_null = if (length(pnull) && is.finite(pnull[1])) med(sub$s_hat_AD) - pnull[1] else NA_real_,
      stringsAsFactors = FALSE)
  }
  d2_spikein_summary <- do.call(rbind, rows)
  W(ord(d2_spikein_summary), "d2_spikein_summary.csv")
}

# =====================================================================
# 5. structure-spike summary (per panel x structure_type x lambda) + lambda=0 ruler
# =====================================================================
structure_spike_summary <- NULL
if (!is.null(strs)) {
  rows <- list()
  for (pn in unique(strs$panel)) for (st in unique(strs$structure_type)) for (lam in sort(unique(strs$lambda))) {
    sub <- strs[strs$panel == pn & strs$structure_type == st & strs$lambda == lam, ]
    if (!nrow(sub)) next
    rows[[length(rows) + 1]] <- data.frame(panel = pn, structure_type = st, lambda = lam, n_rep = nrow(sub),
      mean_S_struct = mean(sub$S_struct, na.rm = TRUE), median_S_struct = med(sub$S_struct),
      q95_S_struct = q95(sub$S_struct),
      mean_S_PC = mean(sub$S_PC, na.rm = TRUE), mean_S_family = mean(sub$S_family, na.rm = TRUE),
      mean_s_hat_AD_raw = mean(sub$s_hat_AD_raw, na.rm = TRUE),
      mean_s_hat_AD_PC10 = mean(sub$s_hat_AD_PC10, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  structure_spike_summary <- do.call(rbind, rows)
  W(ord(structure_spike_summary), "structure_spike_summary.csv")
}
# per-panel structure ruler q95 at lambda=0 (false-positive baseline), max over structure types
struct_ruler_q95 <- function(pn) {
  if (is.null(structure_spike_summary)) return(NA_real_)
  v <- structure_spike_summary$q95_S_struct[structure_spike_summary$panel == pn & structure_spike_summary$lambda == 0]
  if (length(v)) max(v, na.rm = TRUE) else NA_real_
}
struct_high_lambda <- function(pn) {
  if (is.null(structure_spike_summary)) return(NA_real_)
  v <- structure_spike_summary$median_S_struct[structure_spike_summary$panel == pn & structure_spike_summary$lambda == 1]
  if (length(v)) max(v, na.rm = TRUE) else NA_real_
}

# =====================================================================
# 6. prediction summary (per panel)
# =====================================================================
panel_prediction_summary <- NULL
if (!is.null(pred)) {
  rows <- list()
  for (pn in unique(pred$panel)) {
    sub <- pred[pred$panel == pn, ]
    dA <- sub$delta_AD[is.finite(sub$delta_AD)]; dO <- sub$delta_orth[is.finite(sub$delta_orth)]
    seA <- if (length(dA) > 1) sd(dA) / sqrt(length(dA)) else NA_real_
    seO <- if (length(dO) > 1) sd(dO) / sqrt(length(dO)) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(panel = pn, n_traits = length(unique(sub$trait)),
      mean_delta_AD = mean(dA), se_delta_AD = seA,
      ci_low_AD = mean(dA) - 1.96 * ifelse(is.na(seA), 0, seA), ci_high_AD = mean(dA) + 1.96 * ifelse(is.na(seA), 0, seA),
      mean_delta_orth = mean(dO), se_delta_orth = seO,
      ci_low_orth = mean(dO) - 1.96 * ifelse(is.na(seO), 0, seO), ci_high_orth = mean(dO) + 1.96 * ifelse(is.na(seO), 0, seO),
      stringsAsFactors = FALSE)
  }
  panel_prediction_summary <- do.call(rbind, rows)

  W(ord(panel_prediction_summary), "panel_prediction_summary.csv")
}

# =====================================================================
# 7. WORKFLOW SCORECARD  (combine + derived metrics + rule-based labels)
# =====================================================================
get1 <- function(df, pn, col, default = NA_real_) {
  if (is.null(df)) return(default); v <- df[[col]][df$panel == pn]; v <- v[is.finite(v)]
  if (length(v)) v[1] else default
}
mean_by <- function(df, pn, col) { if (is.null(df)) return(NA_real_); v <- df[[col]][df$panel == pn]; mean(v[is.finite(v)]) }

sc <- list()
for (pn in panels) {
  # --- real traits ---
  rsub <- if (!is.null(real)) real[real$panel == pn, ] else NULL
  s_hat_AD_real <- if (!is.null(rsub)) mean(rsub$s_hat_AD, na.rm = TRUE) else NA_real_
  s_hat_AD_range <- if (!is.null(rsub) && nrow(rsub)) diff(range(rsub$s_hat_AD, na.rm = TRUE)) else NA_real_
  P_D2_real  <- if (!is.null(rsub)) mean(rsub$P_D2, na.rm = TRUE) else NA_real_
  L_D2_real  <- if (!is.null(rsub)) mean(rsub$L_D2, na.rm = TRUE) else NA_real_
  CrI_w_real <- if (!is.null(rsub)) mean(rsub$D2_CrI_width, na.rm = TRUE) else NA_real_
  if (!is.null(ortho_loc[[pn]])) {        # single-source override from the cross-crop ortho run
    P_D2_real  <- ortho_loc[[pn]]$P
    L_D2_real  <- ortho_loc[[pn]]$L
    CrI_w_real <- ortho_loc[[pn]]$W
    s_hat_AD_real <- ortho_loc[[pn]]$share
  }
  S_struct_real <- if (!is.null(rsub)) mean(rsub$S_struct_real, na.rm = TRUE) else NA_real_
  AD_cos <- get1(geom, pn, "AD_cosine"); eff_rank <- get1(geom, pn, "effective_rank")
  if (!is.null(ortho_loc[[pn]])) {        # ortho-run kernel geometry (matches tab:crosscrop_cri)
    if (!is.null(ortho_loc[[pn]]$cos))  AD_cos   <- ortho_loc[[pn]]$cos
    if (!is.null(ortho_loc[[pn]]$rank)) eff_rank <- ortho_loc[[pn]]$rank
  }

  # --- perm null: observed and null shares are compared on the SAME traits (the
  #     traits the perm stage ran) with the SAME estimator (share_pt), both from the
  #     workflow's own fits, i.e. the PANEL_SUMMARY row of panel_null_summary. The displayed s_hat_AD_real above may
  #     average more traits (and, for tetraploids, comes from the cross-crop ortho run),
  #     so it is NOT the minuend of delta_null; s_hat_AD_null_traits is. ---
  ps_row <- if (!is.null(panel_null_summary)) panel_null_summary[panel_null_summary$trait_or_panel_summary == "PANEL_SUMMARY", ] else NULL
  pn_mean <- get1(ps_row, pn, "perm_mean")
  pn_sd   <- get1(ps_row, pn, "perm_sd")
  pn_q95  <- get1(ps_row, pn, "perm_q95")
  s_hat_AD_matched <- get1(ps_row, pn, "s_hat_AD_real")
  perm_pctile <- if (!is.null(perm)) { sp <- perm$s_hat_AD_perm[perm$panel == pn]; if (length(sp) && is.finite(s_hat_AD_matched)) mean(sp <= s_hat_AD_matched, na.rm = TRUE) else NA_real_ } else NA_real_
  delta_null <- s_hat_AD_matched - pn_mean
  # z is the standardized distance Delta_null / sd(perm); derive it from the SAME
  # delta_null so the two always share a sign.
  z_null  <- if (is.finite(pn_sd) && pn_sd > 0) delta_null / pn_sd else NA_real_

  # --- additive-only control ---
  asub <- if (!is.null(addc)) addc[addc$panel == pn, ] else NULL
  add_med_s    <- if (!is.null(asub)) med(asub$s_hat_AD) else NA_real_
  add_q95_s    <- if (!is.null(asub)) q95(asub$s_hat_AD) else NA_real_
  add_med_P    <- if (!is.null(asub)) med(asub$P_D2) else NA_real_
  add_med_L    <- if (!is.null(asub)) med(asub$L_D2) else NA_real_
  add_q95_L    <- if (!is.null(asub)) q95(asub$L_D2) else NA_real_

  # --- D2 spike-in at s = 0.50 ---
  d2_rec <- d2_pP <- d2_pL <- d2_w <- NA_real_
  if (!is.null(d2_spikein_summary)) {
    drow <- d2_spikein_summary[d2_spikein_summary$panel == pn & abs(d2_spikein_summary$true_s_D2 - 0.5) < 1e-6, ]
    if (nrow(drow)) { d2_rec <- drow$D2_recovery_rate[1]; d2_pP <- drow$median_P_D2[1]; d2_pL <- drow$median_L_D2[1]; d2_w <- drow$median_D2_CrI_width[1] }
  }

  # --- structure ruler / high-lambda expectation ---
  ruler_q95 <- struct_ruler_q95(pn); high_lambda_S <- struct_high_lambda(pn)

  # --- prediction ---
  dAD  <- if (!is.null(panel_prediction_summary)) get1(panel_prediction_summary, pn, "mean_delta_AD") else NA_real_
  dAD_low <- if (!is.null(panel_prediction_summary)) get1(panel_prediction_summary, pn, "ci_low_AD") else NA_real_
  dOrth<- if (!is.null(panel_prediction_summary)) get1(panel_prediction_summary, pn, "mean_delta_orth") else NA_real_

  # --- derived: resemblance of real localization to spike vs additive-null ---
  dist_P_spike <- abs(P_D2_real - d2_pP);  dist_P_null <- abs(P_D2_real - add_med_P)
  dist_L_spike <- abs(L_D2_real - d2_pL);  dist_L_null <- abs(L_D2_real - add_med_L)
  P_pctile_vs_null <- if (!is.null(asub) && is.finite(P_D2_real)) mean(asub$P_D2 <= P_D2_real, na.rm = TRUE) else NA_real_
  L_pctile_vs_null <- if (!is.null(asub) && is.finite(L_D2_real)) mean(asub$L_D2 <= L_D2_real, na.rm = TRUE) else NA_real_

  # --- rule-based label ---
  structure_sensitive <- (is.finite(S_struct_real) && S_struct_real >= THR$struct_large) ||
                         (is.finite(S_struct_real) && is.finite(ruler_q95) && S_struct_real > ruler_q95)
  exceeds_null <- is.finite(delta_null) && delta_null > THR$null_margin
  recoverable  <- is.finite(d2_rec) && d2_rec >= THR$recover_rate
  resembles_spike <- is.finite(dist_P_spike) && is.finite(dist_P_null) &&
                     (dist_P_spike < dist_P_null) && (is.finite(L_pctile_vs_null) && L_pctile_vs_null > 0.95)
  narrow_ci    <- is.finite(CrI_w_real) && CrI_w_real < THR$ci_wide
  d2_supported <- exceeds_null && recoverable && resembles_spike && narrow_ci && !structure_sensitive
  pred_gain    <- is.finite(dAD) && dAD >= THR$pred_gain && (is.na(dAD_low) || dAD_low > 0)
  unresolved_geom <- !recoverable                      # panel geometry cannot recover a planted D2

  label <- if (structure_sensitive) "structure-sensitive"
           else if (d2_supported) "D2-supported"
           else if (pred_gain) "prediction-useful"
           else "model-basis"
  geom_note <- if (unresolved_geom) "unresolved by panel geometry (D2 spike not recoverable)" else ""

  sc[[length(sc) + 1]] <- data.frame(
    panel = pn, AD_cosine = AD_cos, effective_rank = eff_rank,
    s_hat_AD_real = s_hat_AD_real, s_hat_AD_trait_range = s_hat_AD_range,
    P_D2_real = P_D2_real, L_D2_real = L_D2_real, D2_CrI_width_real = CrI_w_real,
    s_hat_AD_null_traits = s_hat_AD_matched, delta_null_fraction = get1(ps_row, pn, "delta_null_fraction"),
    perm_mean = pn_mean, perm_q95 = pn_q95, delta_null = delta_null, z_null = z_null,
    real_s_hat_AD_percentile_vs_null = perm_pctile,
    add_median_s_hat_AD = add_med_s, add_q95_s_hat_AD = add_q95_s, add_median_P_D2 = add_med_P, add_median_L_D2 = add_med_L,
    D2_recovery_rate_s050 = d2_rec, spike_median_P_D2_s050 = d2_pP, spike_median_L_D2_s050 = d2_pL, spike_median_D2_CrI_width_s050 = d2_w,
    S_struct_real = S_struct_real, structure_ruler_q95_lambda0 = ruler_q95, expected_S_struct_high_lambda = high_lambda_S,
    real_P_D2_percentile_vs_additive_null = P_pctile_vs_null, real_L_D2_percentile_vs_additive_null = L_pctile_vs_null,
    real_P_D2_distance_to_spike = dist_P_spike, real_L_D2_distance_to_spike = dist_L_spike,
    mean_delta_AD = dAD, mean_delta_AD_ci_low = dAD_low, mean_delta_orth = dOrth,
    exceeds_null = exceeds_null, D2_recoverable = recoverable, resembles_spike = resembles_spike,
    narrow_D2_CrI = narrow_ci, structure_sensitive = structure_sensitive, prediction_gain = pred_gain,
    suggested_label = label, geometry_note = geom_note,
    stringsAsFactors = FALSE)
}
scorecard <- do.call(rbind, sc)
num <- sapply(scorecard, is.numeric); scorecard[num] <- lapply(scorecard[num], function(x) round(x, 4))
W(ord(scorecard), "workflow_scorecard.csv")
cat("\n=== WORKFLOW SCORECARD ===\n")
print(scorecard[, c("panel", "s_hat_AD_real", "delta_null", "S_struct_real", "AD_cosine",
                    "D2_recovery_rate_s050", "P_D2_real", "mean_delta_AD", "suggested_label")], row.names = FALSE)

# =====================================================================
# 8. FIGURES
# =====================================================================
suppressMessages(library(ggplot2))
save_both <- function(p, name, w = 12, h = 8) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "in")
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 200)
}
GOOD <- "#1a9850"; MIDc <- "#fee08b"; BADc <- "#d73027"; NEUT <- "#d9d9d9"
BLUEp <- "#4575b4"; BLUE0 <- "#d1e5f0"; GREYn <- "#bdbdbd"
band <- function(v, good, mid) if (!is.finite(v)) NEUT else if (good(v)) GOOD else if (mid(v)) MIDc else BADc
# favourability colour for biological-D2 interpretation (green = favourable)
fav_color <- function(metric, v) switch(metric,
  delta_null            = band(v, function(x) x > 0.10, function(x) x > 0.05),
  S_struct_real         = band(v, function(x) x < 0.05, function(x) x < 0.10),
  AD_cosine             = band(v, function(x) x < 0.85, function(x) x < 0.92),
  effective_rank        = band(v, function(x) x > 3.0,  function(x) x > 2.0),
  D2_recovery_rate_s050 = band(v, function(x) x >= 0.80, function(x) x >= 0.50),
  P_D2_real             = band(v, function(x) x > 0.50, function(x) x > 0.34),
  L_D2_real             = band(v, function(x) x > 0.10, function(x) x > 0.00),
  D2_CrI_width_real     = band(v, function(x) x < 0.25, function(x) x < 0.40),
  mean_delta_AD         = if (!is.finite(v)) NEUT else if (v > 0.01) BLUEp else if (v >= 0) BLUE0 else GREYn,
  NEUT)
lab_color <- function(l) switch(l, "D2-supported" = GOOD, "structure-sensitive" = "#f46d43",
  "prediction-useful" = BLUEp, "#bdbdbd")

# ---- Figure 1: panel-level diagnostic scorecard ----------------------------
metrics <- c("s_hat_AD_real", "delta_null", "S_struct_real", "AD_cosine", "effective_rank",
             "D2_recovery_rate_s050", "P_D2_real", "L_D2_real", "D2_CrI_width_real", "mean_delta_AD")
mlab <- c(s_hat_AD_real = "s_AD", delta_null = "d_null", S_struct_real = "S_struct",
          AD_cosine = "A-D cos", effective_rank = "eff. rank", D2_recovery_rate_s050 = "D2 recov.",
          P_D2_real = "P_D2", L_D2_real = "L_D2", D2_CrI_width_real = "D2 CrI w.", mean_delta_AD = "d A+D pred")
sc_long <- do.call(rbind, lapply(seq_len(nrow(scorecard)), function(i) {
  do.call(rbind, lapply(metrics, function(m) data.frame(
    panel = scorecard$panel[i], metric = mlab[[m]],
    value = scorecard[[m]][i],
    fill = if (m == "s_hat_AD_real") NEUT else fav_color(m, scorecard[[m]][i]),
    txt = ifelse(is.finite(scorecard[[m]][i]), sprintf("%.2f", scorecard[[m]][i]), "NA"),
    stringsAsFactors = FALSE)))
}))
sc_long$panel  <- factor(sc_long$panel, levels = rev(PANEL_ORDER))
sc_long$metric <- factor(sc_long$metric, levels = mlab[metrics])
lab_df <- data.frame(panel = factor(scorecard$panel, levels = rev(PANEL_ORDER)),
                     metric = factor("label", levels = c(levels(sc_long$metric), "label")),
                     txt = scorecard$suggested_label, fill = vapply(scorecard$suggested_label, lab_color, ""),
                     stringsAsFactors = FALSE)
sc_long$metric <- factor(sc_long$metric, levels = c(levels(sc_long$metric), "label"))
fig1 <- ggplot(sc_long, aes(metric, panel)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt), size = 2.7, colour = "grey10") +
  geom_tile(data = lab_df, aes(fill = fill), colour = "white", linewidth = 0.6) +
  geom_text(data = lab_df, aes(label = txt), size = 2.2, colour = "grey10") +
  scale_fill_identity() +
  scale_x_discrete(position = "top", labels = c(mlab[metrics], label = "suggested label")) +
  labs(title = "Panel-conditioned diagnostic scorecard for dominance-labelled variance",
       subtitle = "green = favourable for biological D2 | yellow = borderline | red = unfavourable | blue = prediction utility (separate axis)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 0, size = 8),
        panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 12))
save_both(fig1, "fig_workflow_scorecard", w = 13, h = 7)

# ---- Figure 2: panel-conditioned reference distributions -------------------
build_ctrl_long <- function(metric) {
  out <- list()
  if (!is.null(perm) && metric == "s_hat_AD")
    out$pn <- data.frame(panel = perm$panel, source = "perm null", value = perm$s_hat_AD_perm)
  if (!is.null(addc)) out$ad <- data.frame(panel = addc$panel, source = "additive-only", value = addc[[metric]])
  if (!is.null(d2s)) { d <- d2s[abs(d2s$true_s_D2 - 0.5) < 1e-6, ]; out$d2 <- data.frame(panel = d$panel, source = "D2 spike (s=0.5)", value = d[[metric]]) }
  df <- do.call(rbind, out); if (is.null(df)) return(NULL); df$metric <- metric; df
}
ctrl_metrics <- intersect(c("s_hat_AD", "P_D2", "L_D2"), names(addc))
cl <- do.call(rbind, lapply(ctrl_metrics, build_ctrl_long))
real_pts <- NULL
if (!is.null(real)) {
  rl <- list()
  for (m in ctrl_metrics) { mr <- if (m == "s_hat_AD") "s_hat_AD" else m
    if (mr %in% names(real)) rl[[m]] <- data.frame(panel = real$panel, metric = m, value = real[[mr]]) }
  real_pts <- do.call(rbind, rl)
}
if (!is.null(cl)) {
  cl$panel <- factor(cl$panel, levels = PANEL_ORDER)
  cl$source <- factor(cl$source, levels = c("perm null", "additive-only", "D2 spike (s=0.5)"))
  if (!is.null(real_pts)) real_pts$panel <- factor(real_pts$panel, levels = PANEL_ORDER)
  fig2 <- ggplot(cl, aes(source, value, fill = source)) +
    geom_violin(scale = "width", alpha = 0.7, linewidth = 0.25, colour = "grey40") +
    geom_boxplot(width = 0.16, outlier.size = 0.3, alpha = 0.9, linewidth = 0.25) +
    { if (!is.null(real_pts)) geom_hline(data = aggregate(value ~ panel + metric, real_pts, mean, na.rm = TRUE),
        aes(yintercept = value), colour = "#d73027", linetype = 2, linewidth = 0.4) } +
    facet_grid(metric ~ panel, scales = "free_y") +
    scale_fill_manual(values = c("perm null" = "#bdbdbd", "additive-only" = "#74add1", "D2 spike (s=0.5)" = "#1a9850")) +
    labs(title = "Panel-conditioned reference distributions vs real-trait value (red dashed)",
         subtitle = "each panel's own genotype matrix generates its controls; D2 spike is panel-specific, not a universal ruler",
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6), legend.position = "bottom",
          strip.text.x = element_text(size = 7), plot.title = element_text(face = "bold"))
  save_both(fig2, "fig_panel_conditioned_controls", w = 15, h = 7)
}

# ---- Figure 3: structure ruler (S_struct vs lambda) ------------------------
if (!is.null(structure_spike_summary)) {
  ss <- structure_spike_summary; ss$panel <- factor(ss$panel, levels = PANEL_ORDER)
  real_S <- if (!is.null(real)) aggregate(S_struct_real ~ panel, real, mean, na.rm = TRUE) else NULL
  if (!is.null(real_S)) real_S$panel <- factor(real_S$panel, levels = PANEL_ORDER)
  fig3 <- ggplot(ss, aes(lambda, median_S_struct, colour = structure_type)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
    { if (!is.null(real_S)) geom_hline(data = real_S, aes(yintercept = S_struct_real), linetype = 2, colour = "#d73027", linewidth = 0.4) } +
    facet_wrap(~ panel, scales = "free_y") +
    scale_colour_manual(values = c(PC = "#4575b4", family = "#d73027")) +
    labs(title = "Panel-conditioned structure ruler: removable D-share vs injected structure fraction lambda",
         subtitle = "red dashed = real-trait observed structure sensitivity S_struct (mean over traits)",
         x = "injected structure fraction lambda", y = "S_struct (raw - controlled D-share)", colour = "structure") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  save_both(fig3, "fig_structure_spike_controls", w = 12, h = 8)
}

# ---- Figure 4: D2 recoverability vs true s ---------------------------------
if (!is.null(d2_spikein_summary)) {
  ds <- d2_spikein_summary; ds$panel <- factor(ds$panel, levels = PANEL_ORDER)
  d4 <- rbind(
    data.frame(panel = ds$panel, true_s = ds$true_s_D2, metric = "D2 recovery rate", value = ds$D2_recovery_rate),
    data.frame(panel = ds$panel, true_s = ds$true_s_D2, metric = "median P_D2", value = ds$median_P_D2))
  fig4 <- ggplot(d4, aes(true_s, value, colour = panel)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(title = "Panel-conditioned D2 recoverability vs planted digenic fraction s",
         x = "true digenic fraction s (spike-in)", y = NULL, colour = "panel") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold"))
  save_both(fig4, "fig_d2_recoverability", w = 11, h = 5)
}

# =====================================================================
# 9. WRITTEN SUMMARY  (auto-generated answers to the 7 workflow questions)
# =====================================================================
S <- scorecard
plist <- function(v) if (!length(v)) "_none_" else paste(v, collapse = ", ")
substantial <- S$panel[is.finite(S$s_hat_AD_real) & S$s_hat_AD_real >= 0.30]
exceed_null <- S$panel[S$exceeds_null %in% TRUE]
struct_sens <- S$panel[S$structure_sensitive %in% TRUE]
recoverable <- S$panel[S$D2_recoverable %in% TRUE]
not_recov   <- S$panel[is.finite(S$D2_recovery_rate_s050) & !(S$D2_recoverable %in% TRUE)]
resemble    <- S$panel[S$resembles_spike %in% TRUE]
pred_gain   <- S$panel[S$prediction_gain %in% TRUE]

fmt_row <- function(r) sprintf(
  "- **%s** -> _%s_%s  (s_AD=%.2f, delta_null=%+.2f, z=%.1f, D2 recov=%s, P_D2=%.2f, S_struct=%.2f, dA+D=%s)",
  r$panel, r$suggested_label, if (nzchar(r$geometry_note)) paste0("; ", r$geometry_note) else "",
  r$s_hat_AD_real, r$delta_null, ifelse(is.finite(r$z_null), r$z_null, NA),
  ifelse(is.finite(r$D2_recovery_rate_s050), sprintf("%.2f", r$D2_recovery_rate_s050), "NA"),
  r$P_D2_real, r$S_struct_real, ifelse(is.finite(r$mean_delta_AD), sprintf("%+.3f", r$mean_delta_AD), "NA"))

md <- c(
  "# Panel-conditioned diagnostic workflow - summary",
  "",
  sprintf("Generated from `workflow_scorecard.csv` over %d panels. All reference controls are conditioned on each panel's own genotype matrix (no universal n/M ruler). H2=0.6, n_perm/n_ref per run settings.", nrow(S)),
  "",
  "## Decision rule",
  "Biological digenic (D2) dominance is **not** claimed unless a panel passes BOTH:",
  "(i) panel-conditioned D2 recoverability (a planted s=0.5 D2 spike is the largest non-additive component), AND",
  "(ii) real-trait D2 localization above its own permutation-null / additive-only baseline,",
  "with structure control not removing the signal and a non-wide D2 credible interval. Prediction gain is reported separately and never used to claim biological dominance.",
  "",
  "## 1. Substantial fitted D-labelled share (s_hat_AD >= 0.30)",
  paste0("Panels: ", plist(as.character(substantial))),
  "",
  "## 2. Exceed their own permutation null (delta_null > margin AND real share above perm q95)",
  paste0("Panels: ", plist(as.character(exceed_null))),
  "",
  "## 3. Structure-sensitive (S_struct_real >= 0.10 or above panel structure-ruler q95)",
  paste0("Panels: ", plist(as.character(struct_sens))),
  "",
  "## 4. Can recover a panel-conditioned D2 spike-in (recovery rate >= 0.80 at s=0.5)",
  paste0("Recoverable: ", plist(as.character(recoverable))),
  paste0("NOT recoverable (geometry-limited): ", plist(as.character(not_recov))),
  "",
  "## 5. Real traits resembling the D2 spike-in controls (P_D2 closer to spike than to additive-null AND L_D2 above additive-null q95)",
  paste0("Panels: ", plist(as.character(resemble))),
  "",
  sprintf("## 6. Material prediction gain from standard A+D (mean delta_AD >= %.2f and CI lower bound > 0)", THR$pred_gain),
  paste0("Panels: ", plist(as.character(pred_gain))),
  "",
  "## 7. Recommended reporting language per panel",
  unlist(lapply(seq_len(nrow(S)), function(i) fmt_row(S[i, ]))),
  "",
  "## Notes",
  "- A panel that fails BOTH real-trait D2 localization and D2 spike recovery is **unresolved by panel geometry**, not evidence of absent biological dominance (sanity check 4).",
  "- A panel with material prediction gain but failed localization is **prediction-useful**, not biologically D2-supported (sanity check 5).",
  "- Where real s_hat_AD sits at the permutation null, the fitted D-labelled share is a model-basis artefact, not dominance (sanity check 6).",
  ""
)
writeLines(md, file.path(OUT, "WORKFLOW_SUMMARY.md"))
cat("\nwrote", file.path(OUT, "WORKFLOW_SUMMARY.md"), "\n")

cat("FIGURES + CSVs written to", OUT, "\n")
cat("SCORECARD DONE\n")

## Internal NULL-safe accessors mirroring make_workflow_scorecard.R's get1()/
## mean_by()/med()/q95() helpers, adapted to work panel-by-panel over the
## already-subsetted stage tables passed to workflow_scorecard().
sc_get1 <- function(df, pn, col, default = NA_real_) {
  if (is.null(df) || !(col %in% names(df))) return(default)
  v <- df[[col]][df$panel == pn]; v <- v[is.finite(v)]
  if (length(v)) v[1] else default
}

## P_D2 is overridden by s_hat_D2_total when both columns exist, matching
## make_workflow_scorecard.R's single definition of digenic localisation
## (fitted digenic fraction of TOTAL genetic variance, additive included).
sc_override_pd2 <- function(d) {
  if (!is.null(d) && all(c("P_D2", "s_hat_D2_total") %in% names(d))) d$P_D2 <- d$s_hat_D2_total
  d
}

#' Rule-based diagnostic scorecard for the panel-conditioned workflow
#'
#' Combines the stage tables produced by \code{\link{run_workflow_stage}}
#' (real-trait variance fits, permutation null, and the panel-conditioned
#' additive-only / D2 spike-in / structure-spike known-truth controls, plus
#' a paired prediction-CV comparison) into one row per panel: derived
#' distances/percentiles, boolean gates, a \code{suggested_label} and a
#' \code{geometry_note}. The gate definitions and derived-quantity formulas
#' are transcribed verbatim from the decision rule in
#' \code{make_workflow_scorecard.R} (section 7,
#' "WORKFLOW SCORECARD"), which is the script that produced the manuscript's
#' published \code{workflow_scorecard.csv} and figures.
#'
#' \strong{This is the generic decision rule only.} It does not reproduce
#' \code{make_workflow_scorecard.R}'s manuscript-specific single-sourcing:
#' that script additionally overrides the tetraploid panels' real-trait
#' \code{P_D2}/\code{L_D2}/\code{D2_CrI_width}/\code{s_hat_AD}/kernel-geometry
#' and permutation-null figures from the curated
#' \code{cross_crop_summary.csv} and the matching \code{*_null.csv}
#' files (its \code{ortho_loc} block), so that the workflow scorecard and the
#' manuscript's cross-crop table report identical shared quantities for those
#' panels; it also writes \code{WORKFLOW_SUMMARY.md} and the workflow
#' figures. None of that overriding, or the summary/figure generation, is
#' applied here. This generic rule does reproduce the manuscript's published
#' \code{suggested_label} on all 11 panels, but not necessarily every
#' intermediate numeric column for the tetraploid panels until the
#' \code{cross_crop_summary.csv} single-sourcing above is applied. Callers
#' who need numbers that match the manuscript's published
#' \code{workflow_scorecard.csv} exactly for those panels should run
#' \code{make_workflow_scorecard.R} itself; this function is for scoring
#' panels (including new ones) from their own stage tables directly.
#'
#' Gate definitions (\code{thr} is \code{thresholds}):
#' \preformatted{
#' structure_sensitive <- (is.finite(S_struct_real) && S_struct_real >= thr$struct_large) ||
#'                        (is.finite(S_struct_real) && is.finite(ruler_q95) && S_struct_real > ruler_q95)
#' exceeds_null        <- is.finite(delta_null) && delta_null > thr$null_margin
#' recoverable         <- is.finite(d2_rec) && d2_rec >= thr$recover_rate
#' resembles_spike     <- is.finite(dist_P_spike) && is.finite(dist_P_null) &&
#'                        (dist_P_spike < dist_P_null) &&
#'                        (is.finite(L_pctile_vs_null) && L_pctile_vs_null > thr$L_pctile)
#' narrow_ci           <- is.finite(CrI_w_real) && CrI_w_real < thr$ci_wide
#' d2_supported        <- exceeds_null && recoverable && resembles_spike && narrow_ci && !structure_sensitive
#' pred_gain           <- is.finite(dAD) && dAD >= thr$pred_gain &&
#'                        is.finite(dAD_low) && dAD_low > 0
#' suggested_label <- if (structure_sensitive) "structure-sensitive"
#'                    else if (d2_supported)   "D2-supported"
#'                    else if (pred_gain)      "prediction-useful"
#'                    else                     "model-basis"
#' geometry_note <- if (!recoverable) "unresolved by panel geometry (D2 spike not recoverable)" else ""
#' }
#' with \code{delta_null = s_hat_AD_real - perm_mean}, \code{z_null = delta_null / perm_sd},
#' \code{dist_P_spike = abs(P_D2_real - spike_median_P_D2_s050)},
#' \code{dist_P_null = abs(P_D2_real - add_median_P_D2)} and
#' \code{L_pctile_vs_null = mean(addctrl$L_D2 <= L_D2_real)} (all per panel).
#' \code{perm_mean}/\code{perm_sd}/\code{perm_q95} are themselves computed the
#' way \code{make_workflow_scorecard.R}'s \code{panel_null_summary}
#' \code{"PANEL_SUMMARY"} row is: \code{mean(sp)}/\code{sd(sp)}/
#' \code{quantile(sp, .95)} on each trait's permutation draws \code{sp}
#' separately, then averaged (unweighted) across the panel's traits -- not a
#' single summary of the panel's permutation draws pooled across traits (the
#' two agree for \code{perm_mean} when every trait has the same permutation
#' count, but not in general, and not for \code{perm_sd}/\code{perm_q95}/
#' \code{z_null} even then). \code{real_s_hat_AD_percentile_vs_null} is the
#' one permutation-derived quantity the script does pool directly across all
#' of a panel's traits.
#'
#' @param tables a named list of stage data frames, using any of the names
#'   \code{"geom"}, \code{"real"}, \code{"perm"}, \code{"addctrl"},
#'   \code{"d2spike"}, \code{"structspike"}, \code{"pred"} (the tables
#'   produced by the correspondingly-named \code{\link{run_workflow_stage}}
#'   stage, \code{rbind}-ed across panels). All are optional: a missing or
#'   \code{NULL} table yields \code{NA} metrics and \code{FALSE} gates for
#'   whatever it would have informed, never an error.
#' @param thresholds gate cut-offs; a list with \code{struct_large},
#'   \code{recover_rate}, \code{ci_wide}, \code{null_margin} and
#'   \code{pred_gain}. Defaults to \code{\link{workflow_thresholds}()}.
#' @param panels character vector of panel names to score. Defaults to
#'   \code{unique(tables$geom$panel)}; when \code{tables$geom} is absent it
#'   falls back to \code{unique(tables$real$panel)}.
#' @param spike_level the \code{true_s_D2} level of the D2 spike-in cell the
#'   recovery gate is evaluated at; must be present in \code{tables$d2spike}
#'   (warns when it is not). Defaults to 0.5, the manuscript's level.
#' @param ruler_lambda the \code{lambda} level of the structure spike-in used
#'   as the panel's false-positive ruler; warns when absent. Defaults to 0.
#' @param structure_ruler_q95 optional named numeric vector, panel name to
#'   the panel's structure-injection false-positive baseline at
#'   \code{lambda == 0} (see \code{\link{run_workflow_stage}}'s
#'   \code{"structspike"} stage). When \code{NULL} (the default) it is
#'   derived per panel from \code{tables$structspike}, matching the script's
#'   \code{struct_ruler_q95()}: \code{S_struct} rows at \code{lambda == 0}
#'   are first grouped by \code{structure_type} (e.g. \code{"PC"} vs
#'   \code{"family"}), a q95 is taken \emph{within} each group, and the
#'   ruler is the \emph{maximum} of those per-type q95s -- not a single q95
#'   pooled across structure types, which understates the ruler whenever a
#'   panel has more than one structure type (most real panels do).
#' @return a data frame with one row per panel and columns \code{panel},
#'   \code{AD_cosine}, \code{effective_rank}, \code{s_hat_AD_real},
#'   \code{s_hat_AD_trait_range}, \code{P_D2_real}, \code{L_D2_real},
#'   \code{D2_CrI_width_real}, \code{perm_mean}, \code{perm_q95},
#'   \code{delta_null}, \code{z_null},
#'   \code{real_s_hat_AD_percentile_vs_null}, \code{add_median_s_hat_AD},
#'   \code{add_q95_s_hat_AD}, \code{add_median_P_D2}, \code{add_median_L_D2},
#'   \code{D2_recovery_rate_s050}, \code{spike_median_P_D2_s050},
#'   \code{spike_median_L_D2_s050}, \code{spike_median_D2_CrI_width_s050},
#'   \code{S_struct_real}, \code{structure_ruler_q95_lambda0},
#'   \code{real_P_D2_percentile_vs_additive_null},
#'   \code{real_L_D2_percentile_vs_additive_null},
#'   \code{real_P_D2_distance_to_spike}, \code{real_L_D2_distance_to_spike},
#'   \code{mean_delta_AD}, \code{mean_delta_AD_ci_low}, \code{mean_delta_orth},
#'   \code{exceeds_null}, \code{D2_recoverable}, \code{resembles_spike},
#'   \code{narrow_D2_CrI}, \code{structure_sensitive}, \code{prediction_gain},
#'   \code{d2_supported}, \code{suggested_label} and \code{geometry_note}.
#'   Numeric columns are rounded to 4 decimals. When \code{panels} is empty
#'   (e.g. \code{tables} is an empty list and no \code{geom}/\code{real}
#'   table is available to infer panels from), returns a zero-row data frame
#'   with these same columns rather than \code{NULL}.
#'
#'   This column set differs from the published \code{workflow_scorecard.csv}
#'   in two ways: it \strong{adds} \code{d2_supported} (the boolean
#'   conjunction gate itself, not just the derived \code{suggested_label}),
#'   and it \strong{drops} \code{expected_S_struct_high_lambda}, a
#'   manuscript-specific diagnostic column that
#'   \code{make_workflow_scorecard.R} computes from its own structure-spike
#'   aggregation and appends after calling this generic gate logic; it is not
#'   part of the generic decision rule reproduced here.
#' @examples
#' geom <- data.frame(panel = "Demo", AD_cosine = 0.95, effective_rank = 1.6)
#' real <- data.frame(panel = "Demo", trait = "t1", s_hat_AD = 0.70,
#'                    s_hat_D2_total = 0.60, P_D2 = 0.60, L_D2 = 0.50,
#'                    D2_CrI_width = 0.10, S_struct_real = 0.01)
#' perm <- data.frame(panel = "Demo", trait = "t1", perm_id = seq_len(20),
#'                    s_hat_AD_perm = seq(0.46, 0.54, length.out = 20))
#' addctrl <- data.frame(panel = "Demo", rep = seq_len(20), s_hat_AD = 0.50,
#'                       s_hat_D2_total = 0.20, P_D2 = 0.20,
#'                       L_D2 = seq(-0.05, 0.05, length.out = 20),
#'                       D2_CrI_width = 0.3)
#' d2spike <- data.frame(panel = "Demo", rep = seq_len(20), true_s_D2 = 0.5,
#'                       s_hat_D2_total = 0.62, P_D2 = 0.62, L_D2 = 0.52,
#'                       D2_CrI_width = 0.2, D2_largest = TRUE)
#' sc <- workflow_scorecard(list(geom = geom, real = real, perm = perm,
#'                               addctrl = addctrl, d2spike = d2spike))
#' sc[, c("panel", "suggested_label", "d2_supported")]
#' @export
workflow_scorecard <- function(tables, thresholds = workflow_thresholds(),
                               panels = NULL, structure_ruler_q95 = NULL,
                               spike_level = 0.5, ruler_lambda = 0) {
  thr <- thresholds
  geom        <- tables$geom
  real        <- sc_override_pd2(tables$real)
  perm        <- tables$perm
  addctrl     <- sc_override_pd2(tables$addctrl)
  d2spike     <- sc_override_pd2(tables$d2spike)
  structspike <- tables$structspike
  pred        <- tables$pred

  if (is.null(panels)) {
    panels <- if (!is.null(geom)) unique(geom$panel)
              else if (!is.null(real)) unique(real$panel)
              else character(0)
  }

  if (!length(panels)) {
    return(data.frame(
      panel = character(0), AD_cosine = numeric(0), effective_rank = numeric(0),
      s_hat_AD_real = numeric(0), s_hat_AD_trait_range = numeric(0),
      P_D2_real = numeric(0), L_D2_real = numeric(0), D2_CrI_width_real = numeric(0),
      perm_mean = numeric(0), perm_q95 = numeric(0), delta_null = numeric(0), z_null = numeric(0),
      real_s_hat_AD_percentile_vs_null = numeric(0),
      add_median_s_hat_AD = numeric(0), add_q95_s_hat_AD = numeric(0),
      add_median_P_D2 = numeric(0), add_median_L_D2 = numeric(0),
      D2_recovery_rate_s050 = numeric(0), spike_median_P_D2_s050 = numeric(0),
      spike_median_L_D2_s050 = numeric(0), spike_median_D2_CrI_width_s050 = numeric(0),
      S_struct_real = numeric(0), structure_ruler_q95_lambda0 = numeric(0),
      real_P_D2_percentile_vs_additive_null = numeric(0),
      real_L_D2_percentile_vs_additive_null = numeric(0),
      real_P_D2_distance_to_spike = numeric(0), real_L_D2_distance_to_spike = numeric(0),
      mean_delta_AD = numeric(0), mean_delta_AD_ci_low = numeric(0), mean_delta_orth = numeric(0),
      exceeds_null = logical(0), D2_recoverable = logical(0), resembles_spike = logical(0),
      narrow_D2_CrI = logical(0), structure_sensitive = logical(0), prediction_gain = logical(0),
      d2_supported = logical(0), suggested_label = character(0), geometry_note = character(0),
      stringsAsFactors = FALSE))
  }

  med <- function(x) stats::median(x, na.rm = TRUE)
  q95 <- function(x) q_named(x, 0.95)

  rows <- vector("list", length(panels))
  for (i in seq_along(panels)) {
    pn <- panels[i]

    # --- real traits ---
    rsub <- if (!is.null(real)) real[real$panel == pn, , drop = FALSE] else NULL
    s_hat_AD_real <- if (!is.null(rsub)) mean(rsub$s_hat_AD, na.rm = TRUE) else NA_real_
    s_hat_AD_trait_range <- if (!is.null(rsub) && nrow(rsub)) diff(range(rsub$s_hat_AD, na.rm = TRUE)) else NA_real_
    P_D2_real  <- if (!is.null(rsub)) mean(rsub$P_D2, na.rm = TRUE) else NA_real_
    L_D2_real  <- if (!is.null(rsub)) mean(rsub$L_D2, na.rm = TRUE) else NA_real_
    CrI_w_real <- if (!is.null(rsub)) mean(rsub$D2_CrI_width, na.rm = TRUE) else NA_real_
    S_struct_real <- if (!is.null(rsub)) mean(rsub$S_struct_real, na.rm = TRUE) else NA_real_
    AD_cos   <- sc_get1(geom, pn, "AD_cosine")
    eff_rank <- sc_get1(geom, pn, "effective_rank")

    # --- permutation null ---
    # Panel-level perm_mean/perm_sd/perm_q95 are the *mean of per-trait*
    # summaries (mean(sp)/sd(sp)/quantile(sp,.95) computed separately for
    # each trait's permutation draws, then averaged across traits), exactly
    # as make_workflow_scorecard.R's panel_null_summary PANEL_SUMMARY row
    # does -- NOT a single summary of the pooled draws across all traits
    # (those differ once traits carry different numbers of permutations, or
    # for perm_sd/perm_q95 even when trait counts are equal).
    psub_panel <- if (!is.null(perm) && "s_hat_AD_perm" %in% names(perm)) perm[perm$panel == pn, , drop = FALSE] else NULL
    if (!is.null(psub_panel) && nrow(psub_panel) && "trait" %in% names(psub_panel)) {
      per_trait <- lapply(unique(psub_panel$trait), function(tr) {
        sp <- psub_panel$s_hat_AD_perm[psub_panel$trait == tr]
        c(mean = mean(sp, na.rm = TRUE), sd = stats::sd(sp, na.rm = TRUE), q95 = q95(sp))
      })
      per_trait <- do.call(rbind, per_trait)
      pn_mean <- mean(per_trait[, "mean"], na.rm = TRUE)
      pn_sd   <- mean(per_trait[, "sd"], na.rm = TRUE)
      pn_q95  <- mean(per_trait[, "q95"], na.rm = TRUE)
    } else if (!is.null(psub_panel) && nrow(psub_panel)) {
      # no trait column to group by: fall back to the pooled summary
      sp <- psub_panel$s_hat_AD_perm
      pn_mean <- mean(sp, na.rm = TRUE); pn_sd <- stats::sd(sp, na.rm = TRUE); pn_q95 <- q95(sp)
    } else {
      pn_mean <- NA_real_; pn_sd <- NA_real_; pn_q95 <- NA_real_
    }
    # real_s_hat_AD_percentile_vs_null pools all of the panel's permutation
    # draws directly (the script does not aggregate this one per trait first).
    psub <- if (!is.null(psub_panel)) psub_panel$s_hat_AD_perm else NULL
    perm_pctile <- if (!is.null(psub) && is.finite(s_hat_AD_real)) mean(psub <= s_hat_AD_real, na.rm = TRUE) else NA_real_
    delta_null <- s_hat_AD_real - pn_mean
    z_null <- if (is.finite(pn_sd) && pn_sd > 0) delta_null / pn_sd else NA_real_

    # --- additive-only control ---
    asub <- if (!is.null(addctrl)) addctrl[addctrl$panel == pn, , drop = FALSE] else NULL
    add_med_s <- if (!is.null(asub)) med(asub$s_hat_AD) else NA_real_
    add_q95_s <- if (!is.null(asub)) q95(asub$s_hat_AD) else NA_real_
    add_med_P <- if (!is.null(asub)) med(asub$P_D2) else NA_real_
    add_med_L <- if (!is.null(asub)) med(asub$L_D2) else NA_real_

    # --- D2 spike-in at true_s_D2 == 0.5 ---
    dsub <- if (!is.null(d2spike) && "true_s_D2" %in% names(d2spike)) {
      d2spike[d2spike$panel == pn & abs(d2spike$true_s_D2 - spike_level) < 1e-6, , drop = FALSE]
    } else NULL
    # `s_spike` is a caller-settable grid, so a legal run (e.g. s_spike = 0.6)
    # can contain no cell at `spike_level`. Silently that reads as "the spike
    # was not recovered" and the panel is labelled geometry-limited, which is
    # indistinguishable from every replicate having crashed.
    if (!is.null(d2spike) && any(d2spike$panel == pn) && (is.null(dsub) || !nrow(dsub))) {
      warning(sprintf("panel '%s': no d2spike rows at true_s_D2 == %g, so the recovery gate cannot be evaluated (available: %s)",
                      pn, spike_level,
                      paste(sort(unique(d2spike$true_s_D2[d2spike$panel == pn])), collapse = ", ")),
              call. = FALSE)
    }
    d2_rec <- if (!is.null(dsub) && nrow(dsub)) mean(dsub$D2_largest, na.rm = TRUE) else NA_real_
    d2_pP  <- if (!is.null(dsub) && nrow(dsub)) med(dsub$P_D2) else NA_real_
    d2_pL  <- if (!is.null(dsub) && nrow(dsub)) med(dsub$L_D2) else NA_real_
    d2_w   <- if (!is.null(dsub) && nrow(dsub)) med(dsub$D2_CrI_width) else NA_real_

    # --- structure ruler at lambda == 0 (false-positive baseline) ---
    # The script's struct_ruler_q95() reads this off structure_spike_summary,
    # which is aggregated per (panel, structure_type, lambda): q95(S_struct)
    # is computed WITHIN each structure_type separately, and the ruler is the
    # MAX of those per-type q95s at lambda == 0 -- not a single q95 pooled
    # across structure_type (PC vs family draws are not exchangeable).
    ruler_q95 <- if (!is.null(structure_ruler_q95)) {
      v <- unname(structure_ruler_q95[pn]); if (length(v)) as.numeric(v) else NA_real_
    } else if (!is.null(structspike) && all(c("S_struct", "lambda", "structure_type") %in% names(structspike))) {
      ssub0 <- structspike[structspike$panel == pn & structspike$lambda == ruler_lambda, , drop = FALSE]
      if (nrow(ssub0)) {
        per_type_q95 <- vapply(unique(ssub0$structure_type), function(st) q95(ssub0$S_struct[ssub0$structure_type == st]), numeric(1))
        if (length(per_type_q95)) max(per_type_q95, na.rm = TRUE) else NA_real_
      } else NA_real_
    } else if (!is.null(structspike) && all(c("S_struct", "lambda") %in% names(structspike))) {
      # no structure_type column to group by: fall back to the pooled q95
      ssub <- structspike$S_struct[structspike$panel == pn & structspike$lambda == ruler_lambda]
      if (length(ssub)) q95(ssub) else NA_real_
    } else NA_real_

    # `lambda` is likewise caller-settable; without a `ruler_lambda` cell the
    # ruler is NA, structure_sensitive falls back to the fixed threshold alone,
    # and "the control was never run" scores the same as "the control passed"
    # for the only vetoing term in d2_supported.
    if (!is.null(structspike) && is.null(structure_ruler_q95) &&
        "lambda" %in% names(structspike) && any(structspike$panel == pn) &&
        !any(structspike$panel == pn & structspike$lambda == ruler_lambda)) {
      warning(sprintf("panel '%s': no structspike rows at lambda == %g, so the panel-conditioned structure ruler is unavailable",
                      pn, ruler_lambda), call. = FALSE)
    }

    # --- prediction ---
    predsub <- if (!is.null(pred)) pred[pred$panel == pn, , drop = FALSE] else NULL
    dAD_v <- if (!is.null(predsub) && "delta_AD" %in% names(predsub)) predsub$delta_AD[is.finite(predsub$delta_AD)] else numeric(0)
    dOr_v <- if (!is.null(predsub) && "delta_orth" %in% names(predsub)) predsub$delta_orth[is.finite(predsub$delta_orth)] else numeric(0)
    dAD   <- if (length(dAD_v)) mean(dAD_v) else NA_real_
    dOrth <- if (length(dOr_v)) mean(dOr_v) else NA_real_
    seA   <- if (length(dAD_v) > 1) stats::sd(dAD_v) / sqrt(length(dAD_v)) else NA_real_
    # Substituting 0 for an unknown SE (a single CV rep) turned the lower
    # confidence bound into dAD itself, making the CI half of pred_gain a
    # tautology exactly when the evidence is weakest. One rep now yields NA,
    # which pred_gain treats as "not demonstrated".
    dAD_low <- if (is.finite(dAD) && is.finite(seA)) dAD - 1.96 * seA else NA_real_

    # --- resemblance of real localisation to spike vs additive-only null ---
    dist_P_spike <- abs(P_D2_real - d2_pP); dist_P_null <- abs(P_D2_real - add_med_P)
    dist_L_spike <- abs(L_D2_real - d2_pL)   # the L-distance to the null is not part of any gate
    P_pctile_vs_null <- if (!is.null(asub) && nrow(asub) && "P_D2" %in% names(asub) && is.finite(P_D2_real)) {
      mean(asub$P_D2 <= P_D2_real, na.rm = TRUE)
    } else NA_real_
    L_pctile_vs_null <- if (!is.null(asub) && nrow(asub) && "L_D2" %in% names(asub) && is.finite(L_D2_real)) {
      mean(asub$L_D2 <= L_D2_real, na.rm = TRUE)
    } else NA_real_

    # --- gates (transcribed verbatim from make_workflow_scorecard.R section 7) ---
    structure_sensitive <- (is.finite(S_struct_real) && S_struct_real >= thr$struct_large) ||
      (is.finite(S_struct_real) && is.finite(ruler_q95) && S_struct_real > ruler_q95)
    exceeds_null <- is.finite(delta_null) && delta_null > thr$null_margin
    recoverable  <- is.finite(d2_rec) && d2_rec >= thr$recover_rate
    resembles_spike <- is.finite(dist_P_spike) && is.finite(dist_P_null) &&
      (dist_P_spike < dist_P_null) &&
      (is.finite(L_pctile_vs_null) && L_pctile_vs_null > thr$L_pctile)
    narrow_ci <- is.finite(CrI_w_real) && CrI_w_real < thr$ci_wide
    d2_supported <- exceeds_null && recoverable && resembles_spike && narrow_ci && !structure_sensitive
    pred_gain <- is.finite(dAD) && dAD >= thr$pred_gain &&
      is.finite(dAD_low) && dAD_low > 0

    label <- if (structure_sensitive) "structure-sensitive"
             else if (d2_supported) "D2-supported"
             else if (pred_gain) "prediction-useful"
             else "model-basis"
    geometry_note <- if (!recoverable) "unresolved by panel geometry (D2 spike not recoverable)" else ""

    rows[[i]] <- data.frame(
      panel = pn, AD_cosine = AD_cos, effective_rank = eff_rank,
      s_hat_AD_real = s_hat_AD_real, s_hat_AD_trait_range = s_hat_AD_trait_range,
      P_D2_real = P_D2_real, L_D2_real = L_D2_real, D2_CrI_width_real = CrI_w_real,
      perm_mean = pn_mean, perm_q95 = pn_q95, delta_null = delta_null, z_null = z_null,
      real_s_hat_AD_percentile_vs_null = perm_pctile,
      add_median_s_hat_AD = add_med_s, add_q95_s_hat_AD = add_q95_s,
      add_median_P_D2 = add_med_P, add_median_L_D2 = add_med_L,
      D2_recovery_rate_s050 = d2_rec, spike_median_P_D2_s050 = d2_pP,
      spike_median_L_D2_s050 = d2_pL, spike_median_D2_CrI_width_s050 = d2_w,
      S_struct_real = S_struct_real, structure_ruler_q95_lambda0 = ruler_q95,
      real_P_D2_percentile_vs_additive_null = P_pctile_vs_null,
      real_L_D2_percentile_vs_additive_null = L_pctile_vs_null,
      real_P_D2_distance_to_spike = dist_P_spike, real_L_D2_distance_to_spike = dist_L_spike,
      mean_delta_AD = dAD, mean_delta_AD_ci_low = dAD_low, mean_delta_orth = dOrth,
      exceeds_null = exceeds_null, D2_recoverable = recoverable, resembles_spike = resembles_spike,
      narrow_D2_CrI = narrow_ci, structure_sensitive = structure_sensitive, prediction_gain = pred_gain,
      d2_supported = d2_supported, suggested_label = label, geometry_note = geometry_note,
      stringsAsFactors = FALSE)
  }

  out <- do.call(rbind, rows)
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(x) round(x, 4))
  out
}

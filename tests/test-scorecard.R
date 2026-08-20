library(polyDx)

thr <- workflow_thresholds()
stopifnot(abs(thr$struct_large - 0.10) < 1e-12, abs(thr$recover_rate - 0.80) < 1e-12,
          abs(thr$ci_wide - 0.25) < 1e-12, abs(thr$null_margin - 0.05) < 1e-12,
          abs(thr$pred_gain - 0.01) < 1e-12)

## --- fixture builder ------------------------------------------------------
mk <- function(panel, s_real, perm_mu, perm_sd, rec, P_real, L_real, ciw,
               add_P, add_L, spike_P, spike_L, S_struct, dAD, dAD_ok = TRUE) {
  list(
    geom = data.frame(panel = panel, AD_cosine = 0.95, effective_rank = 1.6),
    real = data.frame(panel = panel, trait = "t1", s_hat_AD = s_real,
                      s_hat_D2_total = P_real, P_D2 = P_real, L_D2 = L_real,
                      D2_CrI_width = ciw, S_struct_real = S_struct),
    perm = data.frame(panel = panel, trait = "t1", perm_id = seq_len(40),
                      s_hat_AD_perm = seq(perm_mu - 2 * perm_sd, perm_mu + 2 * perm_sd,
                                          length.out = 40)),
    addctrl = data.frame(panel = panel, rep = seq_len(40),
                         s_hat_AD = rep(perm_mu, 40),
                         s_hat_D2_total = rep(add_P, 40), P_D2 = rep(add_P, 40),
                         L_D2 = seq(add_L - 0.05, add_L + 0.05, length.out = 40),
                         D2_CrI_width = rep(0.3, 40)),
    d2spike = data.frame(panel = panel, rep = seq_len(40), true_s_D2 = 0.5,
                         s_hat_D2_total = rep(spike_P, 40), P_D2 = rep(spike_P, 40),
                         L_D2 = rep(spike_L, 40), D2_CrI_width = rep(0.2, 40),
                         D2_largest = rep(c(TRUE, FALSE), c(round(40 * rec), 40 - round(40 * rec)))),
    structspike = data.frame(panel = panel, rep = seq_len(10), structure_type = "PC",
                             lambda = 0, S_struct = rep(0.02, 10)),
    pred = data.frame(panel = panel, trait = "t1", rep = seq_len(10),
                      delta_AD = rep(dAD, 10) + if (dAD_ok) 0 else c(rep(0.5, 5), rep(-0.5, 5)),
                      delta_orth = rep(0, 10))
  )
}

## --- case 1: everything passes -> D2-supported ---------------------------
## real share far above the null, spike recoverable, real localization sits at
## the spike and above the additive-null L_D2 q95, narrow CrI, no structure.
t1 <- mk("Good", s_real = 0.70, perm_mu = 0.50, perm_sd = 0.02, rec = 1.0,
         P_real = 0.60, L_real = 0.50, ciw = 0.10,
         add_P = 0.20, add_L = 0.00, spike_P = 0.62, spike_L = 0.52,
         S_struct = 0.01, dAD = 0.03)
s1 <- workflow_scorecard(t1)
stopifnot(nrow(s1) == 1L, identical(s1$panel, "Good"))
stopifnot(s1$exceeds_null, s1$D2_recoverable, s1$resembles_spike, s1$narrow_D2_CrI)
stopifnot(!s1$structure_sensitive, s1$d2_supported)
stopifnot(identical(s1$suggested_label, "D2-supported"))
stopifnot(identical(s1$geometry_note, ""))
stopifnot(abs(s1$delta_null - (0.70 - 0.50)) < 1e-6)

## --- case 2: structure fires and wins the priority order -----------------
t2 <- t1; t2$real$S_struct_real <- 0.30
s2 <- workflow_scorecard(t2)
stopifnot(s2$structure_sensitive, !s2$d2_supported)
stopifnot(identical(s2$suggested_label, "structure-sensitive"))

## --- case 3: geometry cannot recover a planted D2 ------------------------
t3 <- mk("Flat", s_real = 0.70, perm_mu = 0.50, perm_sd = 0.02, rec = 0.20,
         P_real = 0.60, L_real = 0.50, ciw = 0.10,
         add_P = 0.20, add_L = 0.00, spike_P = 0.62, spike_L = 0.52,
         S_struct = 0.01, dAD = 0.03)
s3 <- workflow_scorecard(t3)
stopifnot(!s3$D2_recoverable, !s3$d2_supported)
stopifnot(identical(s3$suggested_label, "prediction-useful"))   # CV gain still fires
stopifnot(grepl("unresolved by panel geometry", s3$geometry_note))

## --- case 4: nothing fires -> model-basis --------------------------------
t4 <- mk("Null", s_real = 0.50, perm_mu = 0.50, perm_sd = 0.02, rec = 0.20,
         P_real = 0.20, L_real = 0.00, ciw = 0.40,
         add_P = 0.20, add_L = 0.00, spike_P = 0.62, spike_L = 0.52,
         S_struct = 0.01, dAD = 0.000)
s4 <- workflow_scorecard(t4)
stopifnot(!s4$exceeds_null, !s4$narrow_D2_CrI, !s4$prediction_gain)
stopifnot(identical(s4$suggested_label, "model-basis"))

## --- case 5: multiple panels in one call, column contract ----------------
tm <- Map(function(a, b) rbind(a, b), t1, t4)
sm <- workflow_scorecard(tm)
stopifnot(nrow(sm) == 2L, setequal(sm$panel, c("Good", "Null")))
need <- c("panel", "AD_cosine", "effective_rank", "s_hat_AD_real", "P_D2_real", "L_D2_real",
          "D2_CrI_width_real", "perm_mean", "perm_q95", "delta_null", "z_null",
          "D2_recovery_rate_s050", "S_struct_real", "structure_ruler_q95_lambda0",
          "mean_delta_AD", "mean_delta_orth", "exceeds_null", "D2_recoverable",
          "resembles_spike", "narrow_D2_CrI", "structure_sensitive", "prediction_gain",
          "d2_supported", "suggested_label", "geometry_note")
stopifnot(all(need %in% names(sm)))

## --- case 6: missing stage tables degrade to NA, not to an error ---------
s6 <- workflow_scorecard(list(geom = t1$geom, real = t1$real))
stopifnot(nrow(s6) == 1L, is.na(s6$D2_recovery_rate_s050), !s6$d2_supported)
stopifnot(identical(s6$suggested_label, "model-basis"))

## --- case 7: thresholds are honoured ------------------------------------
s7 <- workflow_scorecard(t1, thresholds = modifyList(thr, list(recover_rate = 1.5)))
stopifnot(!s7$D2_recoverable, !s7$d2_supported)

## --- case 8: subgroup aggregation, not pooling -----------------------------
## Regression coverage for two review findings: structure_ruler_q95_lambda0
## must be max(q95 per structure_type), and perm_mean/perm_sd/perm_q95 (and
## the z_null derived from them) must be the mean of per-(panel, trait)
## permutation summaries. Both sub-fixtures are built so the pooled-across-
## subgroup answer and the correct group-then-aggregate answer are
## independently computed here (via base quantile()/sd()/mean(), never by
## calling workflow_scorecard() itself) and differ materially, so a
## regression back to pooling fails loudly.

## -- 8a: structure ruler. PC is a large, tight cluster near 0.01; family is
## a single high rep near 0.09. Deliberately unbalanced replicate counts (39
## vs 1) so a pool-everything ruler is dominated by PC and dilutes family's
## high value out of the pooled top-5% entirely, while the correct
## per-type-q95-then-max answer is unaffected by how many reps the other
## type has.
pc_vals  <- rep(0.01, 39)
fam_vals <- c(0.09)
ruler_correct <- max(quantile(pc_vals, 0.95, names = FALSE), quantile(fam_vals, 0.95, names = FALSE))
ruler_pooled  <- quantile(c(pc_vals, fam_vals), 0.95, names = FALSE)
stopifnot(ruler_correct - ruler_pooled > 0.05)             # fixture actually discriminates
stopifnot(0.05 > ruler_pooled, 0.05 < ruler_correct)       # and straddles the real-value below

t8a <- list(
  geom = data.frame(panel = "RulerTest", AD_cosine = 0.90, effective_rank = 2.0),
  real = data.frame(panel = "RulerTest", trait = "t1", s_hat_AD = 0.50,
                    s_hat_D2_total = 0.20, P_D2 = 0.20, L_D2 = 0.00,
                    D2_CrI_width = 0.30, S_struct_real = 0.05),
  structspike = rbind(
    data.frame(panel = "RulerTest", rep = seq_along(pc_vals), structure_type = "PC",
               lambda = 0, S_struct = pc_vals),
    data.frame(panel = "RulerTest", rep = seq_along(fam_vals), structure_type = "family",
               lambda = 0, S_struct = fam_vals))
)
s8a <- workflow_scorecard(t8a)
stopifnot(abs(s8a$structure_ruler_q95_lambda0 - round(ruler_correct, 4)) < 1e-6)
## S_struct_real = 0.05 is below thr$struct_large (0.10), so only the ruler
## comparison decides this gate: it sits below the correct per-type ruler
## (0.09) but above the buggy pooled ruler (0.01) -- a pooled implementation
## would wrongly flag this panel as structure-sensitive.
stopifnot(!s8a$structure_sensitive)

## -- 8b: permutation null. Two traits with the SAME mean (0.5, so perm_mean
## cannot discriminate the two aggregations -- means are order/grouping
## invariant here) but very different spreads: tA tight (sd ~0.012), tB wide
## (sd ~0.125). perm_sd/perm_q95/z_null are NOT grouping-invariant.
tA <- seq(0.48, 0.52, length.out = 20)
tB <- seq(0.30, 0.70, length.out = 20)
perm_mean_correct <- mean(c(mean(tA), mean(tB)))
perm_sd_correct   <- mean(c(sd(tA), sd(tB)))
perm_q95_correct  <- mean(c(quantile(tA, 0.95, names = FALSE), quantile(tB, 0.95, names = FALSE)))
perm_sd_pooled    <- sd(c(tA, tB))
perm_q95_pooled   <- quantile(c(tA, tB), 0.95, names = FALSE)
stopifnot(abs(perm_sd_correct - perm_sd_pooled) > 0.01)   # fixture actually discriminates
stopifnot(abs(perm_q95_correct - perm_q95_pooled) > 0.03)

s_real <- 0.60
t8b <- list(
  geom = data.frame(panel = "PermTest", AD_cosine = 0.90, effective_rank = 2.0),
  real = data.frame(panel = "PermTest", trait = c("tA", "tB"), s_hat_AD = s_real,
                    s_hat_D2_total = 0.20, P_D2 = 0.20, L_D2 = 0.00,
                    D2_CrI_width = 0.30, S_struct_real = 0.01),
  perm = rbind(
    data.frame(panel = "PermTest", trait = "tA", perm_id = seq_along(tA), s_hat_AD_perm = tA),
    data.frame(panel = "PermTest", trait = "tB", perm_id = seq_along(tB), s_hat_AD_perm = tB))
)
s8b <- workflow_scorecard(t8b)
stopifnot(abs(s8b$perm_mean - round(perm_mean_correct, 4)) < 1e-6)
stopifnot(abs(s8b$perm_q95 - round(perm_q95_correct, 4)) < 1e-4)
z_null_correct <- (s_real - perm_mean_correct) / perm_sd_correct
z_null_pooled  <- (s_real - perm_mean_correct) / perm_sd_pooled
stopifnot(abs(z_null_correct - z_null_pooled) > 0.2)       # fixture actually discriminates
stopifnot(abs(s8b$z_null - round(z_null_correct, 4)) < 1e-3)

## --- case 9 (M10): empty input returns a zero-row data frame, not NULL ----
s9 <- workflow_scorecard(list())
stopifnot(is.data.frame(s9), nrow(s9) == 0L)
stopifnot(all(need %in% names(s9)))
stopifnot(identical(names(s9), names(s1)))

cat("test-scorecard.R OK\n")

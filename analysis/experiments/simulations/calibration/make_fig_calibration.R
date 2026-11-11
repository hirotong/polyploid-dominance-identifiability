#!/usr/bin/env Rscript
# =====================================================================
# fig_calibration — calibration curves for the standard A+D dominance share
# and the orthogonal D2 fraction vs the TRUE digenic fraction s.
#   Panel A  Standard A+D:  std_D_share = V_D/(V_A+V_D)  vs true s
#   Panel B  Orthogonal:    fitted D2 fraction           vs true s
# Both panels: one dot-and-I-beam series per population structure (mean +/- 1 SD
# over replicates), connected by lines, with the y=x identity (perfect
# calibration) reference. Style matches the redesigned Figure 1
# (ggplot2 + patchwork, theme_pub, structure by shape+colour, 300 dpi PNG+PDF).
#
# Also writes:
#   experiments/simulations/calibration/out/calibration_metrics.csv  (alpha,beta,RMSE,bias...)
#   paper/figures/figS_calibration_bias.{png,pdf}         (supplementary bias curves)
#
# Input : experiments/simulations/calibration/out/calibration_primary.csv
# Run from repo root:  Rscript experiments/simulations/calibration/make_fig_calibration.R
# =====================================================================
suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi) — single source of truth

IN <- "experiments/simulations/calibration/out/calibration_primary.csv"
OUT <- "paper/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
MET <- "experiments/simulations/calibration/out"
stopifnot(file.exists(IN))
d <- read.csv(IN, stringsAsFactors = FALSE)

## str_* encodings + theme_pub/theme_axes/save_fig come from experiments/fig_theme.R
## (scatter panels add theme_axes for the grid + axis lines).

## ---- long table of (estimate vs true s) for the two estimators ----------
std <- d |>
  filter(model == "A+D") |>
  transmute(
    structure = population_structure, rep = replicate_id, s = realized_true_D2_fraction,
    est = standard_D_share, estimator = "standard"
  )
orth <- d |>
  filter(model == "A+D2+D3+D4") |>
  transmute(
    structure = population_structure, rep = replicate_id, s = realized_true_D2_fraction,
    est = orth_D2_fraction, estimator = "orthogonal"
  )
L <- bind_rows(std, orth) |> filter(is.finite(est), is.finite(s))
L$structure <- factor(L$structure, str_lv)

## ---- calibration metrics: est ~ alpha + beta * s, per structure x estimator ----
metric_one <- function(g) {
  fit <- lm(est ~ s, data = g)
  b0 <- unname(coef(fit)[1])
  b1 <- unname(coef(fit)[2])
  bias <- g$est - g$s
  at <- function(v) {
    sel <- abs(g$s - v) < 1e-6
    if (any(sel)) mean(g$est[sel] - g$s[sel]) else NA_real_
  }
  data.frame(
    intercept_alpha = b0, slope_beta = b1,
    rmse_from_identity = sqrt(mean(bias^2)), mean_bias = mean(bias),
    bias_at_0 = at(0), bias_at_0.50 = at(0.50), bias_at_0.75 = at(0.75),
    n_obs = nrow(g)
  )
}
metrics <- L |>
  group_by(estimator, structure) |>
  group_modify(~ metric_one(.x)) |>
  ungroup() |>
  as.data.frame()
write.csv(metrics, file.path(MET, "calibration_metrics.csv"), row.names = FALSE)
cat("calibration metrics (est = alpha + beta * true_s):\n")
print(metrics, digits = 3)

## ---- per-cell summary (mean +/- SD) for plotting ------------------------
summ <- L |>
  group_by(estimator, structure, s) |>
  summarise(m = mean(est), sd = sd(est), .groups = "drop")
dx <- c(unrelated = -0.012, halfsib = 0, fullsib = 0.012) # tiny x-dodge so I-beams don't overlap
summ$xj <- summ$s + dx[as.character(summ$structure)]

## ---- pooled calibration regression (all structures) per estimator -------
pooled <- L |>
  group_by(estimator) |>
  group_modify(~ {
    f <- lm(est ~ s, data = .x)
    data.frame(alpha = unname(coef(f)[1]), beta = unname(coef(f)[2]))
  }) |>
  ungroup() |>
  as.data.frame()

## ---- fitted curve + 95% confidence band over the observed s range -------
## Confidence interval for the mean response (the fitted regression line). At
## n = 1800 this is a narrow band. fit/lwr/upr feed geom_line + geom_ribbon in panel().
pred <- L |>
  group_by(estimator) |>
  group_modify(~ {
    f <- lm(est ~ s, data = .x)
    sg <- seq(min(.x$s), max(.x$s), length.out = 200)
    ci <- as.data.frame(predict(f, newdata = data.frame(s = sg), interval = "confidence"))
    data.frame(s = sg, fit = ci$fit, lwr = ci$lwr, upr = ci$upr)
  }) |>
  ungroup() |>
  as.data.frame()

panel <- function(est_key, title, ylab, eq_xy, eq_lhs) {
  g <- summ |> filter(estimator == est_key)
  pe <- pooled[pooled$estimator == est_key, ]
  pr <- pred[pred$estimator == est_key, ]
  eqtxt <- sprintf("%s == %.2f + %.2f * italic(s)", eq_lhs, pe$alpha, pe$beta)
  ggplot(g, aes(xj, m, colour = structure, shape = structure, group = structure)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45) +
    # 95% confidence band for the pooled fitted regression (mean response)
    # geom_ribbon(
    #   data = pr, aes(x = s, ymin = lwr, ymax = upr), inherit.aes = FALSE,
    #   fill = "grey55", alpha = 0.55
    # ) +
    geom_line(data = pr, aes(x = s, y = fit), inherit.aes = FALSE, colour = "#542788", linewidth = 1.2) + # pooled fit
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.018, linewidth = 0.5) +
    geom_point(size = 2.5) +
    annotate("text",
      x = 0.7, y = 0.77, label = expression("perfect calibration " * hat(s) == s),
      colour = "grey55", size = 3, hjust = 1
    ) +
    annotate("text",
      x = eq_xy[1], y = eq_xy[2], label = eqtxt, parse = TRUE,
      colour = "#542788", size = 3, hjust = 0
    ) +
    scale_colour_manual(values = str_cols, labels = str_lab, name = paste0(leg_name, ":")) +
    scale_shape_manual(values = str_shp, labels = str_lab, name = paste0(leg_name, ":")) +
    scale_x_continuous(limits = c(-0.02, 0.8), breaks = c(0, 0.10, 0.25, 0.50, 0.75)) +
    scale_y_continuous(
      limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2),
      expand = expansion(mult = c(0, 0.03))
    ) +
    guides(shape = guide_legend(override.aes = list(size = 3.6))) +
    labs(
      title = title, x = expression("True digenic fraction  " * italic(s)),
      y = ylab
    ) +
    theme_pub +
    theme_axes +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.65))
}

yA <- expression("Fitted dominance-labelled share  " * italic(V)[D] / (italic(V)[A] + italic(V)[D]))
yB <- expression("Fitted orthogonal " * italic(D)[2] * " fraction")
pA <- panel("standard", expression("Standard " * italic(A) * " + " * italic(D) * ": dominance-labelled share"), yA,
  eq_xy = c(0.02, 0.64), eq_lhs = "hat(italic(s))[italic(D)]"
) +
  labs(tag = "A")
pB <- panel("orthogonal", expression("Orthogonal " * italic(D)[2] * " share"), yB,
  eq_xy = c(0.40, 0.10), eq_lhs = "hat(italic(s))[italic(D)[2]]"
) +
  labs(tag = "B")

# Shared horizontal legend collected below both panels:
#   "Population structure:  ● Unrelated   ▲ Half-sib   ■ Full-sib"
fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.title = element_text(size = 8.5, colour = "grey20", vjust = 0.5),
    legend.text = element_text(size = 8.5, colour = "grey20"),
    legend.key.size = unit(15, "pt"),
    legend.key.spacing.x = unit(7, "pt"),
    legend.margin = margin(4, 4, 2, 4)
  )

save_fig(fig, "fig_calibration", 200, 105)

## ---- supplementary: bias curves (est - s) -------------------------------
summ$bias <- summ$m - summ$s
pbias <- ggplot(summ, aes(xj, bias, colour = structure, shape = structure, group = structure)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45) +
  geom_line(aes(x = s), linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(~ factor(
    estimator, c("standard", "orthogonal"),
    c("Standard A+D share", "Orthogonal D2 fraction")
  )) +
  scale_colour_manual(values = str_cols, labels = str_lab, name = leg_name) +
  scale_shape_manual(values = str_shp, labels = str_lab, name = leg_name) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.25, 0.50, 0.75)) +
  labs(
    title = "Calibration bias: fitted share minus true digenic fraction",
    x = expression("True digenic fraction  " * italic(s)), y = expression(hat(s) - s)
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(hjust = 1),
    strip.text = element_text(size = 9.5, colour = "grey20")
  )
save_fig(pbias, "figS_calibration_bias", 200, 95)

cat("CALIBRATION FIG DONE\n")

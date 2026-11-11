#!/usr/bin/env Rscript
# Four-panel supplementary calibration figure for the SECONDARY stress grid.
# Shows BOTH parts of calibration -- the false baseline AND the response to truth:
#   A  no-dominance baseline (fitted standard D share at true s = 0)  heatmap n x M
#   B  calibration slope beta (response of the fitted share to true s) heatmap n x M
#   C  calibration RMSE from identity vs additive-digenic Frobenius cosine (scatter)
#   D  calibration RMSE from identity vs four-kernel effective rank      (scatter)
# A is hot (high = large spurious baseline = bad); B is a separate blue-green ramp
# (high = slope near 1 = good) so the two heatmaps are not read on the same scale.
# Style matches make_fig_simulation.R / make_fig_calibration_secondary.R (theme_pub).
#   in : experiments/simulations/calibration/out/calibration_secondary_metrics.csv
# Run from repo root: Rscript experiments/simulations/calibration/make_fig_calibration_4panel.R
suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi)

OUT <- "paper/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cells <- read.csv("experiments/simulations/calibration/out/calibration_secondary_metrics.csv", stringsAsFactors = FALSE)

nlev <- c(400, 1200, 4800, 19200)
nlab <- c("400", "1,200", "4,800", "19,200")
ncol4 <- c("#fdae61", "#74add1", "#4575b4", "#313695") # sequential by n (scatters)
Mlab2 <- function(M) paste0(M / 1000, "k") # heatmap x-ticks: 5k / 20k / 50k / 100k / 200k (axis title already says M)
Mvals <- c(5000, 20000, 50000, 100000, 200000) # 5-col set (M=1000 dropped: ill-conditioned A/D split at n=19200)
Mfac <- function(M) factor(Mlab2(M), Mlab2(Mvals))
full <- expand.grid(N = nlev, nMarker = Mvals)

## ---------- helper: n x M heatmap of a per-cell quantity (mean over LD) ----------
heat <- function(var, pal, dir, lims, hi_thresh, title, tag) {
  h <- aggregate(stats::reformulate("N + nMarker", var), cells, function(z) mean(z, na.rm = TRUE))
  h <- merge(full, h, all.x = TRUE)
  h$Nf <- factor(h$N, nlev, nlab)
  h$Mf <- Mfac(h$nMarker)
  h$val <- h[[var]]
  h$lab <- ifelse(is.finite(h$val), sprintf("%.2f", h$val), "n/a")
  ggplot(h, aes(Mf, Nf)) +
    geom_tile(aes(fill = val), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = lab, colour = is.finite(val) & val > hi_thresh), size = 3) +
    scale_fill_distiller(
      palette = pal, direction = dir, limits = lims,
      na.value = "grey88", name = NULL
    ) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
    labs(title = title, tag = tag, x = "Marker number  M", y = "Sample size  n") +
    theme_pub +
    theme(
      panel.grid = element_blank(), axis.ticks = element_blank(),
      legend.key.height = unit(20, "pt"), legend.key.width = unit(7, "pt")
    )
}

## A: no-dominance baseline  (fitted standard share at true s = 0); hot = bad
pA <- heat(
  "std_bias0", "YlOrRd", 1, c(0, 0.5), 0.28,
  expression("No-dominance baseline fitted " * hat(italic(s))[italic(A) + italic(D)] * " at true " * italic(s) * " = 0"), "A"
)

## B: calibration slope beta  (response of fitted share to true s); dark = good
pB <- heat(
  "std_beta", "GnBu", 1, c(0, 1), 0.55,
  expression("Calibration slope  " * hat(beta) * "  (response to true " * italic(s) * ")"), "B"
)

## ---------- C, D: calibration RMSE from identity vs kernel geometry ----------
cg <- cells[is.finite(cells$std_rmse) & cells$nMarker %in% Mvals, ]
cg$Nf <- factor(cg$N, nlev, nlab)
rC <- sprintf("%.2f", cor(cg$std_rmse, cg$AD_cos))
rD <- sprintf("%.2f", cor(cg$std_rmse, cg$eff_rank))
base_sc <- function(xv, xlab, lab, tag, ttl, ann_x) {
  ggplot(cg, aes(.data[[xv]], std_rmse, colour = Nf)) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey60", linewidth = 0.5, linetype = 2) +
    geom_point(aes(shape = factor(nMarker)), size = 2) +
    scale_colour_manual(values = setNames(ncol4, nlab), name = "Sample size  n") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8), labels = c("5k", "20k", "50k", "100k", "200k"), name = "Markers  M") +
    annotate("text", x = ann_x, y = 0.31, label = lab, size = 3.1, colour = "grey20", hjust = 0) +
    labs(title = ttl, tag = tag, x = xlab, y = "Calibration RMSE from identity") +
    theme_pub +
    theme_axes
}
pC <- base_sc(
  "AD_cos", "Additive–digenic Frobenius cosine", paste0("r = ", rC), "C",
  "Error vs A–D collinearity", 0.70
) + theme(legend.position = "none")
pD <- base_sc(
  "eff_rank", "Four-kernel effective rank", paste0("r = ", rD), "D",
  "Error vs effective rank", 3.05 
) +
  theme(
    legend.position = "right", legend.key.size = unit(9, "pt"),
    legend.title = element_text(size = 7), legend.text = element_text(size = 6.5)
  )

## ===================== assemble 2x2 + write ================================
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))
save_fig(fig, "fig_calibration_geometry", 210, 170)
cat("4-PANEL CALIBRATION FIG DONE\n")

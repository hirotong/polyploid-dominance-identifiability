#!/usr/bin/env Rscript
# =====================================================================
# Figure 1 (known-truth simulation) — single, self-contained generator.
# Side-by-side, publication-styled.
#   Panel A  Standard A+D model: dot-and-linerange of the fitted
#            dominance-labelled share V_D/(V_A+V_D) by simulated
#            architecture (x position only). Population structure
#            (unrelated / half-sib / full-sib) encoded by shape + colour.
#   Panel B  Orthogonal degree model: heatmap (full-sib) of the fitted
#            variance fraction per degree (A,D2,D3,D4) x architecture,
#            YlOrRd colourmap (0-0.7); full-height colourbar, title on
#            the right rotated 90 deg.
# Inputs (produced by poly_ortho_framework.R):
#   experiments/simulations/poly_ortho/out/std_vd_share.csv      (structure,truth,rep,Vd_share)
#   experiments/simulations/poly_ortho/out/ortho_attribution.csv (structure,truth,rep,fA..fD4)
# Output: paper/figures/fig_simulation.png (300 dpi) + paper/figures/fig_simulation.pdf
# Run from repo root:  Rscript experiments/simulations/poly_ortho/make_fig_simulation.R
# =====================================================================
suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi) — single source of truth

OUT <- "paper/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sv <- read.csv("experiments/simulations/poly_ortho/out/std_vd_share.csv")
oa <- read.csv("experiments/simulations/poly_ortho/out/ortho_attribution.csv")

## ---- shared encodings ---------------------------------------------------
arch_lv <- c("additive", "digenic", "trigenic", "quadrigenic", "hetero")
arch_lab <- c("Additive", "Digenic", "Trigenic", "Quadrigenic", "Heterozygosity")
## str_lv/str_lab/str_shp/str_cols/leg_name come from experiments/fig_theme.R
y_lab <- expression("Fitted dominance share  " * italic(V)[D] / (italic(V)[A] + italic(V)[D]))

## theme_pub / theme_axes / save_fig come from experiments/fig_theme.R

## ===================== PANEL A: standard A+D ============================
A <- sv |>
  mutate(truth = factor(truth, arch_lv), structure = factor(structure, str_lv)) |>
  group_by(structure, truth) |>
  summarise(m = mean(Vd_share), s = sd(Vd_share), .groups = "drop")

pd <- position_dodge(width = 0.6)
pA <- ggplot(A, aes(truth, m, group = structure, colour = structure, shape = structure)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), colour = "grey88", linewidth = 0.3) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), position = pd, width = 0.16, linewidth = 0.55) +
  geom_point(position = pd, size = 2.6) +
  scale_colour_manual(values = str_cols, labels = str_lab, name = leg_name) +
  scale_shape_manual(values = str_shp, labels = str_lab, name = leg_name) +
  scale_x_discrete(labels = arch_lab) +
  scale_y_continuous(y_lab,
    limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = expression("Standard " * italic(A) * " + " * italic(D) * ": inflated dominance-labelled share"),
    tag = "A", x = "Simulated architecture"
  ) +
  theme_pub +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    axis.line = element_line(colour = "grey40", linewidth = 0.4),
    axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
    axis.text.x = element_text(angle = 28, hjust = 0.5, vjust = 0.65),
    legend.position = "inside", legend.position.inside = c(0.985, 0.03),
    legend.justification.inside = c(1, 0),
    legend.title = element_text(size = 7.8, colour = "grey20"),
    legend.text = element_text(size = 7.3, colour = "grey20"),
    legend.key.size = unit(10, "pt"),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = "grey85", linewidth = 0.2),
    legend.margin = margin(2, 4, 3, 4)
  )

## ===================== PANEL B: orthogonal heatmap (full-sib) ===========
comp_lv <- c("A", "D2", "D3", "D4")
comp_par <- c("italic(A)", "italic(D)[2]", "italic(D)[3]", "italic(D)[4]")

fs <- oa |>
  filter(structure == "fullsib") |>
  mutate(truth = factor(truth, arch_lv)) |>
  group_by(truth) |>
  summarise(A = mean(fA), D2 = mean(fD2), D3 = mean(fD3), D4 = mean(fD4), .groups = "drop")
B <- data.frame(
  truth = factor(rep(fs$truth, times = 4), arch_lv),
  comp  = factor(rep(comp_lv, each = nrow(fs)), comp_lv),
  frac  = c(fs$A, fs$D2, fs$D3, fs$D4)
)
B <- B |>
  group_by(truth) |>
  mutate(colmax = frac == max(frac)) |>
  ungroup() |>
  as.data.frame()

pB <- ggplot(B, aes(truth, comp)) +
  geom_tile(aes(fill = frac)) + # no gap between cells
  geom_tile(data = B[B$colmax & B$truth != "hetero", ], fill = NA, colour = "black", linewidth = 0.3) + # box per-column max (single-degree truths only)
  geom_text(aes(label = sprintf("%.2f", frac), colour = frac > 0.45), size = 3.2) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 0.7), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  scale_x_discrete(labels = arch_lab) +
  scale_y_discrete(
    limits = rev(comp_lv),
    labels = function(l) parse(text = comp_par[match(l, comp_lv)])
  ) +
  coord_fixed(ratio = 1) + # strictly square cells (1 x-unit = 1 y-unit)
  labs(
    title = expression("Orthogonal model: degree-specific recovery"), tag = "B",
    x = "Simulated architecture", y = "Fitted component"
  ) +
  theme_pub +
  theme(
    plot.title = element_text(size = 10, hjust = 0, colour = "grey20", margin = margin(b = 5, l = 0)),
    panel.grid = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 0.5, vjust = 0.7),
    legend.position = "none"
  )

## ---- standalone colourbar: vertical YlOrRd bar with a thin frame, tick
## marks + labels on the right, and a horizontal two-line "Variance fraction"
## title sitting on top. The bar spans values 0-0.7 over a y-range that leaves
## headroom for the title; bar x-extent [bx0,bx1] fixes its (slim) width.
ybr <- seq(0, 0.7, 0.1)
barcols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(256)
barmat <- matrix(rev(barcols), ncol = 1) # row 1 = top = high value (dark red)
bx0 <- 0.
bx1 <- 0.8 # bar x-extent (width)
cbar <- ggplot() +
  annotation_raster(barmat, xmin = bx0, xmax = bx1, ymin = 0, ymax = 0.7, interpolate = TRUE) +
  annotate("rect",
    xmin = bx0, xmax = bx1, ymin = 0, ymax = 0.7,
    fill = NA, colour = "grey55", linewidth = 0.3
  ) +
  geom_segment(
    data = data.frame(y = ybr), aes(x = bx1, xend = bx1 + 0.33, y = y, yend = y),
    colour = "grey25", linewidth = 0.2
  ) +
  geom_text(
    data = data.frame(y = ybr), aes(x = bx1 + 0.55, y = y, label = sprintf("%.2f", y)),
    hjust = 0, size = 2.2, colour = "grey20"
  ) +
  annotate("text",
    x = bx0, y = 0.8, label = "Variance\nfraction",
    hjust = 0, vjust = 0, lineheight = 0.95, size = 2.6, colour = "grey20"
  ) +
  coord_cartesian(xlim = c(0, 6), ylim = c(-0.1, 0.98), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

## ===================== assemble + write ================================
fig <- (pA | pB | cbar) +
  plot_layout(widths = unit(c(1, 1, 2.2), c("null", "null", "cm"))) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

save_fig(fig, "fig_simulation", 224, 110)

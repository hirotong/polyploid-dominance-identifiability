#!/usr/bin/env Rscript
# =====================================================================
# fig_crosscrop — cross-crop diagnostic pattern (horizontal dot/range plots).
#   Crops on the y axis; three calm diagnostic panels read left to right.
#   Panel A  Standard dominance-labelled share V_D/(V_A+V_D): observed (black
#            circle), trait range (grey line), phenotype-permutation null
#            (red diamond). Shows the observed share is not obviously separated
#            from the no-signal permutation baseline.
#   Panel B  Additive-digenic Frobenius cosine (black circle); vertical
#            reference at 1.0 = "identical direction". Kernel-collinearity
#            diagnostic: the D share is uninterpretable if A and D point the
#            same empirical direction.
#   Panel C  Effective rank of the A + D2 + D3 + D4 kernel design (black
#            circle); dashed reference at the nominal maximum of 4.
# Red is used ONLY for the permutation null, so the rest stays calm. Crops are
# identified by the y-axis labels; no per-crop colour is used.
# Input : tables/cross_crop_summary.csv  (built by experiments/empirical/multi/make_figures.R)
#         panel_null_summary.csv (make_workflow_scorecard.R) for panel A: observed
#         and null shares on the SAME traits with the SAME estimator
#         (env PANEL_NULL_SUMMARY overrides the path)
# Output: <FIG_OUT or paper/figures>/fig_crosscrop.png (1200 dpi) + .pdf
# Run from repo root:  Rscript experiments/empirical/multi/make_fig_crosscrop.R
# =====================================================================
suppressMessages({
  library(ggplot2)
  library(patchwork)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi)


S <- read.csv("tables/cross_crop_summary.csv", check.names = FALSE, stringsAsFactors = FALSE)
# Display labels + top-to-bottom row order, keyed by the stable `name` column.
# Order matches tab:crosscrop / tab:datasets_supp (alphabetical by crop, then dataset).
lab_map <- c(
  alfalfa = "Alfalfa", blueberry_gc = "BB-UF-call", blueberry2 = "BB-UF-fruit",
  blueberry = "BB-UF-GP", miscanthus = "Miscanthus", panicum = "Panicum",
  potato2 = "Potato-SLU", potato = "Potato-SolCAP", rose = "Rose"
)
S <- S[match(names(lab_map), S$name), ] # reorder rows top -> bottom
# First factor level renders at the bottom, so reverse to put alfalfa at the top.
S$label <- factor(unname(lab_map[S$name]), levels = rev(unname(lab_map)))
mid_lab <- levels(S$label)[ceiling(nlevels(S$label) / 2)] # anchor for reference text

# Panel A values: matched comparison from the workflow's own fits. Observed share =
# point estimate on the traits used for the permutation null; null = mean of the
# permutation shares on those traits. Replaces the cross-crop summary columns, which
# averaged the observed share over all traits and took the null from one trait.
NS_PATH <- Sys.getenv("PANEL_NULL_SUMMARY",
  "experiments/empirical/multi/out/panel_workflow/panel_null_summary.csv")
NS <- read.csv(NS_PATH, check.names = FALSE, stringsAsFactors = FALSE)
ns_key <- c(
  alfalfa = "Alfalfa", blueberry_gc = "BB-UF-call*", blueberry2 = "BB-UF-fruit",
  blueberry = "BB-UF-GP", miscanthus = "Miscanthus", panicum = "Panicum",
  potato2 = "Potato-SLU", potato = "Potato-SolCAP", rose = "Rose"
)
for (i in seq_len(nrow(S))) {
  sub <- NS[NS$panel == ns_key[[S$name[i]]], ]
  tr <- sub[sub$trait_or_panel_summary != "PANEL_SUMMARY", ]
  pl <- sub[sub$trait_or_panel_summary == "PANEL_SUMMARY", ]
  stopifnot(nrow(pl) == 1, nrow(tr) >= 1)
  S$std_Vd[i] <- pl$s_hat_AD_real
  S$perm_null[i] <- pl$perm_mean
  S$std_Vd_trait_min[i] <- min(tr$s_hat_AD_real)
  S$std_Vd_trait_max[i] <- max(tr$s_hat_AD_real)
}

obs_col <- "#0072B2" # observed point / cosine / rank
range_col <- "grey65" # trait-range line
null_col <- "#d55e00" # permutation null — the ONLY red in the figure
ref_col <- "grey55" # dashed reference lines

# Grid on the value (x) axis only; calm horizontal dot-plot frame.
theme_dot <- theme_pub + theme_axes +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
  )
# B and C drop the repeated crop labels (rows align across panels).
theme_noy <- theme(axis.text.y = element_blank(), axis.title.y = element_blank())

## ---- Panel A: dominance share + trait range + permutation null ----------
key_lv <- c("Observed value", "Permutation null", "Trait range")
pA <- ggplot(S) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = ref_col, linewidth = 0.4) +
  annotate("text",
    x = 0.45, y = "Alfalfa", label = "hat(italic(s))[italic(D)] == 0.5", parse = TRUE,
    vjust = 0., size = 2.7, colour = range_col
  ) +
  geom_segment(aes(
    x = std_Vd_trait_min, xend = std_Vd_trait_max, y = label, yend = label,
    colour = "Trait range"
  ), linewidth = 1.1, lineend = "round") +
  geom_point(aes(x = std_Vd, y = label, colour = "Observed value"), shape = 16, size = 2.5) +
  geom_point(aes(x = perm_null, y = label, colour = "Permutation null"), shape = 18, size = 3.1) +
  scale_colour_manual(values = c(
    "Observed value" = obs_col, "Trait range" = range_col, "Permutation null" = null_col
  ), breaks = key_lv, name = NULL) +
  scale_x_continuous(limits = c(0.24, 0.70), breaks = seq(0.3, 0.7, 0.1)) +
  scale_y_discrete(labels = function(x) {
    exprs <- sapply(x, function(val) {
      if (val == "BB-UF-call") {
        "\"BB-UF-call\"^\"†\""
      } else {
        paste0("\"", val, "\"")
      }
    })
    parse(text = exprs)
  }) +
  guides(colour = guide_legend(override.aes = list(
    shape = c(16, 18, NA), linetype = c(0, 0, 1),
    linewidth = c(0, 0, 1.1), size = c(2.5, 3.1, 0)
  ))) +
  labs(
    title = "D-kernel share vs permutation null", tag = "A",
    x = expression("Standard " * italic(D) * "-kernel share, " * hat(italic(s))[italic(D)]), y = NULL
  ) +
  theme_dot


## ---- Panel B: additive-digenic Frobenius cosine -------------------------
pB <- ggplot(S, aes(AD_cos, label)) +
  geom_point(colour = obs_col, shape = 16, size = 2.5) +
  geom_vline(xintercept = 1, linetype = 2, colour = ref_col, linewidth = 0.4) +
  annotate("text",
    x = 1.01, y = mid_lab, label = "identical direction",
    angle = 90, hjust = 0.5, vjust = 1.4, size = 2.7, colour = range_col
  ) +
  scale_x_continuous(limits = c(0.85, 1.02), breaks = seq(0.85, 1.0, 0.05)) +
  labs(
    title = "A–D kernel collinearity", tag = "B",
    x = "A–D Frobenius cosine", y = NULL
  ) +
  theme_dot +
  theme_noy

## ---- Panel C: four-kernel effective rank --------------------------------
pC <- ggplot(S, aes(eff_rank, label)) +
  geom_point(colour = obs_col, shape = 16, size = 2.5) +
  geom_vline(xintercept = 4, linetype = 2, colour = ref_col, linewidth = 0.4) +
  annotate("text",
    x = 3.65, y = mid_lab, label = "nominal maximum = 4",
    angle = 90, hjust = 0.5, vjust = 1.4, size = 2.7, colour = range_col
  ) +
  scale_x_continuous(limits = c(1, 4.25), breaks = 1:4) +
  labs(
    title = "Four-kernel effective rank", tag = "C",
    x = "Effective rank", y = NULL
  ) +
  theme_dot +
  theme_noy

fig <- (pA | pB | pC) +
  plot_layout(widths = c(1.32, 1, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
  theme(legend.position = "bottom", legend.margin = margin(2, 4, 0, 4))
save_fig(fig, "fig_crosscrop", 250, 105, out = Sys.getenv("FIG_OUT", "paper/figures"))

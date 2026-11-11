#!/usr/bin/env Rscript
# =====================================================================
# Cross-crop PREDICTION figure (Fig. 4) — dominance-labelled variance is
# DECOUPLED from predictive gain.
#   Panel A  apparent standard dominance share  vs  realized A+D predictive
#            gain over additive (one point per panel, uniform marker size since
#            every panel now uses a comparable handful of traits;
#            alfalfa highlighted as the strongest positive case).
#   Panel B  per-panel forest plot of two paired gains: standard A+D minus A
#            (circle) and orthogonal A+D2+...+Dm minus standard A+D (triangle).
# A per-PANEL layout is used instead of a pooled trait-level boxplot because the
# panels carry different trait counts, so a pooled distribution would
# misrepresent the per-crop evidence. The current raw prediction outputs define
# the plotted trait set, keeping the x (Vd-share) and y (CV gain) consistent.
# Reads current per-repeat outputs in experiments/empirical/multi/out/panel_workflow/raw
# first, with historical
# per-dataset *_ortho.csv and *_PA.csv files as fallbacks.
# Output: paper/figures/fig_pred.png (1200 dpi) + .pdf
# Run from repo root:  Rscript experiments/empirical/multi/make_pred_figure.R
# =====================================================================
suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi)

OUT <- "experiments/empirical/multi/out"
RAW <- "experiments/empirical/multi/out/panel_workflow/raw"
rd <- function(p) if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL
finite <- function(x) x[is.finite(x)]
sd_or_na <- function(x) {
  x <- finite(x)
  if (length(x) > 1) sd(x) else NA_real_
}

# Canonical fitted standard D-kernel share (s_hat_AD_real) from the panel
# scorecard -- the SAME source behind Tables tab:pred / tab:crosscrop_cri /
# tab:workflow -- so panel A's x-axis matches the tables exactly. This uses the
# full-data diagnostic fit (alfalfa collapsed to its single representative
# yi_year trait), not the raw per-repeat, per-trait mean.
SCORE <- rd(file.path(OUT, "panel_workflow", "workflow_scorecard.csv"))
score_key <- c(
  alfalfa = "Alfalfa", blueberry_gc = "BB-UF-call*", blueberry2 = "BB-UF-fruit",
  blueberry = "BB-UF-GP", miscanthus = "Miscanthus", panicum = "Panicum",
  potato2 = "Potato-SLU", potato = "Potato-SolCAP", rose = "Rose",
  sugarcane = "Sugarcane", sweetpotato = "Sweetpotato"
)
canon_sd <- function(key) {
  if (is.null(SCORE) || is.na(score_key[key])) return(NA_real_)
  v <- SCORE$s_hat_AD_real[SCORE$panel == score_key[[key]]]
  if (length(v)) as.numeric(v[1]) else NA_real_
}

raw_pa <- function(key) {
  x <- rd(file.path(RAW, paste0(key, "__pred.csv")))
  if (is.null(x)) return(NULL)
  rows <- lapply(unique(x$trait), function(tr) {
    s <- x[x$trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      PA_A = mean(s$r_A, na.rm = TRUE),
      PA_A_sd = sd_or_na(s$r_A),
      PA_AD = mean(s$r_AD, na.rm = TRUE),
      PA_AD_sd = sd_or_na(s$r_AD),
      PA_orth = mean(s$r_orth, na.rm = TRUE),
      PA_orth_sd = sd_or_na(s$r_orth),
      d_AD_vs_A = mean(s$delta_AD, na.rm = TRUE),
      d_AD_vs_A_sd = sd_or_na(s$delta_AD),
      d_orth_vs_AD = mean(s$delta_orth, na.rm = TRUE),
      d_orth_vs_AD_sd = sd_or_na(s$delta_orth),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pa_for <- function(key) {
  pa <- raw_pa(key)
  if (!is.null(pa)) return(pa)
  rd(file.path(OUT, paste0(key, "_PA.csv")))
}

real_for <- function(key) {
  real <- rd(file.path(RAW, paste0(key, "__real.csv")))
  if (!is.null(real) && "s_hat_AD" %in% names(real)) {
    return(data.frame(trait = real$trait, std_Vd = real$s_hat_AD))
  }
  rd(file.path(OUT, paste0(key, "_ortho.csv")))
}

# Panel order top -> bottom, in the shortened Fig. 3 style. BB-UF-call carries
# the dagger (same individuals as BB-UF-GP, different dosage pipeline). The two
# higher-ploidy panels are listed last and separated by a divider.
cfg <- list(
  list(key = "alfalfa", lab = "Alfalfa", ref = TRUE, ploidy = "tetra"),
  list(key = "blueberry_gc", lab = "BB-UF-call†", ref = FALSE, ploidy = "tetra"),
  list(key = "blueberry2", lab = "BB-UF-fruit", ref = FALSE, ploidy = "tetra"),
  list(key = "blueberry", lab = "BB-UF-GP", ref = FALSE, ploidy = "tetra"),
  list(key = "miscanthus", lab = "Miscanthus", ref = FALSE, ploidy = "tetra"),
  list(key = "panicum", lab = "Panicum", ref = FALSE, ploidy = "tetra"),
  list(key = "potato2", lab = "Potato-SLU", ref = FALSE, ploidy = "tetra"),
  list(key = "potato", lab = "Potato-SolCAP", ref = FALSE, ploidy = "tetra"),
  list(key = "rose", lab = "Rose", ref = FALSE, ploidy = "tetra"),
  list(key = "sugarcane", lab = "Sugarcane (var.)", ref = FALSE, ploidy = "higher"),
  list(key = "sweetpotato", lab = "Sweetpotato (6×)", ref = FALSE, ploidy = "higher")
)

panel <- data.frame()
for (d in cfg) {
  pa <- pa_for(d$key)
  if (is.null(pa)) next
  gAD <- pa$d_AD_vs_A
  gOR <- pa$d_orth_vs_AD
  app_vd <- canon_sd(d$key)
  # uncertainty of the panel-mean gain: between-trait SE for multi-trait panels,
  # or the per-trait CV-rep SD for a single-trait panel.
  se2 <- function(x, sdcol) {
    x <- finite(x)
    if (length(x) > 1) sd(x) / sqrt(length(x)) else mean(sdcol, na.rm = TRUE)
  }
  panel <- rbind(panel, data.frame(
    key = d$key, lab = d$lab, ref = d$ref, ploidy = d$ploidy, app_vd = app_vd,
    gAD = mean(gAD, na.rm = TRUE), gAD_se = se2(gAD, pa$d_AD_vs_A_sd),
    gOR = mean(gOR, na.rm = TRUE), gOR_se = se2(gOR, pa$d_orth_vs_AD_sd),
    ntr = length(finite(gAD))
  ))
}
panel$who <- ifelse(panel$ref, "alfalfa", "other")
# rank 1 = top row; ypos puts rank 1 at the top of the y axis.
panel$rank <- seq_len(nrow(panel))
panel$ypos <- nrow(panel) + 1 - panel$rank

panel$lab_expr <- sapply(seq_len(nrow(panel)), function(i) {
  val <- panel$lab[i]
  who <- panel$who[i]
  if (grepl("†", val) || val == "BB-UF-call") {
    clean_val <- gsub("†", "", val)
    expr <- paste0("\"", clean_val, "\"^\"†\"")
  } else {
    expr <- paste0("\"", val, "\"")
  }
  if (who == "alfalfa") {
    paste0("bold(", expr, ")")
  } else {
    expr
  }
})

## ---- palette --------------------------------------------------------------
RED <- "#d6604d" # alfalfa highlight
BLU <- "#6aaed6" # other panels (Panel A fill)
STD <- "#0072B2" # standard A+D - A   (Panel B, Okabe-Ito blue)
ORT <- "#D55E00" # orthogonal - (A+D) (Panel B, Okabe-Ito vermillion)

## ---- Panel A: apparent dominance share vs realized A+D gain --------------
lab_col <- ifelse(panel$who == "alfalfa", RED, "grey25")
set.seed(7)
pA <- ggplot(panel, aes(app_vd, gAD)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(ymin = gAD - gAD_se, ymax = gAD + gAD_se),
    width = 0, linewidth = 0.4, colour = "grey55"
  ) +
  geom_point(aes(fill = who), size = 3.4, shape = 21, colour = "grey20", stroke = 0.35) +
  geom_text_repel(aes(label = lab_expr), parse = TRUE,
    colour = lab_col, size = 2.5,
    seed = 7, box.padding = 0.45, point.padding = 0.3, min.segment.length = 0.2,
    segment.colour = "grey70", segment.size = 0.3, max.overlaps = Inf
  ) +
  # annotate("text",
  #   x = 0.40, y = 0.083, hjust = 0, vjust = 1, size = 2.6, fontface = 3, colour = "grey40",
  #   label = "~0.4–0.6 dominance share in every panel,\nbut gain is near zero except alfalfa"
  # ) +
  scale_fill_manual(values = c(alfalfa = RED, other = BLU), guide = "none") +
  scale_x_continuous(limits = c(0.385, 0.63), breaks = seq(0.4, 0.6, 0.05)) +
  scale_y_continuous(limits = c(-0.024, 0.088)) +
  labs(
    title = "Standard D-kernel share vs predictive gain", tag = "A",
    x = expression("Fitted standard D-kernel share  " * italic(s)[D]),
    y = expression("Predictive gain of " * italic(A) + italic(D) * " over " * italic(A) * "  (" * Delta * italic(r) * ")")
  ) +
  theme_pub +
  theme_axes +
  theme(legend.position = "none")

## ---- Panel B: per-panel forest plot of the two paired gains --------------
off <- 0.17
fb <- rbind(
  data.frame(
    ypos = panel$ypos + off, val = panel$gAD, se = panel$gAD_se,
    series = "std", who = panel$who
  ),
  data.frame(
    ypos = panel$ypos - off, val = panel$gOR, se = panel$gOR_se,
    series = "orth", who = panel$who
  )
)
ser_lv <- c("std", "orth")
ser_lab <- c(
  std = expression(italic(A) + italic(D) ~ "−" ~ italic(A) ~ "  (Standard)"),
  orth = expression("Orthogonal" ~ "−" ~ (italic(A) + italic(D)))
)
fb$series <- factor(fb$series, levels = ser_lv)
ylabs_expr <- sapply(panel$lab[order(panel$ypos)], function(val) {
  if (grepl("†", val) || val == "BB-UF-call") {
    clean_val <- gsub("†", "", val)
    paste0("\"", clean_val, "\"^\"†\"")
  } else {
    paste0("\"", val, "\"")
  }
})
ylabs <- parse(text = ylabs_expr)
div_y <- max(panel$ypos[panel$ploidy == "higher"]) + 0.5 # divider above higher-ploidy block
alf_y <- panel$ypos[panel$who == "alfalfa"]

pB <- ggplot(fb, aes(val, ypos)) +
  annotate("rect",
    xmin = -Inf, xmax = Inf, ymin = alf_y - 0.5, ymax = alf_y + 0.5,
    fill = alpha(RED, 0.09)
  ) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = div_y, linetype = 3, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = val - se, xend = val + se, y = ypos, yend = ypos, colour = series),
    linewidth = 0.45
  ) +
  geom_point(aes(colour = series, shape = series), size = 2.1, stroke = 0.5, fill = "white") +
  annotate("text",
    x = 0.088, y = div_y - 0.34, hjust = 1, vjust = 1,
    size = 2.4, fontface = 3, colour = "grey55", label = "higher ploidy"
  ) +
  scale_colour_manual(values = c(std = STD, orth = ORT), labels = ser_lab, name = NULL) +
  scale_shape_manual(values = c(std = 16, orth = 17), labels = ser_lab, name = NULL) +
  scale_y_continuous(
    breaks = sort(panel$ypos), labels = ylabs,
    limits = c(0.4, max(panel$ypos) + 0.6), expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(-0.082, 0.092), breaks = seq(-0.08, 0.08, 0.02)) +
  labs(
    title = "Prediction gains are panel-specific", tag = "B",
    x = expression("Change in predictive ability  (" * Delta * italic(r) * ")"), y = NULL
  ) +
  theme_pub +
  theme_axes +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "inside", legend.position.inside = c(0.985, 0.03),
    legend.justification.inside = c(1, 0),
    legend.text = element_text(size = 7.6), legend.key.size = unit(11, "pt"),
    legend.background = element_rect(fill = alpha("white", 0.78), colour = "grey85", linewidth = 0.2),
    legend.margin = margin(2, 5, 3, 4),
    legend.spacing.y = unit(1, "pt")
  )

fig <- (pA | pB) +
  plot_layout(widths = c(1, 1.04)) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))
save_fig(fig, "fig_pred", 240, 110)

with(panel, cat(sprintf(
  "fig_pred: %d panels, %d traits | std A+D-A: alfalfa %+.3f, others median %+.3f | orth-(A+D) all <= 0 (median %+.3f)\n",
  nrow(panel), sum(ntr), gAD[who == "alfalfa"], median(gAD[who != "alfalfa"]), median(gOR)
)))

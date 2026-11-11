#!/usr/bin/env Rscript
# =====================================================================
# Structure-control figures (ggplot2 rewrite, shared theme).
#   fig_structure_vd     : standard V_D/(V_A+V_D) per panel under no control /
#                          +10 genomic PCs / +family-group kernel (grouped bars),
#                          dashed 0.5 reference.
#   fig_structure_cosine : Frobenius-cosine heatmap (A-D, A-PC, D-PC, A-F, D-F,
#                          D2-D3, D3-D4) across datasets; A-D column emphasised;
#                          n/a for diversity panels without family grouping.
# Reads experiments/empirical/multi/out/*_structure_{std,cosines}.csv
# Output: paper/figures/fig_structure_vd.png + fig_structure_cosine.png (1200 dpi) + .pdf
# Run from repo root:  Rscript experiments/empirical/multi/make_structure_figures.R
# =====================================================================
suppressMessages({
  library(ggplot2)
})
source("experiments/fig_theme.R") # theme_pub, theme_axes, save_fig (1200 dpi)

O <- "experiments/empirical/multi/out"
rd <- function(p) if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL
sets <- c("alfalfa", "potato", "blueberry", "blueberry2", "potato2", "panicum", "rose", "miscanthus", "blueberry_gc")
labs <- c(
  alfalfa = "Alfalfa", potato = "Potato SolCAP", blueberry = "Blueberry UF-GP", blueberry2 = "Blueberry UF-fruit",
  potato2 = "Potato SLU", panicum = "Panicum", rose = "Rose", miscanthus = "Miscanthus", blueberry_gc = "Blueberry UF-call*"
)

## ===================== (1) V_D before vs after structure control =========
rows <- list()
for (n in sets) {
  s <- rd(file.path(O, paste0(n, "_structure_std.csv")))
  if (is.null(s)) next
  # alfalfa's three yield slices are near-identical aggregations of one trait
  # (r >= 0.92); use the overall yi_year slice as the single representative.
  if (n == "alfalfa" && "trait" %in% names(s)) s <- s[s$trait == "yi_year", , drop = FALSE]
  ag <- tapply(s$std_Vd, s$control, mean, na.rm = TRUE)
  rows[[n]] <- data.frame(
    dataset = n, PC0 = ag["PC0"], PC10 = ag["PC10"],
    F = if ("F" %in% names(ag)) ag["F"] else NA_real_
  )
}
V <- do.call(rbind, rows)
# x-axis order matches Supplementary Table tab:datasets_supp (alphabetical: crop, then dataset)
panel_order <- c("alfalfa", "blueberry_gc", "blueberry2", "blueberry",
                 "miscanthus", "panicum", "potato2", "potato", "rose")
ord <- intersect(panel_order, V$dataset)
ctrl_lv <- c("PC0", "PC10", "F")
ctrl_lab <- c("No control (A+D)", "+ 10 genomic PCs", "+ family/group kernel")
Vl <- do.call(rbind, lapply(ctrl_lv, function(cc) data.frame(dataset = V$dataset, control = cc, vd = V[[cc]])))
Vl$dataset <- factor(labs[Vl$dataset], levels = labs[ord])
Vl$control <- factor(ctrl_lab[match(Vl$control, ctrl_lv)], levels = ctrl_lab)
ctrl_cols <- setNames(c("#6aaed6", "#3182bd", "#08519c"), ctrl_lab)

p_vd <- ggplot(Vl, aes(dataset, vd, fill = control)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74, na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.4) +
  scale_fill_manual(values = ctrl_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Standard dominance share before vs after structure control",
    x = NULL, y = expression(italic(V)[D] / (italic(V)[A] + italic(V)[D]))
  ) +
  theme_pub + theme_axes +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "inside", legend.position.inside = c(0.99, 0.99),
    legend.justification.inside = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = "grey85", linewidth = 0.2),
    legend.margin = margin(2, 4, 2, 4)
  )
save_fig(p_vd, "fig_structure_vd", 220, 105)

## ===================== (2) kernel-cosine heatmap =========================
ccols <- c("AD", "A_PC", "D_PC", "A_F", "D_F", "D2D3", "D3D4")
clabs <- c(AD = "A-D", A_PC = "A-PC", D_PC = "D-PC", A_F = "A-F", D_F = "D-F", D2D3 = "D2-D3", D3D4 = "D3-D4")
Crows <- lapply(sets, function(n) {
  c <- rd(file.path(O, paste0(n, "_structure_cosines.csv")))
  if (is.null(c)) NULL else data.frame(dataset = n, c[, ccols, drop = FALSE])
})
C <- do.call(rbind, Crows)
Cl <- do.call(rbind, lapply(ccols, function(cc) data.frame(dataset = C$dataset, pair = cc, cos = C[[cc]])))
# y-axis order matches Supplementary Table tab:datasets_supp; rev() puts the first (Alfalfa) at top
ord_c <- intersect(panel_order, C$dataset)
Cl$dataset <- factor(labs[Cl$dataset], levels = rev(labs[ord_c]))
Cl$pair <- factor(clabs[Cl$pair], levels = clabs[ccols])
ndat <- nlevels(Cl$dataset)
Cl$lab <- ifelse(is.na(Cl$cos), "n/a", sprintf("%.2f", Cl$cos))
Cl$lt <- !is.na(Cl$cos) & Cl$cos > 0.6 # light (white) text on dark cells

p_cos <- ggplot(Cl, aes(pair, dataset, fill = cos)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  annotate("rect",
    xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = ndat + 0.5,
    fill = NA, colour = "grey20", linewidth = 0.7
  ) +
  geom_text(data = Cl[Cl$pair != "A-D", ], aes(label = lab, colour = lt), size = 2.8) +
  geom_text(data = Cl[Cl$pair == "A-D", ], aes(label = lab, colour = lt), fontface = "bold", size = 2.8) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 1), na.value = "grey92", name = "Cosine") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  labs(title = "Kernel collinearity (Frobenius cosine)", x = NULL, y = NULL) +
  theme_pub +
  theme(
    panel.grid = element_blank(), axis.ticks = element_blank(),
    legend.position = "right", legend.key.width = unit(8, "pt"), legend.key.height = unit(1, "null")
  )
save_fig(p_cos, "fig_structure_cosine", 200, 118)

cat("STRUCTURE FIGURES DONE\n")

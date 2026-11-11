## =====================================================================
## Shared figure style for all ggplot figures in the paper.
## Single source of truth: canonical theme, shared encodings and a
## 1200-dpi PNG + vector-PDF writer. Source from the repo root:
##   source("experiments/fig_theme.R")
## Edit the theme HERE and every figure that sources it stays consistent.
## =====================================================================
suppressMessages({
  library(ggplot2)
  library(showtext)
  library(sysfonts)
})

## ---- pin Arial for EVERY figure, identical in PNG and PDF -----------------
## fontconfig on macOS mis-resolves "Arial" (-> Hiragino Sans), so register the
## real Arial .ttf files by path and let showtext draw the glyphs. showtext
## renders the same Arial outlines on both the raster (PNG) and vector (PDF)
## devices, so the two outputs are guaranteed identical and reproducible.
fig_font   <- "Arial"
.arial_dir <- "/System/Library/Fonts/Supplemental"
if (!(fig_font %in% sysfonts::font_families())) {
  sysfonts::font_add(fig_font,
    regular    = file.path(.arial_dir, "Arial.ttf"),
    bold       = file.path(.arial_dir, "Arial Bold.ttf"),
    italic     = file.path(.arial_dir, "Arial Italic.ttf"),
    bolditalic = file.path(.arial_dir, "Arial Bold Italic.ttf"))
}
showtext::showtext_auto()
## theme base_family covers axis/title/legend text; geom_text/geom_label/annotate
## need their default family set too, so every drawn character is Arial.
ggplot2::update_geom_defaults("text",  list(family = fig_font))
ggplot2::update_geom_defaults("label", list(family = fig_font))

## ---- population-structure encodings (simulation + calibration figures) ----
str_lv <- c("unrelated", "halfsib", "fullsib")
str_lab <- c("Unrelated", "Half-sib", "Full-sib")
str_shp <- c(unrelated = 16, halfsib = 17, fullsib = 15) # circle / triangle / square
str_cols <- c(unrelated = "black", halfsib = "#D55E00", fullsib = "#0072B2")
leg_name <- "Population structure"

## ---- canonical publication theme (matches fig_simulation) ----
## Text is grey20 throughout; panel labels (tags) are black. Grid/axis lines
## live in theme_axes so heatmap panels can omit them.
theme_pub <- theme_minimal(base_size = 11, base_family = fig_font) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, colour = "grey20", margin = margin(b = 5)),
    plot.tag = element_text(size = 13, face = "bold", colour = "black"),
    axis.title = element_text(size = 9, colour = "grey20"),
    axis.title.x = element_text(hjust = 0.5, margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.text = element_text(size = 9, colour = "grey20"),
    strip.text = element_text(size = 9, colour = "grey20"),
    legend.title = element_text(size = 9, colour = "grey20"),
    legend.text = element_text(size = 8.2, colour = "grey20"),
    plot.margin = margin(9, 10, 7, 9)
  )

## ---- grid + axis lines for scatter / line panels (add on top of theme_pub) ----
theme_axes <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
  axis.line = element_line(colour = "grey40", linewidth = 0.4),
  axis.ticks = element_line(colour = "grey40", linewidth = 0.3)
)

## ---- compact in-panel legend (bottom-right), matching fig_simulation panel A ----
theme_legend_inset <- theme(
  legend.position = "inside",
  legend.position.inside = c(0.985, 0.03),
  legend.justification.inside = c(1, 0),
  legend.title = element_text(size = 7.8, colour = "grey20"),
  legend.text = element_text(size = 7.3, colour = "grey20"),
  legend.key.size = unit(10, "pt"),
  legend.background = element_rect(fill = alpha("white", 0.75), colour = "grey85", linewidth = 0.2),
  legend.margin = margin(2, 4, 3, 4)
)

## ---- rotated x tick labels (matches fig_simulation) ----
theme_xrot <- theme(axis.text.x = element_text(angle = 28, hjust = 0.5, vjust = 0.65))

## ---- 1200-dpi PNG + vector PDF writer (journal requirement) ----
## w, h in millimetres. Writes <stem>.png (1200 dpi) and <stem>.pdf (vector).
save_fig <- function(plot, stem, w, h, out = "paper/figures") {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  ## showtext sizes text by the dpi it is told; match it to each device
  ## (1200 for the raster PNG, 72 for the vector PDF) so the Arial text is the
  ## same physical size in both.
  showtext::showtext_opts(dpi = 1200)
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = w, height = h, units = "mm", dpi = 1200, bg = "white", type = "cairo")
  ok <- tryCatch(
    {
      showtext::showtext_opts(dpi = 72)
      ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
      TRUE
    },
    error = function(e) tryCatch(
      {
        showtext::showtext_opts(dpi = 72)
        ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = w, height = h, units = "mm", bg = "white")
        TRUE
      },
      error = function(e2) FALSE
    )
  )
  showtext::showtext_opts(dpi = 1200)
  cat(sprintf("wrote %s.png (1200 dpi)%s\n", stem, if (ok) paste0(" + ", stem, ".pdf") else " (pdf failed)"))
}

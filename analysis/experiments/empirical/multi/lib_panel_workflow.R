# =====================================================================
# lib_panel_workflow.R  --  COMPATIBILITY SHIM
# -----------------------------------------------------------------------------
# The panel-conditioned diagnostic workflow now lives in the R package polyDx
# (software/polyDx). This file exists so the existing callers -- run_panel_workflow.R
# and the petrichor submission path -- keep working unchanged: it loads polyDx and
# sources lib_ortho.R, after which every symbol these scripts expect
# (build_panel, kernel_geometry, fit_fraction, share_pt, gen_*, load_panel, ...)
# resolves on the search path exactly as before.
#
# Resolution order for polyDx:
#   1. an installed polyDx  ->  library(polyDx)
#   2. POLYDX_SRC           ->  source() every R/*.R there (cluster / no-install path;
#                               submit_panel_workflow.sh stages them into $RUNDIR/polyDx_R)
#   3. software/polyDx/R relative to this file or the working directory
#
# POLYDX_SRC must point at the directory that directly holds the R/*.R sources
# (software/polyDx/R, NOT the package root software/polyDx -- the root has no
# top-level .R files and would otherwise resolve as if loading succeeded while
# sourcing nothing). A candidate directory with zero matching .R files is
# rejected, not selected; if no candidate yields any .R files this stop()s.
#
# Set POLYDX_FORCE_SRC to any value other than "0" to skip the installed-
# package check and always resolve from source (steps 2-3) -- useful when an
# installed polyDx would otherwise silently shadow freshly staged POLYDX_SRC
# sources on a machine where someone has since run R CMD INSTALL. Whichever
# path is taken, the resolved source is always cat()'d to the log: the
# literal directory when source-loaded, or "installed package polyDx
# <version>" when package-loaded.
#
# Env vars (PWF_NREF, PWF_FRAC_NITER, ...) are read by polyDx at call time, so the
# sbatch --export settings keep working untouched.
# =====================================================================

LIB_ORTHO <- Sys.getenv("LIB_ORTHO", "experiments/simulations/poly_ortho/lib_ortho.R")
source(LIB_ORTHO)
suppressMessages(library(stats))

.polydx_load <- function() {
  force_src <- Sys.getenv("POLYDX_FORCE_SRC", "0") != "0"
  if (!force_src && requireNamespace("polyDx", quietly = TRUE)) {
    suppressMessages(library(polyDx))
    cat(sprintf("polyDx source: installed package polyDx %s\n",
                as.character(utils::packageVersion("polyDx"))))
    return("package")
  }
  cand <- Sys.getenv("POLYDX_SRC", "")
  here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA_character_)
  cand <- c(cand,
            if (!is.na(here)) file.path(here, "polyDx_R"),
            if (!is.na(here)) file.path(here, "..", "..", "..", "software", "polyDx", "R"),
            "software/polyDx/R", "polyDx_R")
  cand <- cand[nzchar(cand) & dir.exists(cand)]
  err_msg <- paste0("polyDx not found: install it (R CMD INSTALL software/polyDx) or set POLYDX_SRC ",
                    "to a directory holding its R/*.R sources")
  if (!length(cand)) stop(err_msg, call. = FALSE)
  # A candidate directory with zero matching .R files (e.g. POLYDX_SRC pointing at the
  # package root instead of its R/ subdirectory) is rejected, not silently selected.
  src <- NULL
  for (c in cand) {
    if (length(list.files(c, pattern = "\\.R$"))) { src <- c; break }
  }
  if (is.null(src)) stop(err_msg, call. = FALSE)
  for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) {
    if (basename(f) == "polyDx-package.R") next
    source(f)
  }
  cat(sprintf("polyDx source: %s\n", normalizePath(src)))
  src
}
.POLYDX_SOURCE <- .polydx_load()

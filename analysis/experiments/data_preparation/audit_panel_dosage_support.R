#!/usr/bin/env Rscript
# =====================================================================
# audit_panel_dosage_support.R -- run the per-individual ploidy check that the
# Miscanthus rebuild introduced across EVERY panel input, so the same defect
# cannot be sitting unnoticed in another panel.
#
# Usage: Rscript experiments/data_preparation/audit_panel_dosage_support.R [DATA_DIR]
# Prints one line per panel; exits non-zero if any fixed-ploidy panel fails.
# =====================================================================
source("experiments/data_preparation/lib_dosage_checks.R")

args <- commandArgs(TRUE)
DIR  <- if (length(args) >= 1) args[1] else "data/processed/panel_inputs"

# panel key -> ploidy; "var" panels carry a per-marker ploidy vector and are
# reported but not asserted (their support is legitimately ragged).
PANELS <- list(alfalfa = 4L, potato = 4L, potato2 = 4L, blueberry = 4L, blueberry2 = 4L,
               blueberry_gc = 4L, panicum = 4L, rose = 4L, miscanthus = 4L,
               sweetpotato = 6L, sugarcane = NA_integer_)

fails <- 0L
for (k in names(PANELS)) {
  f <- file.path(DIR, paste0(k, "_dosage.rds"))
  if (!file.exists(f)) { cat(sprintf("%-14s SKIP   (no %s)\n", k, basename(f))); next }
  M  <- readRDS(f)
  pl <- PANELS[[k]]
  sup <- dosage_support(M)
  if (is.na(pl)) {
    cat(sprintf("%-14s VAR    n=%-5d support {%s} (variable per-marker ploidy; not asserted)\n",
                k, nrow(M), paste(sup, collapse = ",")))
    next
  }
  st  <- individual_support_stride(M)
  nb  <- sum(st > 1L)
  msg <- tryCatch({ assert_polyploid_support(M, pl, label = k); NULL },
                  error = function(e) conditionMessage(e))
  cat(sprintf("%-14s %-6s n=%-5d ploidy=%d support {%s} rank %d/%d  rescaled-ploidy individuals: %d\n",
              k, if (is.null(msg)) "PASS" else "FAIL", nrow(M), pl, paste(sup, collapse = ","),
              length(sup), pl + 1L, nb))
  if (!is.null(msg)) { cat("   ", msg, "\n"); fails <- fails + 1L }
}
cat(sprintf("\n%d panel(s) failed the per-individual ploidy check.\n", fails))
if (fails > 0L) quit(save = "no", status = 1L)

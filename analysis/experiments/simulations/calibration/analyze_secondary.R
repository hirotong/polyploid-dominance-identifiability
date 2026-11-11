#!/usr/bin/env Rscript
# Assemble the secondary calibration stress-grid per-task CSVs, dedupe, and compute
# calibration metrics per (n, M, LD) for the standard A+D share and the orthogonal
# D2 fraction. Answers plan section 9/15: does calibration improve with larger n,
# more markers, different LD; and does calibration error track A-D Frobenius cosine
# / effective rank?
#   in : experiments/simulations/calibration/out/calib_sec_t*.csv
#   out: experiments/simulations/calibration/out/calibration_secondary.csv  (deduped)
#        experiments/simulations/calibration/out/calibration_secondary_metrics.csv
# Run from repo root: Rscript experiments/simulations/calibration/analyze_secondary.R
suppressMessages({library(dplyr)})
IND <- Sys.getenv("CALIB_OUTDIR", "experiments/simulations/calibration/out")
fs <- list.files(IND, pattern = "^calib_sec_t.*\\.csv$", full.names = TRUE)
if (!length(fs)) {
  ready <- file.exists(file.path(IND, "calibration_secondary.csv")) &&
    file.exists(file.path(IND, "calibration_secondary_metrics.csv"))
  if (ready) {
    cat("No secondary calibration shard files found; consolidated outputs already exist in ", IND, "\n", sep = "")
    quit(status = 0)
  }
  stop("No secondary calibration shard files found in ", IND,
       ". If calibration_secondary.csv already exists, the shards have probably been consolidated.")
}
d <- do.call(rbind, lapply(fs, function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
cat(sprintf("read %d per-task files -> %d rows\n", length(fs), nrow(d)))

d <- d[d$failed == 0 & is.finite(d$standard_D_share), ]
# dedupe: one row per (N, M, LD, s, rep) -- the timeprobe rep1 duplicates an lgslowB rep1
key <- with(d, paste(N, nMarker, LD, target_D2, rep, sep = "|"))
d <- d[!duplicated(key), ]
d$LDlab <- ifelse(d$LD == "independent", "low-LD", ifelse(d$LD == "high", "high-LD", d$LD))
write.csv(d, file.path(IND, "calibration_secondary.csv"), row.names = FALSE)
cat(sprintf("deduped -> %d rows; cells (N x M x LD x s):\n", nrow(d)))
print(table(paste0("n", d$N), paste0("M", d$nMarker)))

## reps per cell (sanity)
cat("\nreps per (N,M,LD) cell (median across s):\n")
rp <- d |> group_by(N, nMarker, LD, target_D2) |> summarise(reps = n(), .groups = "drop") |>
  group_by(N, nMarker, LD) |> summarise(med_reps = median(reps), .groups = "drop")
print(as.data.frame(rp))

## ---- calibration metrics per (N, M, LD): est = alpha + beta * s ----
metric <- function(g, col) {
  s <- g$target_D2; est <- g[[col]]; ok <- is.finite(est)
  s <- s[ok]; est <- est[ok]; if (length(unique(s)) < 2) return(c(alpha = NA, beta = NA, rmse = NA, bias0 = NA))
  fit <- lm(est ~ s); b <- coef(fit)
  c(alpha = unname(b[1]), beta = unname(b[2]), rmse = sqrt(mean((est - s)^2)),
    bias0 = mean(est[abs(s) < 1e-6] - s[abs(s) < 1e-6]))
}
cells <- d |> group_by(N, nMarker, LD, LDlab) |> group_modify(~ {
  ms <- metric(.x, "standard_D_share"); mo <- metric(.x, "orth_D2_fraction")
  data.frame(adjLD_r2 = round(mean(.x$adjLD_r2, na.rm = TRUE), 3),
             AD_cos = round(mean(.x$A_D_Frobenius_cosine, na.rm = TRUE), 3),
             eff_rank = round(mean(.x$effective_rank, na.rm = TRUE), 3),
             std_alpha = round(ms["alpha"], 3), std_beta = round(ms["beta"], 3),
             std_rmse = round(ms["rmse"], 3), std_bias0 = round(ms["bias0"], 3),
             orth_alpha = round(mo["alpha"], 3), orth_beta = round(mo["beta"], 3),
             orth_rmse = round(mo["rmse"], 3), orth_bias0 = round(mo["bias0"], 3))
}) |> ungroup() |> as.data.frame()
write.csv(cells, file.path(IND, "calibration_secondary_metrics.csv"), row.names = FALSE)

if (Sys.getenv("CALIB_CLEAN_SHARDS", "0") == "1") {
  removed <- unlink(fs)
  cat(sprintf("removed %d per-task shard files\n", sum(removed == 0)))
}

cat("\n=== standard A+D share calibration by (n, M, LD) ===\n")
print(cells[order(cells$nMarker, cells$N, cells$LD),
            c("N","nMarker","LDlab","adjLD_r2","AD_cos","eff_rank","std_alpha","std_beta","std_rmse","std_bias0")], row.names = FALSE)
cat("\n=== orthogonal D2 calibration by (n, M, LD) (NA at n=19200: orth fit gated for cost) ===\n")
print(cells[order(cells$nMarker, cells$N, cells$LD),
            c("N","nMarker","LDlab","orth_alpha","orth_beta","orth_rmse","orth_bias0")], row.names = FALSE)

## ---- does calibration error track kernel geometry? ----
cat("\n=== calibration error vs kernel geometry (cell-level) ===\n")
cc <- cells[is.finite(cells$std_rmse), ]
cat(sprintf("cor(std_rmse, A-D cosine)  = %+.2f\n", cor(cc$std_rmse, cc$AD_cos, use = "complete.obs")))
cat(sprintf("cor(std_rmse, eff_rank)    = %+.2f\n", cor(cc$std_rmse, cc$eff_rank, use = "complete.obs")))
cat(sprintf("cor(std_bias0, A-D cosine) = %+.2f\n", cor(cc$std_bias0, cc$AD_cos, use = "complete.obs")))

## ---- trend summaries for the paper ----
cat("\n=== std intercept alpha by n (averaged over M,LD) ===\n")
print(aggregate(std_alpha ~ N, cells, function(z) round(mean(z), 3)))
cat("=== std intercept alpha by M (averaged over n,LD) ===\n")
print(aggregate(std_alpha ~ nMarker, cells, function(z) round(mean(z), 3)))
cat("=== std intercept alpha by LD (averaged over n,M) ===\n")
print(aggregate(std_alpha ~ LDlab, cells, function(z) round(mean(z), 3)))
cat("=== std bias0 (share at true s=0) by n x M ===\n")
print(reshape(aggregate(std_bias0 ~ N + nMarker, cells, function(z) round(mean(z), 3)),
              idvar = "N", timevar = "nMarker", direction = "wide"))
cat("\nSECONDARY ANALYSIS DONE\n")

#!/usr/bin/env Rscript
# Summarize the expanded regime sweep: prints the numbers for the Results
# paragraph and answers the 7 questions. Reads regime_map_expanded.csv.
CSV <- Sys.getenv("REGIME_CSV", "experiments/simulations/regime_map/out/regime_map_expanded.csv")
d <- read.csv(CSV, stringsAsFactors = FALSE)
d <- d[is.finite(d$spurious_Vd) & (is.na(d$failed) | d$failed == 0), ]
fs <- d[d$structure == "fullsib", ]
NS <- sort(unique(d$N)); MS <- sort(unique(d$nMarker)); LDS <- c("independent","low","medium","high")
cat(sprintf("rows=%d  cells=%d  fullsib rows=%d\n", nrow(d),
    nrow(unique(d[c("N","nMarker","LD","dr")])), nrow(fs)))
cat(sprintf("failed cells in file: %d\n\n", sum(d$failed == 1, na.rm = TRUE)))

p <- function(x) round(mean(x, na.rm = TRUE), 3)
## spurious share vs N (full-sib, averaged over M, LD, dr, reps)
cat("== spurious V_D/(V_A+V_D) vs N (full-sib, mean +/- SD) ==\n")
for (n in NS) { v <- fs$spurious_Vd[fs$N == n]
  cat(sprintf("  n=%6d : %.3f  +/- %.3f   (reps=%d)\n", n, mean(v), sd(v), length(v))) }
## structure agreement (where all 3 present)
cat("\n== spurious share vs N by structure (N<=4800) ==\n")
for (st in c("unrelated","halfsib","fullsib")) { s <- d[d$structure == st, ]
  cat(sprintf("  %-9s: ", st)); for (n in NS[NS<=4800]) cat(sprintf("n%d=%.2f ", n, p(s$spurious_Vd[s$N==n]))); cat("\n") }
## M effect (full-sib, averaged over N, LD, dr)
cat("\n== spurious share vs marker number M (full-sib) ==\n")
for (m in MS) cat(sprintf("  M=%7d : %.3f\n", m, p(fs$spurious_Vd[fs$nMarker == m])))
## M effect within fixed large N (does denser help at large n?)
cat("\n== share vs M at n=19200 (full-sib) ==\n")
for (m in MS) cat(sprintf("  M=%7d : %.3f\n", m, p(fs$spurious_Vd[fs$nMarker==m & fs$N==19200])))
## LD effect (full-sib): share + realized r2 by LD level, and vs N
cat("\n== LD effect (full-sib): realized r2 and spurious share by level ==\n")
for (l in LDS) { s <- fs[fs$LD == l, ]
  cat(sprintf("  %-11s realized r2=%.3f  share=%.3f\n", l, mean(s$adjLD_r2, na.rm=TRUE), p(s$spurious_Vd))) }
cat("  share vs N x LD:\n")
for (l in LDS) { s <- fs[fs$LD==l,]; cat(sprintf("   %-11s ", l)); for (n in NS) cat(sprintf("n%d=%.2f ", n, p(s$spurious_Vd[s$N==n]))); cat("\n") }
## dr effect
cat("\n== double-reduction effect (full-sib, share by dr) ==\n")
for (dv in sort(unique(fs$dr))) cat(sprintf("  dr=%.2f : %.3f\n", dv, p(fs$spurious_Vd[fs$dr==dv])))
cat("  share vs N x dr:\n")
for (dv in sort(unique(fs$dr))) { cat(sprintf("   dr=%.2f ", dv)); for (n in NS) cat(sprintf("n%d=%.2f ", n, p(fs$spurious_Vd[fs$dr==dv & fs$N==n]))); cat("\n") }
## diagnostics consistency
dg <- fs[is.finite(fs$AD_cos) & is.finite(fs$eff_rank), ]
cat(sprintf("\n== diagnostics (full-sib): AD_cos mean=%.3f [%.2f,%.2f], eff_rank mean=%.3f [%.2f,%.2f] ==\n",
    mean(dg$AD_cos), min(dg$AD_cos), max(dg$AD_cos), mean(dg$eff_rank), min(dg$eff_rank), max(dg$eff_rank)))
if (nrow(dg) > 5) cat(sprintf("  cor(spurious_Vd, AD_cos)=%.2f ; cor(spurious_Vd, eff_rank)=%.2f\n",
    cor(dg$spurious_Vd, dg$AD_cos), cor(dg$spurious_Vd, dg$eff_rank)))
## diagnostics vs N (does collinearity persist with n?)
cat("  AD_cos / eff_rank vs N:\n")
for (n in NS) { s <- dg[dg$N==n,]; if(nrow(s)) cat(sprintf("   n=%6d cos=%.3f eff_rank=%.3f\n", n, mean(s$AD_cos), mean(s$eff_rank))) }
cat("\nDONE\n")

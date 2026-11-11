#!/usr/bin/env Rscript
# Concatenate all per-task tier CSVs, DEDUPE on the deterministic key
# (N,nMarker,LD,dr,structure,rep) so reruns/gap-fills are idempotent, write the
# final regime_map_expanded.csv, and report any (N,M,LD,dr) cells below target
# reps (target 50 for N<=4800, 30 for N>=9600; full-sib is the binding structure).
# Usage: Rscript assemble_expanded.R <indir> <outcsv>
a <- commandArgs(TRUE)
indir <- if (length(a) >= 1) a[1] else "experiments/simulations/regime_map/out/expanded"
outf  <- if (length(a) >= 2) a[2] else "experiments/simulations/regime_map/out/regime_map_expanded.csv"
files <- list.files(indir, pattern = "^regime_exp_t.*\\.csv$", full.names = TRUE)  # originals (tg*) + top-up tiers; originals sort first so dedupe keeps reps 1-30, adds 31-50
cat(sprintf("reading %d task CSVs from %s\n", length(files), indir))
d <- do.call(rbind, lapply(files, function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)))
cat(sprintf("raw rows: %d\n", nrow(d)))
d <- d[!is.na(d$N) & d$structure != "ALL", ]                       # drop failed-cell placeholder rows
key <- paste(d$N, d$nMarker, d$LD, d$dr, d$structure, d$rep, sep = "|")
d <- d[!duplicated(key), ]                                         # idempotent dedupe
cat(sprintf("deduped rows: %d  (cells: %d)\n", nrow(d), nrow(unique(d[c("N","nMarker","LD","dr")]))))
write.csv(d, outf, row.names = FALSE); cat("wrote", outf, "\n")

## gap report: full-sib reps per (N,M,LD,dr) vs target
fs <- d[d$structure == "fullsib", ]
tgt <- function(N) 50   # uniform 50 reps/cell (top-up target; was 50 for N<=4800, 30 above)
ag <- aggregate(rep ~ N + nMarker + LD + dr, fs, function(x) length(unique(x)))
ag$target <- tgt(ag$N); ag$deficit <- ag$target - ag$rep
def <- ag[ag$deficit > 0, ]
cat(sprintf("\n== cells below target: %d (of %d) ==\n", nrow(def), nrow(ag)))
if (nrow(def)) {
  def <- def[order(-def$deficit), ]
  print(head(def[c("N","nMarker","LD","dr","rep","target","deficit")], 40))
  cat(sprintf("\ntotal missing full-sib reps: %d\n", sum(def$deficit)))
  ## emit the distinct (N,M) needing gap-fill (compact rerun targets)
  nm <- unique(def[c("N","nMarker")]); nm <- nm[order(nm$N, nm$nMarker), ]
  cat("gap-fill (N,M) targets:\n"); for (i in seq_len(nrow(nm))) cat(sprintf("  N=%d M=%d\n", nm$N[i], nm$nMarker[i]))
} else cat("no gaps — all cells at target.\n")

#!/usr/bin/env Rscript
# =====================================================================
# Petrichor SLURM array task: updog dosage calling on a contiguous marker chunk.
# Each array task processes markers chunk SLURM_ARRAY_TASK_ID of NCHUNK total.
# Input : $WD/blueberry_counts.rds = list(RO, AO), individuals x markers.
# Output: $WD/chunk_<task>.rds = list(idx, geno) with geno = markers x individuals.
# =====================================================================
suppressMessages(library(updog))
WD   <- Sys.getenv("WD", "/path/to/work/code/bb_dosage")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
K    <- as.integer(Sys.getenv("NCHUNK", "32"))
nc   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

cnt <- readRDS(file.path(WD, "blueberry_counts.rds"))   # list(RO, AO) indiv x markers
RO <- cnt$RO; AO <- cnt$AO; m <- ncol(RO)
chunks <- split(seq_len(m), ceiling(seq_len(m) / ceiling(m / K)))
idx <- chunks[[task]]
ref  <- t(RO[, idx, drop = FALSE])            # markers x individuals
size <- t((RO + AO)[, idx, drop = FALSE])     # markers x individuals
cat(sprintf("task %d/%d: %d markers, nc=%d\n", task, K, length(idx), nc))

md   <- multidog(refmat = ref, sizemat = size, ploidy = 4L, model = "norm", nc = nc)
geno <- format_multidog(md, varname = "geno") # markers x individuals (0..4)
saveRDS(list(idx = idx, geno = geno), file.path(WD, sprintf("chunk_%03d.rds", task)))
cat("CHUNK DONE\n")

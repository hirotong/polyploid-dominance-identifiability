library(polyDx)

## --- panel_registry -------------------------------------------------------
reg <- panel_registry()
stopifnot(length(reg) == 11L)
stopifnot(identical(names(reg)[1:3], c("alfalfa", "potato", "potato2")))
stopifnot(identical(reg$alfalfa$name, "Alfalfa"), identical(reg$alfalfa$ploidy, "4"))
stopifnot(identical(reg$sweetpotato$ploidy, "6"))
stopifnot(grepl("^var:", reg$sugarcane$ploidy))
stopifnot(identical(reg$rose$group, "none"))
stopifnot(all(vapply(reg, function(r) all(c("name", "ploidy", "group") %in% names(r)), logical(1))))

## --- a tiny in-memory panel ----------------------------------------------
set.seed(31)
n <- 50L; m <- 200L
M <- matrix(rbinom(n * m, 4, 0.4), n, m)
rownames(M) <- paste0("ind", seq_len(n))
phe <- data.frame(id = rownames(M), t1 = rnorm(n), t2 = rnorm(n), stringsAsFactors = FALSE)
grp <- rep(c("g1", "g2", "g3"), length.out = n)
P <- panel_from_matrix(M, phe, "4", group = grp, min_n = 40L)
P$name <- "TestPanel"; P$key <- "test"

## --- geom stage is phenotype-free and needs no BGLR ----------------------
KP <- build_panel(P$M, P$pcfg, need_X = FALSE)
g <- run_workflow_stage("geom", P, KP, panel_name = "TestPanel", panel_key = "test")
stopifnot(is.data.frame(g$table), nrow(g$table) == 1L)
stopifnot(all(c("panel", "panel_key", "ploidy", "n", "markers", "AD_cosine",
                "effective_rank", "ideal_rank", "rank_ratio", "condition_number",
                "max_abs_degree_cosine", "valid_markers") %in% names(g$table)))
stopifnot(identical(g$table$panel, "TestPanel"))
stopifnot("geom_cosine" %in% names(g$extra))
gc_tab <- g$extra$geom_cosine
stopifnot(nrow(gc_tab) == 16L)                    # 4 degrees x 4 degrees, long form
stopifnot(all(c("component_1", "component_2", "cosine") %in% names(gc_tab)))

## (M9) same seed -> identical output, for one cheap deterministic stage: geom
## needs no BGLR (see above) and every quantity it reports is a pure function
## of KP, so this exercises run_workflow_stage()'s reproducibility contract
## directly rather than relying on it being trivially true.
g1 <- run_workflow_stage("geom", P, KP, panel_name = "TestPanel", panel_key = "test", seed = 777L)
g2 <- run_workflow_stage("geom", P, KP, panel_name = "TestPanel", panel_key = "test", seed = 777L)
stopifnot(identical(g1$table, g2$table))
stopifnot(identical(g1$extra$geom_cosine, g2$extra$geom_cosine))

## an unknown stage is an error
stopifnot(inherits(try(run_workflow_stage("nope", P, KP), silent = TRUE), "try-error"))

if (!requireNamespace("BGLR", quietly = TRUE)) {
  cat("test-workflow.R OK (BGLR absent; fitting stages skipped)\n")
} else {
  Sys.setenv(PWF_FRAC_NITER = "300", PWF_FRAC_BURN = "100", PWF_FRAC_THIN = "5",
             PWF_PT_NITER = "300", PWF_PT_BURN = "100")
  KPX <- build_panel(P$M, P$pcfg, need_X = TRUE)

  ## real
  r <- run_workflow_stage("real", P, KP, panel_name = "TestPanel", panel_key = "test")
  stopifnot(is.data.frame(r$table), nrow(r$table) == length(P$traits))
  stopifnot(all(c("trait", "s_hat_AD", "P_D2", "L_D2", "D2_CrI_width",
                  "S_PC_real", "S_family_real", "S_struct_real") %in% names(r$table)))
  ## S_struct_real = max(0, S_PC_real, S_family_real) (workflow.R); recompute
  ## from the row's own columns rather than just >= 0, which max(0, ...) makes
  ## true by construction and so can never fail
  expected_S_struct_real <- with(r$table, mapply(function(pc, fam) {
    max(0, pc, if (is.finite(fam)) fam else -Inf)
  }, S_PC_real, S_family_real))
  stopifnot(max(abs(r$table$S_struct_real - expected_S_struct_real)) < 1e-12)

  ## perm
  p <- run_workflow_stage("perm", P, KP, panel_name = "TestPanel", panel_key = "test",
                          nperm = 2L, perm_traits = 1L)
  stopifnot(nrow(p$table) == 2L)
  stopifnot(all(c("trait", "perm_id", "s_hat_AD_perm") %in% names(p$table)))

  ## addctrl (needs X)
  a <- run_workflow_stage("addctrl", P, KPX, panel_name = "TestPanel", panel_key = "test", nref = 2L)
  stopifnot(nrow(a$table) == 2L, all(a$table$true_s_D2 == 0))
  stopifnot("D2_largest" %in% names(a$table))

  ## d2spike (needs X)
  d <- run_workflow_stage("d2spike", P, KPX, panel_name = "TestPanel", panel_key = "test",
                          nref = 1L, s_spike = c(0.5))
  stopifnot(nrow(d$table) == 1L, abs(d$table$true_s_D2 - 0.5) < 1e-12)

  ## structspike (needs X; PC + family because the panel has a group)
  s <- run_workflow_stage("structspike", P, KPX, panel_name = "TestPanel", panel_key = "test",
                          nref = 1L, lambda = c(0, 1))
  stopifnot(nrow(s$table) == 4L)                  # 2 types x 2 lambdas x 1 rep
  stopifnot(setequal(unique(s$table$structure_type), c("PC", "family")))
  ## S_struct = max(0, S_PC, S_family) (workflow.R); recompute from the row's
  ## own columns rather than just >= 0, which max(0, ...) makes true by
  ## construction and so can never fail
  expected_S_struct <- with(s$table, mapply(function(pc, fam) {
    max(0, pc, if (is.finite(fam)) fam else -Inf)
  }, S_PC, S_family))
  stopifnot(max(abs(s$table$S_struct - expected_S_struct)) < 1e-12)

  ## pred
  pr <- run_workflow_stage("pred", P, KP, panel_name = "TestPanel", panel_key = "test", cv_reps = 1L)
  stopifnot(nrow(pr$table) == length(P$traits))
  stopifnot(all(c("r_A", "r_AD", "r_orth", "delta_AD", "delta_orth") %in% names(pr$table)))
  stopifnot(all(abs(pr$table$delta_AD - (pr$table$r_AD - pr$table$r_A)) < 1e-12))

  ## (M9) pred's paired-folds property: the A, AD and orth models must see
  ## IDENTICAL fold partitions at each rep -- that pairing is the entire
  ## meaning of delta_AD/delta_orth (a model-vs-model comparison on the same
  ## held-out individuals); if it broke, every prediction-gain number would
  ## silently become noise. run_unit()'s fold assignment and BGLR fit are both
  ## internal to run_workflow_stage(), so this verifies the property from the
  ## outside by temporarily stubbing BGLR::BGLR to record, for every
  ## (rep, fold) call it receives, which individuals were held out (NA'd) in
  ## y -- independent of which model (A/AD/orth) made the call -- and to skip
  ## the real (slow) MCMC fit entirely. The stub is restored in `finally` (a
  ## top-level on.exit() would be a no-op here, per the note below).
  held <- list()
  bglr_ns <- asNamespace("BGLR")
  orig_bglr <- get("BGLR", envir = bglr_ns)
  fake_bglr <- function(y, ETA, nIter, burnIn, verbose, saveAt, ...) {
    ## key = "<rep>_<fold>_", i.e. saveAt's basename with the model prefix
    ## stripped, so A/AD/orth calls for the same (rep, fold) collide on key
    key <- sub("^[A-Za-z]+_", "", basename(saveAt))
    held[[key]][[length(held[[key]]) + 1L]] <<- which(is.na(y))
    list(yHat = rep(0, length(y)))
  }
  tryCatch({
    assignInNamespace("BGLR", fake_bglr, ns = bglr_ns)
    invisible(run_workflow_stage("pred", P, KP, panel_name = "TestPanel", panel_key = "test",
                                 cv_reps = 2L))
  }, finally = {
    assignInNamespace("BGLR", orig_bglr, ns = bglr_ns)
  })
  stopifnot(length(held) > 0L)
  for (k in names(held)) {
    grp <- held[[k]]
    stopifnot(length(grp) >= 2L)                          # >=2 models called for this (rep, fold)
    stopifnot(all(vapply(grp[-1], function(x) identical(x, grp[[1]]), logical(1))))
  }

  ## a simulator stage without contrast matrices is a clear error
  stopifnot(inherits(try(run_workflow_stage("addctrl", P, KP, nref = 1L), silent = TRUE), "try-error"))

  ## a top-level on.exit() in a script is a no-op (nothing left to trigger it),
  ## so the env vars set above are cleared explicitly here instead.
  Sys.unsetenv(c("PWF_FRAC_NITER", "PWF_FRAC_BURN", "PWF_FRAC_THIN",
                 "PWF_PT_NITER", "PWF_PT_BURN"))
  cat("test-workflow.R OK\n")
}

## Internal: robust parallel map for run_workflow_stage(). Ported verbatim
## (tryCatch -> NULL-on-error behaviour) from run_panel_workflow.R's pmap(),
## with mclapply's worker count derived from `cores` and the item count.
pwf_pmap <- function(items, fn, cores = 1L) {
  parallel::mclapply(items, function(it) tryCatch(fn(it),
                     error = function(e) { cat("  ERR:", conditionMessage(e), "\n"); NULL }),
                     mc.cores = max(1L, min(cores, length(items))))
}

#' Panel registry for the manuscript's cross-crop diagnostic workflow
#'
#' The eleven panels the manuscript's panel-conditioned diagnostic workflow
#' runs over, each with its display name, ploidy specification (a
#' \code{\link{parse_ploidy}} argument string) and group/family RDS filename
#' (relative to a data directory) or \code{"none"}. Ported verbatim from
#' \code{run_panel_workflow.R}'s \code{REG} object.
#'
#' @return a named list, one element per panel key, each element a list with
#'   \code{name} (display name), \code{ploidy} (a \code{\link{parse_ploidy}}
#'   specification string; the sugarcane entry keeps an unresolved \code{\%s}
#'   placeholder for the data directory, filled in with \code{sprintf()} by
#'   the caller) and \code{group} (a group/family RDS filename, or
#'   \code{"none"}).
#' @examples
#' reg <- panel_registry()
#' names(reg)
#' reg$alfalfa
#' @export
panel_registry <- function() {
  list(
    alfalfa      = list(name = "Alfalfa",      ploidy = "4", group = "alfalfa_group.rds"),
    potato       = list(name = "Potato-SolCAP",ploidy = "4", group = "potato_group.rds"),
    potato2      = list(name = "Potato-SLU",   ploidy = "4", group = "potato2_group.rds"),
    blueberry    = list(name = "BB-UF-GP",     ploidy = "4", group = "none"),
    blueberry2   = list(name = "BB-UF-fruit",  ploidy = "4", group = "none"),
    blueberry_gc = list(name = "BB-UF-call*",  ploidy = "4", group = "none"),
    panicum      = list(name = "Panicum",      ploidy = "4", group = "panicum_group.rds"),
    rose         = list(name = "Rose",         ploidy = "4", group = "none"),
    miscanthus   = list(name = "Miscanthus",   ploidy = "4", group = "miscanthus_group.rds"),
    sweetpotato  = list(name = "Sweetpotato",  ploidy = "6", group = "none"),
    sugarcane    = list(name = "Sugarcane",    ploidy = "var:%s/sugarcane_marker_ploidy.rds:6", group = "none")
  )
}

#' Run one stage of the panel-conditioned diagnostic workflow
#'
#' Runs a single stage of the manuscript's panel-conditioned diagnostic
#' workflow battery on an already-loaded panel and its already-built
#' kernels: the phenotype-free kernel-geometry diagnostics (\code{geom}),
#' real-trait variance fits (\code{real}), a phenotype-permutation null
#' (\code{perm}), and three panel-conditioned known-truth controls fit on
#' the panel's own genotype matrix — additive-only (\code{addctrl}), digenic
#' (D2) spike-in (\code{d2spike}) and structure spike-in (\code{structspike})
#' — plus a paired cross-validation comparison of additive, additive+
#' dominance and orthogonal-degree prediction (\code{pred}). Ported from the
#' corresponding \code{if (stage == ...)} block of
#' \code{run_panel_workflow.R}, with identical arithmetic, seeds and column
#' layout; only the CLI-driver plumbing (argument parsing, CSV writing,
#' \code{quit()}) is removed in favour of returning the result.
#'
#' @param stage one of \code{"geom"}, \code{"real"}, \code{"perm"},
#'   \code{"addctrl"}, \code{"d2spike"}, \code{"structspike"} or
#'   \code{"pred"}; any other value is an error.
#' @param panel a \code{\link{load_panel}}/\code{\link{panel_from_matrix}}
#'   result (needs \code{M}, \code{phe}, \code{traits}, \code{group},
#'   \code{pcfg}, \code{n}, \code{markers}, \code{qc}; \code{name}/\code{key}
#'   are used only as defaults for \code{panel_name}/\code{panel_key}).
#' @param kernels a \code{\link{build_panel}} result. Stages \code{addctrl},
#'   \code{d2spike} and \code{structspike} additionally need it built with
#'   \code{need_X = TRUE} (or \code{kernels$X} added via
#'   \code{\link{build_contrasts12}}).
#' @param panel_name display name recorded in the output table's
#'   \code{panel} column; defaults to \code{panel$name}, falling back to
#'   \code{"panel"}.
#' @param panel_key key recorded in the output table's \code{panel_key}
#'   column; defaults to \code{panel$key}, falling back to \code{"panel"}.
#' @param nperm number of permutations per trait for the \code{perm} stage.
#' @param nref number of simulated replicates per cell for the
#'   \code{addctrl}/\code{d2spike}/\code{structspike} stages.
#' @param seed base random seed; every stage's per-unit seed is derived
#'   deterministically from this plus the unit's index, and \code{seed}
#'   itself is set once at entry (matching the driver's global
#'   \code{set.seed(SEED)} before stage dispatch). Exception: the \code{pred}
#'   stage's fold assignment is seeded independently of \code{seed}, as
#'   \code{set.seed(1000 + r)} for cv rep \code{r} -- deliberately, so the
#'   same fold partition is reused for the \code{A}, \code{AD} and
#'   \code{orth} models within a rep regardless of \code{seed}, keeping
#'   \code{delta_AD}/\code{delta_orth} a paired comparison across models.
#' @param perm_traits number of (leading) traits the \code{perm} stage runs
#'   the permutation null on.
#' @param perm_ortho if \code{TRUE}, the \code{perm} stage also fits the
#'   orthogonal-degree model on each permutation and records
#'   \code{P_D2_perm}/\code{L_D2_perm}.
#' @param cv_reps number of repeats for the \code{pred} stage's paired CV.
#' @param h2 broad-sense heritability for the simulated known-truth
#'   controls; defaults to \code{\link{PWF_H2_DEFAULT}}.
#' @param s_spike numeric vector of target true digenic (D2) fractions for
#'   the \code{d2spike} stage.
#' @param lambda numeric vector of target structure fractions for the
#'   \code{structspike} stage.
#' @param cores number of parallel workers (via \code{\link[parallel]{mclapply}}).
#' @return a list with \code{table} (a \code{data.frame}, or \code{NULL}
#'   when a stage has nothing to report) and \code{extra} (a named list of
#'   auxiliary \code{data.frame}s; empty except for the \code{geom} stage's
#'   \code{geom_cosine}).
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(40 * 60, 4, 0.4), 40, 60)
#' rownames(M) <- paste0("ind", seq_len(40))
#' P <- panel_from_matrix(M, phe = NULL, "4")
#' KP <- build_panel(P$M, P$pcfg, need_X = FALSE)
#' g <- run_workflow_stage("geom", P, KP, panel_name = "Demo", panel_key = "demo")
#' g$table$AD_cosine
#' names(g$extra)
#' @export
run_workflow_stage <- function(stage, panel, kernels,
                               panel_name = panel$name %||% "panel",
                               panel_key  = panel$key  %||% "panel",
                               nperm      = pwf_opt("PWF_NPERM", 30L),
                               nref       = pwf_opt("PWF_NREF", 30L),
                               seed       = pwf_opt("PWF_SEED", 2026L),
                               perm_traits = pwf_opt("PWF_PERM_TRAITS", 3L),
                               perm_ortho  = pwf_opt("PWF_PERM_ORTHO", 0L) == 1L,
                               cv_reps    = pwf_opt("PWF_REPS", 5L),
                               h2         = PWF_H2_DEFAULT,
                               s_spike    = c(0.25, 0.50, 0.75),
                               lambda     = c(0, 0.25, 0.50, 0.75, 1.00),
                               cores      = 1L) {
  set.seed(seed)

  need_X <- function() {
    if (is.null(kernels$X)) {
      stop(sprintf("run_workflow_stage(): stage '%s' needs kernels built with need_X = TRUE", stage), call. = FALSE)
    }
  }

  # =====================================================================
  # STAGE: geom  (kernel-geometry diagnostics row + full degree cosine matrix)
  # =====================================================================
  if (stage == "geom") {
    geo <- kernel_geometry(kernels)
    df <- data.frame(panel = panel_name, panel_key = panel_key, ploidy = panel$pcfg$label,
                     n = panel$n, raw_markers = panel$qc$raw_markers,
                     markers_pre_qc = panel$qc$markers_pre_qc,
                     markers = panel$markers,
                     markers_dropped_qc = panel$qc$markers_dropped_qc,
                     valid_markers = kernels$valid,
                     AD_cosine = geo$AD_cosine, max_abs_degree_cosine = geo$max_abs_degree_cosine,
                     effective_rank = geo$eff_rank, ideal_rank = geo$ideal_rank,
                     rank_ratio = geo$rank_ratio, condition_number = geo$condition_number,
                     stringsAsFactors = FALSE)
    cm <- as.data.frame(as.table(geo$cosine_matrix)); names(cm) <- c("component_1", "component_2", "cosine")
    cm <- cbind(panel = panel_name, panel_key = panel_key, cm)
    return(list(table = df, extra = list(geom_cosine = cm)))
  }

  # =====================================================================
  # STAGE: real  (real-trait variance fits)
  # =====================================================================
  if (stage == "real") {
    if (!length(panel$traits)) return(list(table = data.frame(), extra = list()))
    geo <- kernel_geometry(kernels)
    ids <- rownames(panel$M); PC <- panel_pc(kernels, 10L); Fk <- family_kernel(panel$group, ids)
    rows <- pwf_pmap(panel$traits, function(tr) {
      y <- suppressWarnings(as.numeric(panel$phe[[tr]])); ok <- is.finite(y)
      if (sum(ok) < 40 || stats::sd(y[ok]) == 0) return(NULL)
      fr <- fit_fraction(y, kernels)            # full y (NAs handled by BGLR) + full kernels (manuscript convention)
      # real-trait structure sensitivity (same 3000/750 point machinery as the structure ruler)
      raw  <- share_pt(y, kernels, PC = PC, Fk = Fk, what = "raw")
      pc10 <- share_pt(y, kernels, PC = PC, Fk = Fk, what = "PC10")
      fam  <- if (!is.null(Fk)) share_pt(y, kernels, PC = PC, Fk = Fk, what = "family") else NA_real_
      S_PC <- raw - pc10; S_fam <- if (is.finite(fam)) raw - fam else NA_real_
      cbind(data.frame(panel = panel_name, panel_key = panel_key, trait = tr, ploidy = panel$pcfg$label,
                       n = sum(ok), markers = panel$markers,
                       AD_cosine = geo$AD_cosine, effective_rank = geo$eff_rank, stringsAsFactors = FALSE), fr,
            data.frame(s_hat_AD_raw_pt = raw, s_hat_AD_PC10 = pc10, s_hat_AD_family = fam,
                       S_PC_real = S_PC, S_family_real = S_fam,
                       S_struct_real = max(0, S_PC, if (is.finite(S_fam)) S_fam else -Inf), stringsAsFactors = FALSE))
    }, cores = cores)
    res <- do.call(rbind, Filter(is.data.frame, rows))
    return(list(table = res, extra = list()))
  }

  # =====================================================================
  # STAGE: perm  (phenotype-permutation null on representative trait[s])
  # =====================================================================
  if (stage == "perm") {
    if (!length(panel$traits)) return(list(table = data.frame(), extra = list()))
    use_tr <- panel$traits[seq_len(min(perm_traits, length(panel$traits)))]
    grid <- expand.grid(trait = use_tr, perm = seq_len(nperm), stringsAsFactors = FALSE)
    rows <- pwf_pmap(seq_len(nrow(grid)), function(i) {
      tr <- grid$trait[i]; pid <- grid$perm[i]
      y <- suppressWarnings(as.numeric(panel$phe[[tr]]))
      if (sum(is.finite(y)) < 40) return(NULL)
      set.seed(seed + 101 * pid + which(use_tr == tr))
      yp <- sample(y)                          # permute full vector (manuscript convention)
      sperm <- share_pt(yp, kernels, what = "raw")
      out <- data.frame(panel = panel_name, panel_key = panel_key, trait = tr, perm_id = pid,
                        s_hat_AD_perm = sperm, stringsAsFactors = FALSE)
      if (perm_ortho) { fr <- fit_fraction(yp, kernels); out$P_D2_perm <- fr$P_D2; out$L_D2_perm <- fr$L_D2 }
      out
    }, cores = cores)
    res <- do.call(rbind, Filter(is.data.frame, rows))
    return(list(table = res, extra = list()))
  }

  # =====================================================================
  # STAGE: addctrl  (additive-only control; true D2 = 0)
  # =====================================================================
  if (stage == "addctrl") {
    need_X()
    geo <- kernel_geometry(kernels)
    rows <- pwf_pmap(seq_len(nref), function(rep) {
      set.seed(seed + 211 * rep)
      y <- gen_additive_only(kernels, H2 = h2); fr <- fit_fraction(y, kernels)
      cbind(data.frame(panel = panel_name, panel_key = panel_key, rep = rep, H2 = h2, true_s_D2 = 0,
                       AD_cosine = geo$AD_cosine, effective_rank = geo$eff_rank, stringsAsFactors = FALSE), fr)
    }, cores = cores)
    res <- do.call(rbind, Filter(is.data.frame, rows))
    return(list(table = res, extra = list()))
  }

  # =====================================================================
  # STAGE: d2spike  (D2 spike-in positive control)
  # =====================================================================
  if (stage == "d2spike") {
    need_X()
    geo <- kernel_geometry(kernels)
    grid <- expand.grid(s = s_spike, rep = seq_len(nref))
    rows <- pwf_pmap(seq_len(nrow(grid)), function(i) {
      s <- grid$s[i]; rep <- grid$rep[i]
      set.seed(seed + 311 * rep + round(1000 * s))
      y <- gen_d2_spikein(kernels, s, H2 = h2); fr <- fit_fraction(y, kernels)
      cbind(data.frame(panel = panel_name, panel_key = panel_key, rep = rep, true_s_D2 = s, H2 = h2,
                       AD_cosine = geo$AD_cosine, effective_rank = geo$eff_rank, stringsAsFactors = FALSE), fr)
    }, cores = cores)
    res <- do.call(rbind, Filter(is.data.frame, rows))
    return(list(table = res, extra = list()))
  }

  # =====================================================================
  # STAGE: structspike  (structure-driven false-signal control)
  # =====================================================================
  if (stage == "structspike") {
    need_X()
    geo <- kernel_geometry(kernels)
    ids <- rownames(panel$M); PC <- panel_pc(kernels, 10L); Fk <- family_kernel(panel$group, ids)
    stypes <- c("PC", if (!is.null(panel$group)) "family")
    grid <- expand.grid(structure_type = stypes, lambda = lambda, rep = seq_len(nref), stringsAsFactors = FALSE)
    rows <- pwf_pmap(seq_len(nrow(grid)), function(i) {
      st <- grid$structure_type[i]; lam <- grid$lambda[i]; rep <- grid$rep[i]
      set.seed(seed + 411 * rep + 37 * round(100 * lam) + ifelse(st == "PC", 1, 2))
      y <- gen_structure(kernels, PC, panel$group, lam, st, H2 = h2)
      raw  <- share_pt(y, kernels, PC = PC, Fk = Fk, what = "raw")
      pc10 <- share_pt(y, kernels, PC = PC, Fk = Fk, what = "PC10")
      fam  <- if (!is.null(Fk)) share_pt(y, kernels, PC = PC, Fk = Fk, what = "family") else NA_real_
      S_PC <- raw - pc10; S_fam <- if (is.finite(fam)) raw - fam else NA_real_
      data.frame(panel = panel_name, panel_key = panel_key, rep = rep, structure_type = st, lambda = lam, H2 = h2,
                 s_hat_AD_raw = raw, s_hat_AD_PC10 = pc10, s_hat_AD_family = fam,
                 S_PC = S_PC, S_family = S_fam, S_struct = max(0, S_PC, if (is.finite(S_fam)) S_fam else -Inf),
                 AD_cosine = geo$AD_cosine, effective_rank = geo$eff_rank, stringsAsFactors = FALSE)
    }, cores = cores)
    res <- do.call(rbind, Filter(is.data.frame, rows))
    return(list(table = res, extra = list()))
  }

  # =====================================================================
  # STAGE: pred  (paired 5-fold x reps prediction CV; A, A+D, orthogonal)
  # Replicates run_pred_cv.R's run_unit (same seeds -> identical, paired folds).
  # Records per-rep pooled-fold predictive ability r.
  # =====================================================================
  if (stage == "pred") {
    if (!length(panel$traits)) return(list(table = data.frame(), extra = list()))
    if (!requireNamespace("BGLR", quietly = TRUE)) {
      stop("package 'BGLR' is required for the 'pred' stage; install.packages(\"BGLR\")", call. = FALSE)
    }
    mk_eta <- function(model) {
      if (model == "A")  return(list(A = list(K = kernels$A, model = "RKHS")))
      if (model == "AD") return(list(A = list(K = kernels$A, model = "RKHS"), D = list(K = kernels$D, model = "RKHS")))
      e <- list(); for (k in seq_len(kernels$ndeg)) e[[kernels$etanames[k]]] <- list(K = kernels$O[[k]], model = "RKHS"); e
    }
    run_unit <- function(tr, model, r, folds = 5, nIter = 1500, burnIn = 400) {
      y <- suppressWarnings(as.numeric(panel$phe[[tr]])); ok <- which(is.finite(y))
      set.seed(1000 + r); fold <- sample(rep(1:folds, length.out = length(ok)))
      ETA <- mk_eta(model); d <- tempfile("cvt"); dir.create(d); pr <- rep(NA_real_, length(y))
      for (f in 1:folds) { te <- ok[fold == f]; yN <- y; yN[te] <- NA
        fm <- BGLR::BGLR(y = yN, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE,
                   saveAt = file.path(d, sprintf("%s_%d_%d_", model, r, f))); pr[te] <- fm$yHat[te] }
      unlink(d, recursive = TRUE)
      data.frame(trait = tr, model = model, rep = r, pa = suppressWarnings(stats::cor(pr[ok], y[ok])))
    }
    units <- expand.grid(trait = panel$traits, model = c("A", "AD", "orth"), r = seq_len(cv_reps), stringsAsFactors = FALSE)
    parts <- pwf_pmap(seq_len(nrow(units)), function(i) run_unit(units$trait[i], units$model[i], units$r[i]), cores = cores)
    raw <- do.call(rbind, Filter(is.data.frame, parts))
    rows <- list()
    for (tr in unique(raw$trait)) for (r in sort(unique(raw$rep))) {
      g <- function(md) { v <- raw$pa[raw$trait == tr & raw$model == md & raw$rep == r]; if (length(v)) v[1] else NA_real_ }
      rA <- g("A"); rAD <- g("AD"); rO <- g("orth")
      rows[[length(rows) + 1]] <- data.frame(panel = panel_name, panel_key = panel_key, trait = tr,
        `repeat` = r, fold = "pooled", r_A = rA, r_AD = rAD, r_orth = rO,
        delta_AD = rAD - rA, delta_orth = rO - rAD, check.names = FALSE, stringsAsFactors = FALSE)
    }
    res <- do.call(rbind, rows)
    return(list(table = res, extra = list()))
  }

  stop(sprintf("unknown stage '%s'", stage), call. = FALSE)
}

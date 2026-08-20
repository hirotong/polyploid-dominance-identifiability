# polyDx — diagnostics for dominance-labelled variance in polyploid genomic prediction

Reusable software for the paper _"Biological identifiability of dominance-labelled variance
components in polyploid genomic prediction."_ `polyDx` is an installable R package that builds
**orthogonal additive / digenic / trigenic / quadrigenic relationship kernels** from allele-dosage
data (generalizing NOIA to arbitrary even ploidy, with explicit higher-order intralocus terms; see
Appendix A of the manuscript), and ships the **panel-conditioned diagnostic
workflow** the manuscript uses to decide, per empirical panel, whether a fitted dominance share is
biologically real or a collinearity artefact.

## Why

Standard autotetraploid dominance kernels (digenic; AGHmatrix Slater/Endelman) report a large
"dominance variance" share `s_hat_AD` of roughly 0.4–0.6 in essentially every real panel — including
when the true dominance is zero — because the additive and digenic kernels built from a real marker
matrix are nearly collinear (A–D cosine typically 0.91–0.98). Under that geometry an A+D model
splits genetic variance ~50/50 between A and D regardless of the true architecture: `s_hat_AD ≈ 0.5`
is **the null**, not evidence of dominance. `polyDx` supplies two things to get past that:

1. the **orthogonal degree decomposition** (additive, digenic, trigenic, quadrigenic, ...), which is
   identifiable in expectation even where the standard A/D split is not; and
2. the **panel-conditioned decision procedure** — known-truth controls (additive-only, digenic
   spike-in, structure spike-in, permutation null) regenerated from each empirical panel's own
   genotype matrix, and a rule-based scorecard that reads the real-trait signal against those
   controls before licensing a dominance claim.

## Install

```r
# install.packages("remotes")
remotes::install_github("hirotong/polyploid-dominance-identifiability")

install.packages("BGLR")   # optional; needed for variance-component fitting
library(polyDx)
```

Or from a clone: `R CMD INSTALL .` at the repository root. The package depends only on base R
(`stats`, `parallel`); `BGLR` is suggested and needed only for the fitting and workflow functions
(`run_rkhs()`, `fit_fraction()`, `share_pt()`, `partition_variance()`, `run_workflow_stage()`).

## Minimal worked example (kernels)

```r
library(polyDx)
set.seed(1)
n <- 200; m <- 400
M <- matrix(rbinom(n * m, 4, rep(runif(m, 0.2, 0.8), each = n)), n, m)
rownames(M) <- paste0("ind", seq_len(n))
y <- as.numeric(scale(rowSums(M[, 1:50]))) + rnorm(n)

KO <- orthogonal_kernels(M)   # list(A, D2, D3, D4) -- orthogonal degree kernels
KS <- standard_kernels(M)     # list(A, D)          -- VanRaden + digenic comparator

fr_orth <- partition_variance(y, KO)   # identifiable fraction per degree (needs BGLR)
fr_std  <- partition_variance(y, KS)   # the non-identifiable standard A+D share
```

## Worked example (the diagnostic workflow)

This is the example that backs the abstract's claim that the package implements the paper's
diagnostic workflow. It uses small panel and MCMC-chain settings purely so the example finishes in
under a minute on a laptop; the printed verdict is illustrative, not a real diagnostic result.

```r
library(polyDx)
set.seed(2)

# demo-scale MCMC chains; the manuscript uses 12000/3000/5 (fraction fits) and 3000/750 (point
# estimates) -- see PWF_FRAC_NITER/PWF_FRAC_BURN/PWF_FRAC_THIN/PWF_PT_NITER/PWF_PT_BURN in fit.R
Sys.setenv(PWF_FRAC_NITER = "600", PWF_FRAC_BURN = "150", PWF_FRAC_THIN = "5",
           PWF_PT_NITER = "600", PWF_PT_BURN = "150")

n <- 80L; m <- 300L
M <- matrix(rbinom(n * m, 4, 0.4), n, m); rownames(M) <- paste0("ind", seq_len(n))
phe <- data.frame(id = rownames(M), trait1 = rnorm(n))

P  <- panel_from_matrix(M, phe, ploidy_arg = "4", min_n = 20L)  # QC + trait selection
KP <- build_panel(P$M, P$pcfg, need_X = TRUE)                   # A, D, degree kernels + contrasts

kernel_geometry(KP)[c("AD_cosine", "eff_rank", "condition_number")]

stages <- c("geom", "real", "perm", "addctrl", "d2spike", "structspike")
tabs <- lapply(stages, function(s)
  run_workflow_stage(s, P, KP, panel_name = "MyPanel", panel_key = "mypanel",
                     nperm = 5L, nref = 5L)$table)   # manuscript runs use 30; kept tiny here for speed
names(tabs) <- stages

workflow_scorecard(tabs)[, c("panel", "delta_null", "D2_recovery_rate_s050",
                             "S_struct_real", "suggested_label")]
```

`workflow_scorecard()` assigns each panel exactly one of four `suggested_label` values:

- **`structure-sensitive`** — the real-trait dominance signal is at or below what population
  structure / relatedness alone produces on this panel; the fit is confounded, not dominance.
- **`D2-supported`** — the real signal exceeds the panel's own permutation null, the panel's
  geometry can recover a planted digenic (D2) spike, the real localisation resembles that spike
  rather than an additive-only null, the recovered credible interval is narrow, and the panel is
  not structure-sensitive. **This is the only verdict that licenses a biological dominance claim.**
- **`prediction-useful`** — not D2-supported, but adding the dominance kernel gives a material,
  reliably positive cross-validated prediction-accuracy gain over additive alone (useful in
  practice, without implying the gain is _biological_ dominance).
- **`model-basis`** — none of the above; the fitted dominance share is a modelling artefact of the
  panel's own kernel geometry, not evidence of anything biological or predictively useful.

## API

### Kernels & contrasts

| function                                 | purpose                                                                                                 |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `ortho_contrasts(p, ploidy)`             | orthonormal dosage contrasts (degree 0..ploidy) under binomial genotype frequencies                     |
| `impute_mean_dosage(M)`                  | mean-impute missing dosage calls by marker                                                              |
| `validate_dosage_matrix(M, ploidy, ...)` | validate and filter a dosage matrix before kernel construction                                          |
| `ortho_etanames(ndeg)`                   | names of the orthogonal degree-wise components, e.g. `c("A","D2","D3","D4")`                            |
| `orthogonal_kernels(M, ploidy)`          | orthogonal degree-wise relationship kernels `A, D2, D3, D4, ...` from dosage                            |
| `standard_kernels(M, ploidy)`            | standard additive (VanRaden) and digenic dominance kernels from dosage                                  |
| `build_panel(M, pcfg, need_X)`           | standard A/D plus orthogonal degree kernels and contrasts, from a panel's dosage matrix and ploidy spec |
| `build_contrasts12(...)`                 | degree-1 and degree-2 orthonormal contrast matrices only                                                |

### Panel & QC

| function                                             | purpose                                                                      |
| ---------------------------------------------------- | ---------------------------------------------------------------------------- |
| `parse_ploidy(arg, nmarkers)`                        | parse a fixed or per-marker variable ploidy specification string             |
| `subset_ploidy(pcfg, keep)`                          | subset a ploidy specification to a retained marker set                       |
| `filter_markers(M, pcfg, ...)`                       | marker QC (call rate, MAF, polymorphism) then mean imputation                |
| `load_panel(dosage_rds, pheno_rds, ploidy_arg, ...)` | load a panel from RDS files: align dosage/phenotype/group, QC, select traits |
| `panel_from_matrix(M, phe, ploidy_arg, ...)`         | build a panel from in-memory dosage, phenotype and group objects             |

### Geometry

| function                      | purpose                                                                                       |
| ----------------------------- | --------------------------------------------------------------------------------------------- |
| `kernel_geometry(KP)`         | phenotype-free kernel-geometry diagnostics: `AD_cosine`, `eff_rank`, `condition_number`, etc. |
| `check_kernels(KP, ids, tol)` | validate panel kernels before fitting (ID order, missingness, symmetry, PSD)                  |

### Fitting

| function                              | purpose                                                                   |
| ------------------------------------- | ------------------------------------------------------------------------- |
| `run_rkhs(y, ETA, etanames, ...)`     | fit a BGLR RKHS model and return posterior variance-fraction intervals    |
| `fit_fraction(y, KP, ...)`            | standard A+D and orthogonal-degree dominance-share fit metrics            |
| `share_pt(y, KP, PC, Fk, what, ...)`  | fast point-estimate dominance share, optionally controlled for structure  |
| `partition_variance(y, kernels, ...)` | partition fitted variance across an arbitrary set of kernels (needs BGLR) |
| `panel_pc(KP, npc)`                   | top genomic PCs from a panel's additive kernel                            |
| `family_kernel(group, ids)`           | family/group incidence kernel                                             |
| `q_named(x, p)`                       | unnamed quantile helper                                                   |

### Simulators (known-truth controls)

| function                                                   | purpose                                                                                 |
| ---------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `sim_comp(Xd, valid)`                                      | standardized degree-_d_ genetic value drawn in a panel's own orthonormal contrast space |
| `gen_additive_only(KP, H2)`                                | additive-only known-truth phenotype (true digenic fraction 0)                           |
| `gen_d2_spikein(KP, s, H2)`                                | digenic (D2) spike-in known-truth phenotype at a target true D2 fraction                |
| `gen_structure(KP, PC, group, lambda, structure_type, H2)` | population-structure known-truth phenotype                                              |
| `simulate_pure_degree(M, causal, degree)`                  | construct a pure degree-_k_ genotypic value (falsification-check simulator)             |

### Workflow & scorecard

| function                                         | purpose                                                                                                                             |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `panel_registry()`                               | the manuscript's 11-panel registry: display name, ploidy spec, group/family file                                                    |
| `run_workflow_stage(stage, panel, kernels, ...)` | run one stage (`geom`/`real`/`perm`/`addctrl`/`d2spike`/`structspike`/`pred`) of the diagnostic workflow on an already-loaded panel |
| `workflow_scorecard(tables, thresholds, ...)`    | rule-based scorecard turning the stage tables into one verdict row per panel                                                        |
| `workflow_thresholds()`                          | the gate cut-offs used by `workflow_scorecard()`                                                                                    |
| `pwf_opt(name, default)`                         | read a polyDx run-time option from an environment variable, coerced to `default`'s type                                             |
| `PWF_H2_DEFAULT`                                 | default broad-sense heritability (0.6) for the panel-conditioned simulators                                                         |

## The decision rule

`workflow_scorecard()` combines five boolean gates, each compared against a threshold from
`workflow_thresholds()`:

| gate                  | condition                                                                                                      | threshold             |
| --------------------- | -------------------------------------------------------------------------------------------------------------- | --------------------- |
| `structure_sensitive` | real structure-injection signal at or above the panel's own structure ruler                                    | `struct_large = 0.10` |
| `exceeds_null`        | real dominance share exceeds the panel's permutation-null mean by more than a margin                           | `null_margin = 0.05`  |
| `recoverable` (D2)    | fraction of digenic spike-in replicates where D2 is correctly identified as the largest non-additive component | `recover_rate = 0.80` |
| `narrow_ci`           | the real-trait D2 credible interval is narrower than a width cut-off                                           | `ci_wide = 0.25`      |
| `pred_gain`           | mean paired cross-validated prediction gain of A+D over A alone, with the lower CI bound above zero            | `pred_gain = 0.01`    |

`resembles_spike` (real localisation closer to the D2-spike distribution than to the additive-only
null, and above the 95th percentile of the additive-only null) is required alongside `exceeds_null`,
`recoverable` and `narrow_ci`. The conjunction

```
d2_supported <- exceeds_null && recoverable && resembles_spike && narrow_ci && !structure_sensitive
```

is the only condition under which `suggested_label` is `"D2-supported"`; everything else falls
through to `"structure-sensitive"` (checked first), `"prediction-useful"` (if `pred_gain` holds),
or `"model-basis"`.

`workflow_scorecard()` implements only this **generic** decision rule from stage tables the caller
supplies. The manuscript's published `workflow_scorecard.csv` and figures were produced by an
aggregation script in the authors' research repository (`make_workflow_scorecard.R`) that calls the
same gate logic but additionally single-sources some real-trait columns (`P_D2`, `L_D2`,
`D2_CrI_width`, `s_hat_AD`, kernel-geometry and permutation-null figures) for the nine tetraploid
panels from the manuscript's curated cross-crop table, so that the scorecard and the manuscript's
cross-crop table report identical numbers; it also writes the workflow summary and figures. None of
that overriding or figure generation happens inside the package. The generic rule does reproduce
the manuscript's published `suggested_label` on all 11 panels (Alfalfa, Potato-SolCAP, Potato-SLU
and Sugarcane land on `"prediction-useful"`; the other seven on `"model-basis"`), though individual
intermediate numeric columns for the tetraploid panels can differ from the published table where
that single-sourcing was applied.

## Relation to the manuscript

This repository ships the reusable methods: the orthogonal contrasts and kernels, the
kernel-geometry diagnostics, the panel-conditioned known-truth simulators, the per-stage workflow
runner and the rule-based scorecard. The manuscript's full analysis pipeline -- the simulation
studies, the per-panel SLURM drivers, and the aggregation script that produces the published
scorecard, summary and figures -- lives in the authors' research repository and is described in the
paper. All 11 empirical panels analysed in the manuscript are public; sources are listed in the
paper's Data and software availability statement.

## Limitations

Orthogonality holds in expectation under linkage equilibrium; finite-sample LD and strong structure
reintroduce some collinearity. Double reduction is handled by re-orthonormalizing under
double-reduction-adjusted genotype frequencies (Appendix A.7). Identifiability is fundamentally
limited at small n; report uncertainty.

The public helper validates dosage range, missingness and polymorphism before building kernels. Use
`impute_mean_dosage()` or call kernel builders with `impute = TRUE` for matrices containing `NA`; use
`min_maf` and `min_call_rate` to match manuscript-style QC filters.

## License

MIT (c) 2026 Jinguang Tong. See [`LICENSE.md`](LICENSE.md).

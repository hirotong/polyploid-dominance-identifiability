# polyDx — diagnostics for dominance-labelled variance in polyploid genomic prediction

Software and analysis code for the manuscript _"Limited biological identifiability of
dominance-labelled variance components in polyploid genomic prediction."_ `polyDx` is an installable
R package that builds **orthogonal additive / digenic / trigenic / quadrigenic relationship
kernels** from allele-dosage data (orthogonal dosage contrasts for any even ploidy, with explicit
higher-order intralocus terms; see the Materials and Methods of the manuscript), and ships the
**panel-conditioned diagnostic workflow** the manuscript uses to decide, per empirical panel,
whether a fitted dominance share can be read as biological digenic dominance.

The repository has two parts:

- the **R package** at the repository root (`R/`, `man/`, `tests/`);
- [`analysis/`](analysis/README.md): the scripts and summary tables behind every figure and table
  of the manuscript.

## Why

A standard additive-plus-dominance model reports a dominance share `s_hat_AD` that is not a
calibrated estimate of digenic dominance. In the manuscript's autotetraploid simulations the share
was about 0.4 when true digenic dominance was zero, and it rose by only 0.28 per unit of true
digenic fraction. In the eleven public panels the additive and dominance kernels built from the
marker matrix were nearly collinear (A–D cosine 0.90 to 0.99). A large fitted share is therefore
not, by itself, evidence of dominance. `polyDx` supplies two things for reading such a share:

1. the **orthogonal degree decomposition** (additive, digenic, trigenic, quadrigenic, ...). The
   degrees are orthogonal in expectation under linkage equilibrium and binomial genotype
   frequencies; in finite panels this improves calibration but does not restore it; and
2. the **panel-conditioned decision procedure**: known-truth controls (additive-only, digenic
   spike-in, structure spike-in, permutation reference) generated from each empirical panel's own
   genotype matrix, and a rule-based scorecard that compares the real-trait result with those
   controls.

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

This example runs every stage of the workflow on a small simulated panel. It uses small panel and
MCMC-chain settings so that it finishes in under a minute on a laptop; the printed label is
illustrative, not a real diagnostic result.

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

- **`structure-sensitive`**: the fitted dominance share changes materially when population
  structure or family membership is controlled, so it cannot be separated from structure on this
  panel.
- **`D2-supported`**: the real share exceeds the panel's own permutation reference, the panel can
  recover a planted digenic (D2) spike, the real localization resembles that spike rather than an
  additive-only control, the credible interval of the D2 fraction is narrow, and the panel is not
  structure-sensitive. This is the only label under which the workflow reads the share as digenic
  dominance.
- **`prediction-useful`**: not D2-supported, but adding the dominance kernel gives a material
  cross-validated gain in prediction accuracy over the additive model, with the lower interval
  bound above zero. The label says nothing about whether the gain reflects biological dominance.
- **`model-basis`**: none of the above. This label is not evidence that digenic dominance is
  absent. It records that the criteria for support were not met, which can also result from weak
  signal, limited power or model misspecification.

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

`workflow_scorecard()` implements this **generic** decision rule from stage tables the caller
supplies. The manuscript's scorecard (Supplemental Table S6,
`analysis/experiments/empirical/multi/out/panel_workflow/workflow_scorecard.csv`) was produced by
the aggregation script `analysis/experiments/empirical/multi/make_workflow_scorecard.R`. That
script applies the same gates and differs from the package function in two ways. For the nine
autotetraploid panels it takes some real-trait columns (`P_D2`, `L_D2`, `D2_CrI_width`, `s_hat_AD`
and the kernel geometry) from the cross-crop fits, so that the scorecard and the manuscript's
cross-crop table report the same numbers. It also compares the observed share with the permutation
reference on the same traits and with the same estimator. The generic rule reproduces the
manuscript's `suggested_label` on all eleven panels and on the simulated positive control (Alfalfa,
Potato-SolCAP, Potato-SLU and Sugarcane are `"prediction-useful"`; the other seven panels are
`"model-basis"`; the positive control is `"D2-supported"`). Intermediate numeric columns can differ
from the published table for the reasons above.

## Relation to the manuscript

The package holds the reusable methods: the orthogonal contrasts and kernels, the kernel-geometry
diagnostics, the panel-conditioned known-truth simulators, the per-stage workflow runner and the
rule-based scorecard. The manuscript's analysis is in [`analysis/`](analysis/README.md): the
simulation studies, the per-panel drivers and job scripts, the aggregation scripts, and the summary
tables from which every figure and table was made. All eleven empirical panels are public; their
sources are listed in the manuscript's Data Availability Statement. The genotype and phenotype
files themselves are not redistributed here.

## Limitations

Orthogonality holds in expectation under linkage equilibrium and binomial genotype frequencies;
finite-sample LD and population structure reintroduce collinearity. Double reduction is not
modelled: it moves genotype frequencies away from the binomial reference, so the contrasts are
then not orthogonal even at the population level. The decision thresholds are conventions, not
calibrated error rates, and the workflow's labels rest on one estimator with default priors
(BGLR). Identifiability is limited at small n; report uncertainty.

The public helper validates dosage range, missingness and polymorphism before building kernels. Use
`impute_mean_dosage()` or call kernel builders with `impute = TRUE` for matrices containing `NA`; use
`min_maf` and `min_call_rate` to match manuscript-style QC filters.

## License

MIT (c) 2026 Jinguang Tong. See [`LICENSE.md`](LICENSE.md).

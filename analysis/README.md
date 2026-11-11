# Analysis scripts and summary tables

This folder holds the scripts and the summary tables behind the figures and tables of the
manuscript _"Limited biological identifiability of dominance-labelled variance components in
polyploid genomic prediction."_ The reusable methods are in the `polyDx` R package at the
repository root; the scripts here call it.

## How the folder is organised

Paths inside the scripts are relative to this folder, so run every command from `analysis/`:

```bash
cd analysis
Rscript experiments/simulations/poly_ortho/make_fig_simulation.R   # writes paper/figures/
```

| Folder | Content |
| --- | --- |
| `experiments/simulations/` | known-truth simulations: degree attribution, calibration grids, stress grid |
| `experiments/empirical/multi/` | per-panel fits, permutation reference, structure control, prediction, the panel-conditioned workflow, estimator sensitivity |
| `experiments/empirical/multi/out/` | result tables of those analyses (one file per panel and analysis) |
| `experiments/aggregation/` | aggregation of the cross-validated prediction results |
| `experiments/data_preparation/` | construction of panel inputs from the public sources |
| `experiments/slurm/`, `experiments/empirical/multi/petrichor*/` | job scripts that record the settings of the cluster runs |
| `tables/` | cross-crop summary table |
| `paper/figure_source_data/` | one workbook per figure with the plotted values |

Panel keys used in file names: `alfalfa`, `potato` (Potato-SolCAP), `potato2` (Potato-SLU),
`blueberry` (BB-UF-GP), `blueberry2` (BB-UF-fruit), `blueberry_gc` (BB-UF-call), `panicum`,
`rose`, `miscanthus`, `sweetpotato`, `sugarcane`, and `poscontrol` (the simulated positive
control).

## Figures and tables

All paths below start at `experiments/`.

| Manuscript item | Scripts | Summary tables |
| --- | --- | --- |
| Figure 1 (degree attribution in simulation) | `simulations/poly_ortho/poly_ortho_framework.R`, `make_fig_simulation.R` | `simulations/poly_ortho/out/std_vd_share.csv`, `ortho_attribution.csv` |
| Supplemental Table S3 (alternative dominance kernels) | `simulations/poly_ortho/poly_ortho_framework.R` | `simulations/poly_ortho/out/alt_kernel_share.csv` |
| Figure 2, Supplemental Table S2 (calibration) | `simulations/calibration/calibration_sim.R`, `make_fig_calibration.R` | `simulations/calibration/out/calibration_primary.csv`, `calibration_metrics.csv` |
| Figure 3 (calibration across sample size, marker number and LD) | `simulations/calibration/calibration_secondary.R`, `analyze_secondary.R`, `make_fig_calibration_4panel.R` | `simulations/calibration/out/calibration_secondary.csv`, `calibration_secondary_metrics.csv` |
| Double-reduction sensitivity grid | `simulations/regime_map/regime_map_expanded.R`, `assemble_expanded.R`, `analyze_expanded.R` | `simulations/regime_map/out/regime_map_expanded.csv` |
| Table 1, Figure 4, Supplemental Table S4 (cross-crop diagnostics) | `empirical/multi/run_ortho_dataset.R`, `run_ortho_perm.R`, `make_figures.R`, `make_fig_crosscrop.R` | `empirical/multi/out/<panel>_{ortho,null,ortho_perm,diagnostics}.csv`, `../tables/cross_crop_summary.csv`, `empirical/multi/out/panel_workflow/panel_null_summary.csv` |
| Table 2, Figure 5 (prediction) | `empirical/multi/run_pred_cv.R`, `aggregation/prediction_cv/aggregate_predcv.py`, `empirical/multi/make_pred_figure.R` | `empirical/multi/out/<panel>_PA.csv` (alfalfa: `simulations/poly_ortho/out/real_ortho_roza_PA_paired.csv`), `empirical/multi/out/panel_workflow/raw/<panel>__pred.csv`, `aggregation/prediction_cv/*.csv` |
| Supplemental Table S1 (panels after quality control) | `empirical/multi/lib_marker_qc.R`, `data_preparation/` | `empirical/multi/out/panel_workflow/raw/<panel>__geom.csv`, `<panel>__qcaudit.csv` |
| Supplemental Figure S1 (structure control) | `empirical/multi/run_structure_dataset.R`, `make_structure_figures.R` | `empirical/multi/out/<panel>_structure_*.csv` |
| Supplemental Table S5 (estimator and prior sensitivity) | `empirical/multi/estimator_sensitivity.R`, `lib_reml.R`, `aggregate_estimator_sensitivity.py` | `empirical/multi/out/estimator_sensitivity/` |
| Figure 6 labels, Supplemental Table S6 (workflow scorecard) | `empirical/multi/run_panel_workflow.R`, `build_poscontrol_panel.R`, `make_workflow_scorecard.R` | `empirical/multi/out/panel_workflow/` (`raw/` holds one file per panel and stage) |
| Supplemental Table S7 and threshold sensitivity | `empirical/multi/threshold_sensitivity.py` | `empirical/multi/out/panel_workflow/threshold_sensitivity.csv` |
| Source data of the figures | `make_figure_source_data.py` | `../paper/figure_source_data/*.xlsx` |

Figure 6 itself is a hand-drawn schematic and has no plotted data.

## What can be rerun from this repository

- **Figures and aggregated tables** can be regenerated from the summary tables included here
  (the `make_*` and `aggregate_*` scripts and `threshold_sensitivity.py`).
- **Simulations** can be rerun from scratch. They need `AlphaSimR`, `AGHmatrix` and `BGLR`, and
  were run as cluster jobs; the job scripts give the replicate numbers and seeds.
- **Empirical fits** need the panel inputs. All eleven panels are public and are cited in the
  manuscript's Data Availability Statement, but the genotype and phenotype files are not
  redistributed here. The scripts in `experiments/data_preparation/` and the `build_*` scripts
  show how each panel input was built from its source.

## Software

R (4.4), with `polyDx` (this repository), `BGLR`, `AlphaSimR`, `AGHmatrix`, `data.table`, `updog`
(blueberry dosage calling), and `ggplot2`, `dplyr`, `patchwork`, `ggrepel`, `showtext` for the
figures; Python 3 with `numpy`, `pandas` and `openpyxl` for the aggregation scripts.

Install the package before running the workflow scripts:

```r
remotes::install_github("hirotong/polyploid-dominance-identifiability")
```

## Notes

- In the job scripts, cluster-specific paths and the scheduler account are replaced with
  placeholders (`/path/to/work`, `/path/to/scratch`, `YOUR_ACCOUNT`).
- The simulated population in `poly_ortho_framework.R` is not fixed by the seed across runs, so a
  rerun reproduces the reported means within replicate noise, not the individual replicate values.

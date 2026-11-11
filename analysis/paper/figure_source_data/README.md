# Figure source data

One Excel workbook per paper figure, holding the raw data used to plot it.
Regenerate with: `python3 experiments/make_figure_source_data.py` (from repo root).

Each workbook contains the underlying rows plus the aggregated quantities that
become the points / tiles / bars in the published figure.

| Figure | Workbook | Source data | Plotting script |
|---|---|---|---|
| Fig 1 — known-truth simulation | `Figure1_simulation.xlsx` | `experiments/simulations/poly_ortho/out/{std_vd_share,ortho_attribution}.csv` | `experiments/simulations/poly_ortho/make_fig_simulation.R` |
| Fig 2 — calibration curves | `Figure2_calibration.xlsx` | `experiments/simulations/calibration/out/calibration_primary.csv` (+ `calibration_metrics.csv`) | `experiments/simulations/calibration/make_fig_calibration.R` |
| Fig 3 — calibration geometry (4-panel) | `FigS_calibration_4panel.xlsx`$^\dagger$ | `experiments/simulations/calibration/out/calibration_secondary_metrics.csv` | `experiments/simulations/calibration/make_fig_calibration_4panel.R` |
| Fig 4 — cross-crop diagnostics | `Figure3_crosscrop.xlsx` | `tables/cross_crop_summary.csv`; panel A (matched share and null) from `experiments/empirical/multi/out/panel_workflow/panel_null_summary.csv` | `experiments/empirical/multi/make_fig_crosscrop.R` |
| Fig 5 — prediction | `Figure4_prediction.xlsx` | `experiments/empirical/multi/out/panel_workflow/raw/*__{pred,real}.csv` (miscanthus: `out/miscanthus_PA.csv`) | `experiments/empirical/multi/make_pred_figure.R` |
| Fig 6 — diagnostic workflow | *(none)* | hand-authored draw.io schematic (`figures/fig_diagnostic_workflow.drawio`), no plotted data | — |
| Supp — structure control | `FigS_structure_vd.xlsx` | `experiments/empirical/multi/out/*_structure_std.csv` | `experiments/empirical/multi/make_structure_figures.R` |
| Supp — regime map | `FigS_regime.xlsx` | `experiments/simulations/regime_map/out/regime_map_expanded.csv` | `experiments/simulations/regime_map/make_fig_regime_expanded.R` |

$^\dagger$ Workbook filenames still carry the pre-renumbering labels (`Figure3_crosscrop.xlsx` is now Fig 4, `Figure4_prediction.xlsx` is now Fig 5, `FigS_calibration_4panel.xlsx` is now Fig 3). The **Figure** column above is authoritative.

Sheet naming convention: `Panel*` sheets hold the exact plotted values
(means, SDs, fitted lines, heatmap/bar values); `Raw_*` sheets hold the
per-replicate / per-trait / per-cell rows they were aggregated from.

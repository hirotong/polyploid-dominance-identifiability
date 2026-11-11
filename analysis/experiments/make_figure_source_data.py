#!/usr/bin/env python3
# =====================================================================
# Collect the raw data behind each paper figure into one Excel file per
# figure. Each workbook reconstructs exactly the quantities plotted by the
# corresponding R figure script (raw underlying rows + the aggregated values
# that become points/tiles/bars in the figure).
#
# Figures covered (active in paper/main.tex + paper/appendix.tex):
#   Fig 1  fig_simulation              -> Figure1_simulation.xlsx
#   Fig 2  fig_calibration             -> Figure2_calibration.xlsx
#   Fig 3  fig_crosscrop               -> Figure3_crosscrop.xlsx
#   Fig 4  fig_pred                    -> Figure4_prediction.xlsx
#   Fig 5  fig_diagnostic_workflow     -> (AI-generated schematic, no data)
#   Supp   figS_calibration_4panel     -> FigS_calibration_4panel.xlsx
#   Supp   fig_structure_vd            -> FigS_structure_vd.xlsx
#   Supp   fig_regime                  -> FigS_regime.xlsx
#
# Run from repo root:  python3 experiments/make_figure_source_data.py
# =====================================================================
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "paper", "figure_source_data")
os.makedirs(OUT, exist_ok=True)


def p(*a):
    return os.path.join(ROOT, *a)


def write_xlsx(fname, sheets):
    """sheets: list of (sheet_name, DataFrame)."""
    path = os.path.join(OUT, fname)
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets:
            df.to_excel(xw, sheet_name=name[:31], index=False)
    print(f"wrote {fname}: {', '.join(n for n, _ in sheets)}")


STR_LV = ["unrelated", "halfsib", "fullsib"]
STR_LAB = {"unrelated": "Unrelated", "halfsib": "Half-sib", "fullsib": "Full-sib"}
ARCH_LV = ["additive", "digenic", "trigenic", "quadrigenic", "hetero"]
ARCH_LAB = {"additive": "Additive", "digenic": "Digenic", "trigenic": "Trigenic",
            "quadrigenic": "Quadrigenic", "hetero": "Heterozygosity"}


def ordcat(s, order):
    return pd.Categorical(s, categories=order, ordered=True)


# ===================== Figure 1: known-truth simulation =====================
def fig1():
    sv = pd.read_csv(p("experiments/simulations/poly_ortho/out/std_vd_share.csv"))
    oa = pd.read_csv(p("experiments/simulations/poly_ortho/out/ortho_attribution.csv"))

    # Panel A: dot = mean Vd_share, error bar = +/-1 SD, per structure x architecture
    a_summary = (sv.groupby(["structure", "truth"])["Vd_share"]
                 .agg(mean="mean", sd="std", n="count").reset_index())
    a_summary["truth"] = ordcat(a_summary["truth"], ARCH_LV)
    a_summary["structure"] = ordcat(a_summary["structure"], STR_LV)
    a_summary = a_summary.sort_values(["truth", "structure"]).reset_index(drop=True)
    a_summary.insert(1, "architecture", a_summary["truth"].map(ARCH_LAB))
    a_summary.insert(1, "structure_label", a_summary["structure"].map(STR_LAB))

    a_raw = sv.copy()

    # Panel B: full-sib heatmap, mean fitted fraction per (architecture x component)
    fs = oa[oa["structure"] == "fullsib"].copy()
    b_heat = (fs.groupby("truth")[["fA", "fD2", "fD3", "fD4"]].mean().reset_index()
              .rename(columns={"fA": "A", "fD2": "D2", "fD3": "D3", "fD4": "D4"}))
    b_heat["truth"] = ordcat(b_heat["truth"], ARCH_LV)
    b_heat = b_heat.sort_values("truth").reset_index(drop=True)
    b_heat.insert(1, "architecture", b_heat["truth"].map(ARCH_LAB))

    write_xlsx("Figure1_simulation.xlsx", [
        ("PanelA_summary", a_summary),
        ("PanelA_raw_replicates", a_raw),
        ("PanelB_heatmap_fullsib", b_heat),
        ("PanelB_raw_replicates", oa.copy()),
    ])


# ===================== Figure 2: calibration curves =========================
def fig2():
    d = pd.read_csv(p("experiments/simulations/calibration/out/calibration_primary.csv"))
    std = d[d["model"] == "A+D"][["population_structure", "replicate_id",
                                  "realized_true_D2_fraction", "standard_D_share"]].copy()
    std.columns = ["structure", "rep", "s", "est"]
    std["estimator"] = "standard"
    orth = d[d["model"] == "A+D2+D3+D4"][["population_structure", "replicate_id",
                                          "realized_true_D2_fraction", "orth_D2_fraction"]].copy()
    orth.columns = ["structure", "rep", "s", "est"]
    orth["estimator"] = "orthogonal"
    L = pd.concat([std, orth], ignore_index=True)
    L = L[np.isfinite(L["est"]) & np.isfinite(L["s"])]

    # per-cell summary (mean +/- SD) that becomes the plotted points/I-beams
    summ = (L.groupby(["estimator", "structure", "s"])["est"]
            .agg(mean="mean", sd="std", n="count").reset_index())
    summ["structure"] = ordcat(summ["structure"], STR_LV)
    summ = summ.sort_values(["estimator", "structure", "s"]).reset_index(drop=True)

    # pooled regression line per estimator (est = alpha + beta*s)
    rows = []
    for est, g in L.groupby("estimator"):
        beta, alpha = np.polyfit(g["s"], g["est"], 1)
        rows.append({"estimator": est, "alpha": alpha, "beta": beta, "n_obs": len(g)})
    pooled = pd.DataFrame(rows)

    metrics = pd.read_csv(p("experiments/simulations/calibration/out/calibration_metrics.csv"))

    write_xlsx("Figure2_calibration.xlsx", [
        ("Panel_summary_points", summ),
        ("Pooled_fit_lines", pooled),
        ("Calibration_metrics", metrics),
        ("Raw_per_replicate", L.reset_index(drop=True)),
    ])


# ===================== Figure 3: cross-crop diagnostics =====================
def fig3():
    S = pd.read_csv(p("tables/cross_crop_summary.csv"))
    lab_order = ["alfalfa", "blueberry_gc", "blueberry2", "blueberry", "miscanthus",
                 "panicum", "potato2", "potato", "rose"]
    lab_map = {"alfalfa": "Alfalfa", "blueberry_gc": "BB-UF-call", "blueberry2": "BB-UF-fruit",
               "blueberry": "BB-UF-GP", "miscanthus": "Miscanthus", "panicum": "Panicum",
               "potato2": "Potato-SLU", "potato": "Potato-SolCAP", "rose": "Rose"}
    S = S.set_index("name").loc[[k for k in lab_order if k in S["name"].values]].reset_index()
    S["display_label"] = S["name"].map(lab_map)
    # Panel A plots the matched comparison (same traits and same estimator as the
    # permutations), as in make_fig_crosscrop.R: override the summary columns from
    # the workflow's panel_null_summary.csv.
    NS = pd.read_csv(p("experiments/empirical/multi/out/panel_workflow/panel_null_summary.csv"))
    ns_key = dict(lab_map, blueberry_gc="BB-UF-call*")
    for i, key in enumerate(S["name"]):
        sub = NS[NS["panel"] == ns_key[key]]
        tr = sub[sub["trait_or_panel_summary"] != "PANEL_SUMMARY"]
        pl = sub[sub["trait_or_panel_summary"] == "PANEL_SUMMARY"]
        assert len(pl) == 1 and len(tr) >= 1, key
        S.loc[i, "std_Vd"] = pl["s_hat_AD_real"].iloc[0]
        S.loc[i, "perm_null"] = pl["perm_mean"].iloc[0]
        S.loc[i, "std_Vd_trait_min"] = tr["s_hat_AD_real"].min()
        S.loc[i, "std_Vd_trait_max"] = tr["s_hat_AD_real"].max()
        S.loc[i, "n_null_traits"] = len(tr)
    cols = ["name", "display_label", "clade", "n", "n_traits", "n_null_traits",
            "std_Vd", "std_Vd_trait_min", "std_Vd_trait_max", "perm_null",
            "AD_cos", "eff_rank",
            "ortho_A", "ortho_D2", "ortho_D3", "ortho_D4"]
    plotted = S[[c for c in cols if c in S.columns]].copy()
    write_xlsx("Figure3_crosscrop.xlsx", [
        ("Crosscrop_diagnostics", plotted),
    ])


# ===================== Figure 4: prediction =================================
def sd1(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return x.std(ddof=1) if x.size > 1 else np.nan


def fig4():
    RAW = "experiments/empirical/multi/out/panel_workflow/raw"
    OUTD = "experiments/empirical/multi/out"
    cfg = [("alfalfa", "Alfalfa", True, "tetra"),
           ("blueberry_gc", "BB-UF-call", False, "tetra"),
           ("blueberry2", "BB-UF-fruit", False, "tetra"),
           ("blueberry", "BB-UF-GP", False, "tetra"),
           ("miscanthus", "Miscanthus", False, "tetra"),
           ("panicum", "Panicum", False, "tetra"),
           ("potato2", "Potato-SLU", False, "tetra"),
           ("potato", "Potato-SolCAP", False, "tetra"),
           ("rose", "Rose", False, "tetra"),
           ("sugarcane", "Sugarcane (var.)", False, "higher"),
           ("sweetpotato", "Sweetpotato (6x)", False, "higher")]

    def raw_pa(key):
        f = p(RAW, f"{key}__pred.csv")
        if not os.path.exists(f):
            return None
        x = pd.read_csv(f)
        rows = []
        for tr, s in x.groupby("trait"):
            rows.append({
                "trait": tr,
                "PA_A": s["r_A"].mean(), "PA_A_sd": sd1(s["r_A"]),
                "PA_AD": s["r_AD"].mean(), "PA_AD_sd": sd1(s["r_AD"]),
                "PA_orth": s["r_orth"].mean(), "PA_orth_sd": sd1(s["r_orth"]),
                "d_AD_vs_A": s["delta_AD"].mean(), "d_AD_vs_A_sd": sd1(s["delta_AD"]),
                "d_orth_vs_AD": s["delta_orth"].mean(), "d_orth_vs_AD_sd": sd1(s["delta_orth"]),
            })
        return pd.DataFrame(rows)

    def pa_for(key):
        legacy = p(OUTD, f"{key}_PA.csv")
        pa = raw_pa(key)
        if pa is not None:
            return pa
        return pd.read_csv(legacy) if os.path.exists(legacy) else None

    def real_for(key):
        f = p(RAW, f"{key}__real.csv")
        if os.path.exists(f):
            r = pd.read_csv(f)
            if "s_hat_AD" in r.columns:
                return pd.DataFrame({"trait": r["trait"], "std_Vd": r["s_hat_AD"]})
        leg = p(OUTD, f"{key}_ortho.csv")
        return pd.read_csv(leg) if os.path.exists(leg) else None

    def se2(x, sdcol):
        x = np.asarray([v for v in x if np.isfinite(v)], float)
        if x.size > 1:
            return x.std(ddof=1) / np.sqrt(x.size)
        return np.nanmean(np.asarray(sdcol, float))

    # Panel A's x-axis in fig_pred is the scorecard's full-data diagnostic share
    # (s_hat_AD_real), the same source behind tab:pred / tab:workflow -- mirror
    # make_pred_figure.R exactly rather than re-deriving trait means here.
    score = pd.read_csv(p("experiments/empirical/multi/out/panel_workflow/workflow_scorecard.csv"))
    score_key = {"alfalfa": "Alfalfa", "blueberry_gc": "BB-UF-call*", "blueberry2": "BB-UF-fruit",
                 "blueberry": "BB-UF-GP", "miscanthus": "Miscanthus", "panicum": "Panicum",
                 "potato2": "Potato-SLU", "potato": "Potato-SolCAP", "rose": "Rose",
                 "sugarcane": "Sugarcane", "sweetpotato": "Sweetpotato"}

    def canon_sd(key):
        v = score.loc[score["panel"] == score_key.get(key, ""), "s_hat_AD_real"]
        return float(v.iloc[0]) if len(v) else np.nan

    panel_rows, pertrait = [], []
    for key, lab, ref, ploidy in cfg:
        pa = pa_for(key)
        if pa is None:
            continue
        vd = real_for(key)
        app_vd = canon_sd(key)
        gAD = pa["d_AD_vs_A"]
        gOR = pa["d_orth_vs_AD"]
        panel_rows.append({
            "key": key, "label": lab, "reference_panel": ref, "ploidy": ploidy,
            "n_traits": int(np.isfinite(gAD).sum()),
            "fitted_std_D_share": app_vd,
            "gain_AD_vs_A": np.nanmean(gAD), "gain_AD_vs_A_se": se2(gAD, pa["d_AD_vs_A_sd"]),
            "gain_orth_vs_AD": np.nanmean(gOR), "gain_orth_vs_AD_se": se2(gOR, pa["d_orth_vs_AD_sd"]),
        })
        pt = pa.copy()
        pt.insert(0, "panel", lab)
        pt.insert(0, "key", key)
        if vd is not None and "std_Vd" in vd:
            pt = pt.merge(vd.rename(columns={"std_Vd": "trait_std_D_share"}), on="trait", how="left")
        pertrait.append(pt)

    panel = pd.DataFrame(panel_rows)
    pertrait = pd.concat(pertrait, ignore_index=True)
    write_xlsx("Figure4_prediction.xlsx", [
        ("Panel_summary_plotted", panel),
        ("Per_trait_predictive_ability", pertrait),
    ])


# ===================== Supp: figS_calibration_4panel =======================
def figS_calib4():
    cells = pd.read_csv(p("experiments/simulations/calibration/out/calibration_secondary_metrics.csv"))
    # Plotted M set matches make_fig_calibration_4panel.R (M=1000 dropped there:
    # ill-conditioned A/D split at n=19200). Panel sheets carry plotted cells only;
    # the full grid incl. M=1000 stays in Raw_per_cell_metrics.
    Mvals = [5000, 20000, 50000, 100000, 200000]
    plot_cells = cells[cells["nMarker"].isin(Mvals)]
    # Panel A heatmap: mean no-dominance baseline (std_bias0) over LD, by n x M
    A = plot_cells.groupby(["N", "nMarker"])["std_bias0"].mean().reset_index(name="std_bias0_mean")
    # Panel B heatmap: mean calibration slope (std_beta) over LD, by n x M
    B = plot_cells.groupby(["N", "nMarker"])["std_beta"].mean().reset_index(name="std_beta_mean")
    # Panels C,D scatter: one point per (n,M,LD) cell with finite RMSE; M in Mvals
    CD = cells[np.isfinite(cells["std_rmse"]) & cells["nMarker"].isin(Mvals)][
        ["N", "nMarker", "LD", "LDlab", "adjLD_r2", "std_rmse", "AD_cos", "eff_rank"]].copy()
    write_xlsx("FigS_calibration_4panel.xlsx", [
        ("PanelA_baseline_heatmap", A),
        ("PanelB_slope_heatmap", B),
        ("PanelCD_scatter", CD.reset_index(drop=True)),
        ("Raw_per_cell_metrics", cells),
    ])


# ===================== Supp: fig_structure_vd ==============================
def figS_structvd():
    O = "experiments/empirical/multi/out"
    sets = ["alfalfa", "potato", "blueberry", "blueberry2", "potato2",
            "panicum", "rose", "miscanthus", "blueberry_gc"]
    labs = {"alfalfa": "Alfalfa", "potato": "Potato SolCAP", "blueberry": "Blueberry UF-GP",
            "blueberry2": "Blueberry UF-fruit", "potato2": "Potato SLU", "panicum": "Panicum",
            "rose": "Rose", "miscanthus": "Miscanthus", "blueberry_gc": "Blueberry UF-call"}
    panel_order = ["alfalfa", "blueberry_gc", "blueberry2", "blueberry",
                   "miscanthus", "panicum", "potato2", "potato", "rose"]
    rows, raw = [], []
    for n in sets:
        f = p(O, f"{n}_structure_std.csv")
        if not os.path.exists(f):
            continue
        s = pd.read_csv(f)
        if n == "alfalfa" and "trait" in s.columns:
            s = s[s["trait"] == "yi_year"]
        raw.append(s.assign(dataset_label=labs[n]))
        ag = s.groupby("control")["std_Vd"].mean()
        rows.append({"dataset": n, "dataset_label": labs[n],
                     "PC0_no_control": ag.get("PC0", np.nan),
                     "PC10_genomic_PCs": ag.get("PC10", np.nan),
                     "F_family_kernel": ag.get("F", np.nan)})
    V = pd.DataFrame(rows).set_index("dataset").loc[
        [k for k in panel_order if k in [r["dataset"] for r in rows]]].reset_index()
    write_xlsx("FigS_structure_vd.xlsx", [
        ("Structure_control_Vd_share", V),
        ("Raw_per_trait_control", pd.concat(raw, ignore_index=True)),
    ])


# ===================== Supp: fig_regime ====================================
def figS_regime():
    d = pd.read_csv(p("experiments/simulations/regime_map/out/regime_map_expanded.csv"))
    d = d[np.isfinite(d["spurious_Vd"]) & (d["failed"].isna() | (d["failed"] == 0))]
    fs = d[d["structure"] == "fullsib"].copy()
    # Panel A: heatmap mean +/- SD spurious share over n x M (full-sib)
    A = (fs.groupby(["N", "nMarker"])["spurious_Vd"]
         .agg(mean="mean", sd="std", n="count").reset_index())
    # Panel B: spurious share vs n, by LD level (full-sib)
    B = (fs.groupby(["LD", "N"])["spurious_Vd"]
         .agg(mean="mean", sd="std", n="count").reset_index())
    B["LD"] = ordcat(B["LD"], ["independent", "low", "medium", "high"])
    B = B.sort_values(["LD", "N"]).reset_index(drop=True)
    # realized-LD ladder (mean adjacent r2 by LD level)
    ladder = d.groupby("LD")["adjLD_r2"].mean().reset_index(name="mean_adjLD_r2")
    write_xlsx("FigS_regime.xlsx", [
        ("PanelA_heatmap_nxM", A),
        ("PanelB_spurious_vs_n_by_LD", B),
        ("Realized_LD_ladder", ladder),
        ("Raw_fullsib_cells", fs.reset_index(drop=True)),
    ])


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    figS_calib4()
    figS_structvd()
    figS_regime()
    print("\nAll figure source-data workbooks written to:", OUT)

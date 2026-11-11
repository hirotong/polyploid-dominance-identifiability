#!/usr/bin/env python3
"""Sensitivity of the workflow labels to the decision thresholds.

Re-applies the label rule of make_workflow_scorecard.R to workflow_scorecard.csv
with each threshold moved by -20% and +20% (one at a time, then all jointly), and
reports which panel labels change. No model is refitted.

Usage: python3 threshold_sensitivity.py [workflow_scorecard.csv] [out.csv]
Standard library only.
"""
import csv
import itertools
import sys

SCORE = sys.argv[1] if len(sys.argv) > 1 else "experiments/empirical/multi/out/panel_workflow/workflow_scorecard.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "experiments/empirical/multi/out/panel_workflow/threshold_sensitivity.csv"

# same names and values as THR in make_workflow_scorecard.R, plus the additive-null percentile
BASE = {"struct_large": 0.10, "recover_rate": 0.80, "ci_wide": 0.25,
        "null_margin": 0.05, "pred_gain": 0.01, "L_pctile": 0.95}
REL = 0.20


def num(x):
    try:
        v = float(x)
        return v if v == v else None       # drop NaN
    except (TypeError, ValueError):
        return None


def label(r, thr):
    S, ruler = num(r["S_struct_real"]), num(r["structure_ruler_q95_lambda0"])
    structure_sensitive = (S is not None and S >= thr["struct_large"]) or \
                          (S is not None and ruler is not None and S > ruler)
    dn = num(r["delta_null"])
    exceeds_null = dn is not None and dn > thr["null_margin"]
    rec = num(r["D2_recovery_rate_s050"])
    recoverable = rec is not None and rec >= thr["recover_rate"]
    P, pS, pN = num(r["P_D2_real"]), num(r["spike_median_P_D2_s050"]), num(r["add_median_P_D2"])
    Lp = num(r["real_L_D2_percentile_vs_additive_null"])
    resembles = (None not in (P, pS, pN)) and abs(P - pS) < abs(P - pN) and Lp is not None and Lp > thr["L_pctile"]
    W = num(r["D2_CrI_width_real"])
    narrow = W is not None and W < thr["ci_wide"]
    d, lo = num(r["mean_delta_AD"]), num(r["mean_delta_AD_ci_low"])
    pred = d is not None and d >= thr["pred_gain"] and (lo is None or lo > 0)
    if structure_sensitive:
        return "structure-sensitive"
    if exceeds_null and recoverable and resembles and narrow:
        return "D2-supported"
    return "prediction-useful" if pred else "model-basis"


rows = list(csv.DictReader(open(SCORE)))
# The scorecard rounds to four decimals, which can sit exactly on a moved threshold
# (e.g. a gain of 0.01198 printed as 0.012). Use the unrounded prediction summary
# written next to the scorecard when it is there.
import os
PRED = os.path.join(os.path.dirname(SCORE) or ".", "panel_prediction_summary.csv")
if os.path.exists(PRED):
    pr = {r["panel"]: r for r in csv.DictReader(open(PRED))}
    for r in rows:
        if r["panel"] in pr:
            r["mean_delta_AD"] = pr[r["panel"]]["mean_delta_AD"]
            r["mean_delta_AD_ci_low"] = pr[r["panel"]]["ci_low_AD"]
else:
    print(f"note: {PRED} not found; using rounded scorecard values for prediction gain")
base_lab = {r["panel"]: label(r, BASE) for r in rows}
for r in rows:                                # the re-implemented rule must reproduce the scorecard
    if r.get("suggested_label") and r["suggested_label"] != base_lab[r["panel"]]:
        sys.exit(f"rule mismatch for {r['panel']}: scorecard={r['suggested_label']} here={base_lab[r['panel']]}")

out = []
def run(name, thr):
    for r in rows:
        new = label(r, thr)
        out.append({"scenario": name, "panel": r["panel"], "base_label": base_lab[r["panel"]],
                    "label": new, "changed": int(new != base_lab[r["panel"]])})

for k in BASE:
    for sgn in (-1, +1):
        thr = dict(BASE); thr[k] = BASE[k] * (1 + sgn * REL)
        if k == "L_pctile":
            thr[k] = min(thr[k], 0.999)       # a percentile cannot exceed 1
        run(f"{k} {'+' if sgn > 0 else '-'}{int(REL * 100)}%", thr)
for signs in itertools.product((-1, +1), repeat=len(BASE)):   # all joint corners
    thr = {k: BASE[k] * (1 + s * REL) for k, s in zip(BASE, signs)}
    thr["L_pctile"] = min(thr["L_pctile"], 0.999)
    run("joint " + "".join("+" if s > 0 else "-" for s in signs), thr)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["scenario", "panel", "base_label", "label", "changed"])
    w.writeheader(); w.writerows(out)

single = [o for o in out if not o["scenario"].startswith("joint")]
joint = [o for o in out if o["scenario"].startswith("joint")]
print(f"panels: {len(rows)}; thresholds: {len(BASE)}; relative change: +/-{int(REL * 100)}%")
print("one-at-a-time changes:")
ch = [o for o in single if o["changed"]]
for o in ch:
    print(f"  {o['scenario']:22} {o['panel']:24} {o['base_label']} -> {o['label']}")
print(f"  total: {len(ch)} of {len(single)} panel-scenario cells; panels affected: {sorted({o['panel'] for o in ch})}")
chj = [o for o in joint if o["changed"]]
print(f"joint corners: {len({o['scenario'] for o in joint})} scenarios; panels whose label changes in at least one: "
      f"{sorted({o['panel'] for o in chj})}")
labs = {}
for o in joint + single:
    labs.setdefault(o["panel"], set()).add(o["label"])
for p, l in labs.items():
    print(f"  {p:24} labels seen: {sorted(l)}")
print(f"written: {OUT}")

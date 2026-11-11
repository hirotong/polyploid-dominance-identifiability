#!/usr/bin/env python3
"""Summarise estimator_sensitivity.R output (<panel>__estsens.csv) per panel.

Usage: python3 aggregate_estimator_sensitivity.py <RAW_DIR> [OUT_CSV]
  RAW_DIR  directory holding <panel>__estsens.csv
  OUT_CSV  default experiments/empirical/multi/out/estimator_sensitivity/estimator_sensitivity_summary.csv
Standard library only. Prints a readable table and writes one row per panel.

Columns (standard model, share s_D):
  bglr_default        mean over traits of the three-seed mean
  bglr_seed_sd        mean over traits of the SD across the three seeds (Monte Carlo noise)
  prior_min/max       smallest / largest trait-mean share over the four changed prior settings
  prior_max_shift     largest |changed setting - default| for any single trait
  reml_shares         REML share per trait, with 95% likelihood-ratio interval
  reml_ci_width_mean  mean interval width
  reml_at_boundary    number of traits with REML share at 0 or 1
  reml_perm_sd        mean over traits of the SD of the REML share on permuted phenotypes
Columns (orthogonal model): BGLR default D2 / D3 / D4 fractions, range of the D2
fraction over prior settings, and the REML D2 / D3 / D4 fractions per trait.
"""
import csv, glob, os, statistics as st, sys
from collections import defaultdict

RAW = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "experiments/empirical/multi/out/estimator_sensitivity/estimator_sensitivity_summary.csv"
ORDER = ["alfalfa", "potato", "potato2", "blueberry", "blueberry2", "blueberry_gc", "panicum", "rose",
         "miscanthus", "sweetpotato", "sugarcane", "poscontrol"]
CHANGED = ["R2=0.2", "R2=0.8", "df0=1", "df0=30"]
SEEDS = ["default_seed1", "default_seed2", "default_seed3"]

def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")

summary = []
for key in ORDER:
    fn = os.path.join(RAW, f"{key}__estsens.csv")
    if not os.path.exists(fn):
        continue
    rows = list(csv.DictReader(open(fn)))
    traits = list(dict.fromkeys(r["trait"] for r in rows))
    get = defaultdict(dict)                       # (trait, model, estimator, component) -> setting -> row
    for r in rows:
        get[(r["trait"], r["model"], r["estimator"], r["component"])][r["setting"]] = r
    def bg(tr, model, comp, setting):
        return f(get[(tr, model, "BGLR", comp)][setting]["value"])
    d_def, d_sd, shifts = [], [], []
    pri = {s: [] for s in CHANGED}
    for tr in traits:
        seeds = [bg(tr, "standard", "s_D", s) for s in SEEDS]
        d_def.append(st.mean(seeds)); d_sd.append(st.stdev(seeds))
        for s in CHANGED:
            v = bg(tr, "standard", "s_D", s); pri[s].append(v); shifts.append(abs(v - st.mean(seeds)))
    pri_means = {s: st.mean(v) for s, v in pri.items()}
    reml, widths, nbound, permsd, lr0 = [], [], 0, [], []
    for tr in traits:
        r = get[(tr, "standard", "REML", "s_D")]["profile"]
        v, lo, hi = f(r["value"]), f(r["lo"]), f(r["hi"])
        reml.append(f"{v:.2f} [{lo:.2f}, {hi:.2f}]"); widths.append(hi - lo); nbound += int(v <= 0 or v >= 1)
        lr0.append(f(get[(tr, "standard", "REML", "LR_vs_s0")]["profile"]["value"]))
        pr = [x for k, x in get[(tr, "standard", "REML", "s_D_perm_mean")].items()]
        permsd.append(f(pr[0]["extra"]))
    comps = [c for c in dict.fromkeys(r["component"] for r in rows if r["model"] == "orthogonal" and r["estimator"] == "BGLR")]
    nonadd = [c for c in comps if c != "A"]
    o_def = {c: st.mean(st.mean(bg(tr, "orthogonal", c, s) for s in SEEDS) for tr in traits) for c in comps}
    o_d2_pri = [st.mean(bg(tr, "orthogonal", "D2", s) for tr in traits) for s in CHANGED]
    o_reml = []
    for tr in traits:
        g = get.get((tr, "orthogonal", "REML", "D2"))
        if g:
            o_reml.append("/".join(f"{f(get[(tr, 'orthogonal', 'REML', c)]['optim']['value']):.2f}" for c in nonadd[:3]))
    row = dict(panel=rows[0]["panel"], panel_key=key, n=rows[0]["n"], n_traits=len(traits),
               bglr_default=st.mean(d_def), bglr_seed_sd=st.mean(d_sd),
               prior_min=min(pri_means.values()), prior_max=max(pri_means.values()),
               prior_max_shift=max(shifts),
               reml_shares="; ".join(reml), reml_ci_width_mean=st.mean(widths),
               reml_n_ci_whole_range=sum(w >= 0.999 for w in widths),
               reml_at_boundary=nbound, reml_lr_vs_additive_max=max(lr0), reml_perm_sd=st.mean(permsd),
               orth_bglr_D2=o_def.get("D2"), orth_bglr_D3=o_def.get("D3"), orth_bglr_D4=o_def.get("D4"),
               orth_bglr_D2_prior_min=min(o_d2_pri), orth_bglr_D2_prior_max=max(o_d2_pri),
               orth_reml_D2_D3_D4="; ".join(o_reml))
    summary.append(row)

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader()
    for r in summary:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

print(f"{'panel':22}{'BGLR':>6}{'seedSD':>8}{'prior range':>14}{'maxshift':>9}  REML share [95% LR interval] per trait")
for r in summary:
    print(f"{r['panel']:22}{r['bglr_default']:6.2f}{r['bglr_seed_sd']:8.3f}{r['prior_min']:8.2f}-{r['prior_max']:.2f}"
          f"{r['prior_max_shift']:9.2f}  {r['reml_shares']}   permSD={r['reml_perm_sd']:.2f}")
print()
print(f"{'panel':22}{'BGLR D2/D3/D4':>18}{'D2 prior range':>16}  REML D2/D3/D4 per trait")
for r in summary:
    d4 = r['orth_bglr_D4'] if r['orth_bglr_D4'] is not None else float('nan')
    print(f"{r['panel']:22}{r['orth_bglr_D2']:6.2f}/{r['orth_bglr_D3']:.2f}/{d4:.2f}"
          f"{r['orth_bglr_D2_prior_min']:10.2f}-{r['orth_bglr_D2_prior_max']:.2f}  {r['orth_reml_D2_D3_D4']}")
print(f"\nwritten: {OUT}")


# ---- LaTeX rows for the supplementary table + the counts quoted in the text ----
TEXNAME = {"alfalfa": "Alfalfa", "potato": "Potato (SolCAP)", "potato2": "Potato (SLU)",
           "blueberry": "Blueberry (UF-GP)", "blueberry2": "Blueberry (UF-fruit)",
           "blueberry_gc": "Blueberry (UF-call$^\\dagger$)", "panicum": "\\emph{Panicum}", "rose": "Rose",
           "miscanthus": "\\emph{Miscanthus}", "sweetpotato": "Sweetpotato", "sugarcane": "Sugarcane",
           "poscontrol": "Positive control (sim)"}
TEXORDER = ["alfalfa", "blueberry_gc", "blueberry2", "blueberry", "miscanthus", "panicum", "potato2", "potato",
            "rose", "sugarcane", "sweetpotato", "poscontrol"]
tex, stats = [], dict(n_tr=0, whole=0, excl0=0, excl1=0, seed_sd_max=0.0, shift_max=0.0, perm_bound=[], perm_sd=[],
                      high_small=0, d2_prior_shift_max=0.0)
bykey = {r["panel_key"]: r for r in summary}
for key in TEXORDER:
    if key not in bykey:
        continue
    rows = list(csv.DictReader(open(os.path.join(RAW, f"{key}__estsens.csv"))))
    traits = list(dict.fromkeys(r["trait"] for r in rows))
    val = {(r["trait"], r["model"], r["estimator"], r["setting"], r["component"]): r for r in rows}
    comps = list(dict.fromkeys(r["component"] for r in rows if r["model"] == "orthogonal" and r["estimator"] == "BGLR"))
    high = [c for c in comps if c not in ("A", "D2")]
    reml_s, reml_o = [], []
    for tr in traits:
        r = val[(tr, "standard", "REML", "profile", "s_D")]; v, lo, hi = f(r["value"]), f(r["lo"]), f(r["hi"])
        reml_s.append(f"{v:.2f} [{lo:.2f}, {hi:.2f}]")
        d2 = f(val[(tr, "orthogonal", "REML", "optim", "D2")]["value"])
        hs = sum(f(val[(tr, "orthogonal", "REML", "optim", c)]["value"]) for c in high)
        reml_o.append(f"{d2:.2f} / {hs:.2f}")
        seeds = [f(val[(tr, "standard", "BGLR", s, "s_D")]["value"]) for s in SEEDS]
        d2seeds = st.mean(f(val[(tr, "orthogonal", "BGLR", s, "D2")]["value"]) for s in SEEDS)
        if key != "poscontrol":
            stats["n_tr"] += 1; stats["whole"] += int(lo <= 0 and hi >= 1); stats["excl0"] += int(lo > 0); stats["excl1"] += int(hi < 1)
            stats["seed_sd_max"] = max(stats["seed_sd_max"], st.stdev(seeds))
            stats["high_small"] += int(hs < 0.05)
            for s in CHANGED:
                stats["shift_max"] = max(stats["shift_max"], abs(f(val[(tr, "standard", "BGLR", s, "s_D")]["value"]) - st.mean(seeds)))
                stats["d2_prior_shift_max"] = max(stats["d2_prior_shift_max"], abs(f(val[(tr, "orthogonal", "BGLR", s, "D2")]["value"]) - d2seeds))
            pk = [k for k in val if k[0] == tr and k[4] == "s_D_perm_at_boundary"][0]
            stats["perm_bound"].append(f(val[pk]["value"]))
            pk = [k for k in val if k[0] == tr and k[4] == "s_D_perm_mean"][0]
            stats["perm_sd"].append(f(val[pk]["extra"]))
    b = bykey[key]
    d4 = b["orth_bglr_D4"]
    tex.append(f"      {TEXNAME[key]} & {b['bglr_default']:.2f} & {b['prior_min']:.2f}--{b['prior_max']:.2f} & "
               f"{'; '.join(reml_s)} & {b['orth_bglr_D2']:.2f} / {b['orth_bglr_D3']:.2f} / {d4:.2f} & "
               f"{b['orth_bglr_D2_prior_min']:.2f}--{b['orth_bglr_D2_prior_max']:.2f} & {'; '.join(reml_o)} \\\\")
texfile = os.path.join(os.path.dirname(OUT) or ".", "estimator_sensitivity_rows.tex")
open(texfile, "w").write("\n".join(tex) + "\n")
print("\n".join(tex))
print("\ncounts over empirical traits:")
print(f"  traits={stats['n_tr']} REML interval whole range={stats['whole']} excludes 0={stats['excl0']} excludes 1={stats['excl1']}")
print(f"  max per-trait BGLR seed SD={stats['seed_sd_max']:.4f}  max per-trait prior shift (s_D)={stats['shift_max']:.4f}  (D2)={stats['d2_prior_shift_max']:.4f}")
print(f"  REML permutation share at 0 or 1: mean={st.mean(stats['perm_bound']):.3f}  REML perm SD range={min(stats['perm_sd']):.2f}-{max(stats['perm_sd']):.2f}")
print(f"  traits with REML higher-degree (>=3) fraction below 0.05: {stats['high_small']} of {stats['n_tr']}")
print(f"written: {texfile}")

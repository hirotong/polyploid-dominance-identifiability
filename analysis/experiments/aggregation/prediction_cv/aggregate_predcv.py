#!/usr/bin/env python3
"""Aggregate per-panel prediction-CV outputs into canonical experiments/aggregation/prediction_cv tables.

Inputs:
  experiments/empirical/multi/out/panel_workflow/raw/*__pred.csv  (primary current output)
  experiments/empirical/multi/out/*_PA.csv                  (fallback historical output)
  experiments/simulations/poly_ortho/out/real_ortho_roza_PA_paired.csv  (fallback alfalfa CV)

Outputs:
  experiments/aggregation/prediction_cv/predcv_combined_alltraits.csv
  experiments/aggregation/prediction_cv/predcv_perdataset_summary.csv

The primary raw files store one paired prediction result per repeat and trait.
Rows are aggregated to trait-level means and repeat-level standard deviations.
"""

import csv
from pathlib import Path

OUT = Path("experiments/empirical/multi/out")
RAW = Path("experiments/empirical/multi/out/panel_workflow/raw")
STAGE = Path("experiments/aggregation/prediction_cv")
POLY = Path("experiments/simulations/poly_ortho/out/real_ortho_roza_PA_paired.csv")

ORDER = [
    "alfalfa",
    "potato",
    "potato2",
    "blueberry",
    "panicum",
    "rose",
    "miscanthus",
    "blueberry_gc",
    "blueberry2",
    "sweetpotato",
    "sugarcane",
]

PLOIDY = {
    "alfalfa": "4",
    "potato": "4",
    "potato2": "4",
    "blueberry": "4",
    "blueberry_gc": "4",
    "blueberry2": "4",
    "panicum": "4",
    "rose": "4",
    "miscanthus": "4",
    "sweetpotato": "6",
    "sugarcane": "var",
}

NDEG = {
    "alfalfa": "4",
    "potato": "4",
    "potato2": "4",
    "blueberry": "4",
    "blueberry_gc": "4",
    "blueberry2": "4",
    "panicum": "4",
    "rose": "4",
    "miscanthus": "4",
    "sweetpotato": "6",
    "sugarcane": "6",
}

COMBINED_FIELDS = [
    "dataset", "ploidy", "trait", "n", "ndeg",
    "PA_A", "PA_A_sd", "PA_AD", "PA_AD_sd", "PA_orth", "PA_orth_sd",
    "d_AD_vs_A", "d_AD_vs_A_sd", "d_orth_vs_AD", "d_orth_vs_AD_sd",
]

SUMMARY_FIELDS = [
    "dataset", "ploidy", "n_traits", "A", "AD", "orth",
    "AD_vs_A", "AD_wins", "orth_vs_AD", "orth_wins", "orth_win_1sd",
]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def stdev(xs):
    xs = [x for x in xs if x is not None]
    if len(xs) < 2:
        return None
    mu = sum(xs) / len(xs)
    return (sum((x - mu) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5


def fmt(x, digits=3):
    if x is None:
        return ""
    if abs(x) < 0.5 * 10 ** (-digits):
        x = 0.0
    return f"{x:.{digits}f}".rstrip("0").rstrip(".")


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def alfalfa_rows():
    if not POLY.exists():
        return []
    rows = []
    for r in read_csv(POLY):
        trait = r["trait"]
        pa_a = f(r["PA_additive"])
        pa_ad = f(r["PA_std_AD"])
        pa_o = f(r["PA_orthogonal"])
        rows.append({
            "dataset": "alfalfa",
            "ploidy": "4",
            "trait": trait,
            "n": "",
            "ndeg": "4",
            "PA_A": fmt(pa_a, 4),
            "PA_A_sd": "",
            "PA_AD": fmt(pa_ad, 4),
            "PA_AD_sd": "",
            "PA_orth": fmt(pa_o, 4),
            "PA_orth_sd": "",
            "d_AD_vs_A": fmt(pa_ad - pa_a, 4),
            "d_AD_vs_A_sd": "",
            "d_orth_vs_AD": fmt(pa_o - pa_ad, 4),
            "d_orth_vs_AD_sd": "",
        })
    return rows


def raw_panel_rows(key):
    path = RAW / f"{key}__pred.csv"
    if not path.exists():
        return []

    grouped = {}
    for row in read_csv(path):
        grouped.setdefault(row["trait"], []).append(row)

    rows = []
    for trait, trait_rows in grouped.items():
        pa_a = [f(r.get("r_A")) for r in trait_rows]
        pa_ad = [f(r.get("r_AD")) for r in trait_rows]
        pa_o = [f(r.get("r_orth")) for r in trait_rows]
        d_ad = [f(r.get("delta_AD")) for r in trait_rows]
        d_o = [f(r.get("delta_orth")) for r in trait_rows]
        rows.append({
            "dataset": key,
            "ploidy": PLOIDY.get(key, ""),
            "trait": trait,
            "n": "",
            "ndeg": NDEG.get(key, ""),
            "PA_A": fmt(mean(pa_a), 4),
            "PA_A_sd": fmt(stdev(pa_a), 4),
            "PA_AD": fmt(mean(pa_ad), 4),
            "PA_AD_sd": fmt(stdev(pa_ad), 4),
            "PA_orth": fmt(mean(pa_o), 4),
            "PA_orth_sd": fmt(stdev(pa_o), 4),
            "d_AD_vs_A": fmt(mean(d_ad), 4),
            "d_AD_vs_A_sd": fmt(stdev(d_ad), 4),
            "d_orth_vs_AD": fmt(mean(d_o), 4),
            "d_orth_vs_AD_sd": fmt(stdev(d_o), 4),
        })
    return rows


def panel_rows(key):
    rows = raw_panel_rows(key)
    if rows:
        return rows
    if key == "alfalfa":
        return alfalfa_rows()
    path = OUT / f"{key}_PA.csv"
    if not path.exists():
        return []
    rows = []
    for r in read_csv(path):
        rows.append({field: r.get(field, "") for field in COMBINED_FIELDS})
    return rows


def summarize(key, rows):
    a = [f(r["PA_A"]) for r in rows]
    ad = [f(r["PA_AD"]) for r in rows]
    orth = [f(r["PA_orth"]) for r in rows]
    d_ad = [f(r["d_AD_vs_A"]) for r in rows]
    d_orth = [f(r["d_orth_vs_AD"]) for r in rows]
    orth_sd = [f(r["d_orth_vs_AD_sd"]) for r in rows]
    n = len(rows)
    orth_win_1sd = sum(
        1 for d, sd in zip(d_orth, orth_sd)
        if d is not None and sd is not None and d > sd
    )
    return {
        "dataset": key,
        "ploidy": rows[0]["ploidy"] if rows else "",
        "n_traits": str(n),
        "A": fmt(mean(a)),
        "AD": fmt(mean(ad)),
        "orth": fmt(mean(orth)),
        "AD_vs_A": fmt(mean(d_ad)),
        "AD_wins": f"{sum(x is not None and x > 0 for x in d_ad)}/{n}",
        "orth_vs_AD": fmt(mean(d_orth)),
        "orth_wins": f"{sum(x is not None and x > 0 for x in d_orth)}/{n}",
        "orth_win_1sd": f"{orth_win_1sd}/{n}",
    }


def main():
    STAGE.mkdir(exist_ok=True)
    combined = []
    grouped = {}
    for key in ORDER:
        rows = panel_rows(key)
        if not rows:
            continue
        grouped[key] = rows
        combined.extend(rows)

    with (STAGE / "predcv_combined_alltraits.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMBINED_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(combined)

    summary = [summarize(key, grouped[key]) for key in ORDER if key in grouped]
    with (STAGE / "predcv_perdataset_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)

    non_alfalfa = [r for r in combined if r["dataset"] != "alfalfa"]
    d_ad = [f(r["d_AD_vs_A"]) for r in non_alfalfa]
    d_orth = [f(r["d_orth_vs_AD"]) for r in non_alfalfa]
    print(f"wrote {len(combined)} total trait rows across {len(grouped)} panels")
    print(
        "non-alfalfa: "
        f"{len(non_alfalfa)} traits, "
        f"A+D wins {sum(x is not None and x > 0 for x in d_ad)}/{len(d_ad)}, "
        f"mean d_AD={fmt(mean(d_ad), 4)}, mean d_orth={fmt(mean(d_orth), 4)}"
    )


if __name__ == "__main__":
    main()

"""Find stars that faded or dropped out of NEOWISE visit light curves.

Usage: python faders.py EPOCHS.csv OUT.csv [--survey-end MJD]

Two selections are reported for sources with at least eight W1 visits.

fade     the median of the last three visits is below half the median of the
         first five visits in W1 (and in W2 when W2 is measured), and each of
         the last three visits is below 0.6 of the early level.
dropout  the source was measured in at least eight visits at a median W1
         brighter than 14.5 mag, and its last measured visit is more than 300
         days before the end of the survey, so at least one full sky pass is
         missing.

Single-exposure catalogs hold detections only, so a star that fades below
about W1 = 15.5 leaves no row. Dropouts therefore need image-level checks.
"""
import argparse
import csv

import numpy as np
import pandas as pd

W1_ZERO_MJY = 309.54e3


def examine(group, survey_end):
    group = group.dropna(subset=["f1"]).sort_values("mjd")
    n = len(group)
    if n < 8:
        return None
    mjd, f1, f2 = group.mjd.values, group.f1.values, group.f2.values
    early = np.median(f1[:5])
    late = np.median(f1[-3:])
    row = {"n_visits": n, "first_mjd": mjd[0], "last_mjd": mjd[-1],
           "early_w1_mag": -2.5 * np.log10(early / W1_ZERO_MJY), "late_over_early_w1": late / early,
           "peak_over_early_w1": float(np.max(f1) / early), "peak_mjd": float(mjd[np.argmax(f1)]),
           "late_over_early_w2": np.nan, "kind": ""}
    have2 = np.isfinite(f2)
    if have2[:5].sum() >= 3 and have2[-3:].sum() >= 2:
        row["late_over_early_w2"] = float(np.nanmedian(f2[-3:]) / np.nanmedian(f2[:5]))
    fade = late < 0.5 * early and np.all(f1[-3:] < 0.6 * early) and not (row["late_over_early_w2"] > 0.6)
    dropout = (survey_end - mjd[-1] > 300) and row["early_w1_mag"] < 14.5
    if fade:
        row["kind"] = "fade"
    if dropout:
        row["kind"] = (row["kind"] + "+dropout").strip("+")
    return row if row["kind"] else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("epochs")
    parser.add_argument("out")
    parser.add_argument("--survey-end", type=float, default=60520.0)
    args = parser.parse_args()
    table = pd.read_csv(args.epochs, dtype={"id": str})
    rows, examined = [], 0
    for name, group in table.groupby("id", sort=False):
        examined += group.f1.notna().sum() >= 8
        result = examine(group, args.survey_end)
        if result:
            rows.append({"id": name, **result})
    rows.sort(key=lambda r: r["late_over_early_w1"])
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id"] + list(rows[0].keys())[1:] if rows else ["id"])
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in row.items()})
    print("sources with >=8 W1 visits:", examined, "selected:", len(rows))


if __name__ == "__main__":
    main()

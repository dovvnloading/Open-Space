"""Light curves of the unreported flares that pass every cut, most luminous first.

Cuts: class "unreported", one peak, at least two visits after onset, baseline
scatter at 4.6 um of 3% or less, a rise of at least 0.4 mag at 4.6 um and a
peak blackbody luminosity of at least 3e42 erg/s.

Usage: python figure.py SUMMARY.csv EPOCHS.csv OUT.png [--count N]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def prime(summary):
    keep = ((summary["class"] == "unreported") & (summary.n_peaks == 1) & (summary.visits_after_onset >= 2)
            & (summary.baseline_scatter_w2 <= 0.03) & (summary.peak_rise_w2_mag >= 0.4)
            & (summary.excess_bb_luminosity_erg_s >= 3e42))
    return summary[keep].sort_values("excess_bb_luminosity_erg_s", ascending=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary")
    parser.add_argument("epochs")
    parser.add_argument("out")
    parser.add_argument("--count", type=int, default=12)
    args = parser.parse_args()
    summary = pd.read_csv(args.summary, dtype={"id": str})
    chosen = prime(summary).head(args.count)
    epochs = pd.read_csv(args.epochs, dtype={"id": str})
    columns = 3
    rows = int(np.ceil(len(chosen) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(13.5, 3.3 * rows), squeeze=False)
    for axis, (_, row) in zip(axes.flat, chosen.iterrows()):
        visits = epochs[epochs.id == row.id].sort_values("mjd")
        year = 2000 + (visits.mjd - 51544.5) / 365.25
        axis.errorbar(year, visits.f1, yerr=visits.s1, fmt="o-", ms=3.5, lw=1, color="#1f5fa8", label="3.4 um")
        axis.errorbar(year, visits.f2, yerr=visits.s2, fmt="s-", ms=3.5, lw=1, color="#c0392b", label="4.6 um")
        axis.axvline(row.last_quiet_year, color="0.6", lw=0.8, ls=":")
        axis.set_ylim(bottom=0)
        name = {"Z 539-25": "CGCG 539-025", "MCG+06-08-012": "UGC 2653"}.get(row.host, row.host)
        axis.set_title(f"{name}   {row.distance_mpc:.0f} Mpc   "
                       f"{row.excess_bb_luminosity_erg_s:.1e} erg/s", fontsize=9)
        axis.set_ylabel("nuclear flux (mJy)", fontsize=8)
        axis.tick_params(labelsize=8)
    axes.flat[0].legend(fontsize=7, loc="upper left")
    for axis in axes.flat[len(chosen):]:
        axis.set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=110)
    print(len(chosen), "panels")


if __name__ == "__main__":
    main()

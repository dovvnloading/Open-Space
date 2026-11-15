"""Find mid-infrared brightenings after a quiet baseline in visit light curves.

Usage: python flares.py EPOCHS.csv OUT.csv

For every source with at least 12 visits in both bands, each visit k >= 6 is
tested as an onset. The baseline is every earlier visit: level = median,
scatter = max(1.4826 MAD, median visit uncertainty). An onset needs
  * a baseline scatter below 10% of the level in both bands,
  * visit k more than 5 sigma and 10% above the level in both W1 and W2,
  * the next visit (when one exists) more than 3 sigma above in both bands.
The earliest qualifying onset is reported with the peak excess that follows
and the number of separate maxima in the W2 excess.
"""
import csv
import sys

import numpy as np
import pandas as pd


def mad_sigma(values):
    return 1.4826 * np.median(np.abs(values - np.median(values)))


def count_peaks(excess):
    """Number of separate maxima above 30% of the highest, split by a dip below
    60% of the lower neighbouring maximum."""
    top = excess.max()
    peaks = []
    for i, value in enumerate(excess):
        left = excess[i - 1] if i > 0 else 0.0
        right = excess[i + 1] if i + 1 < len(excess) else -np.inf
        if value > 0.3 * top and value > left and value >= right:
            if peaks and excess[peaks[-1]:i + 1].min() > 0.6 * min(excess[peaks[-1]], value):
                if value > excess[peaks[-1]]:
                    peaks[-1] = i
                continue
            peaks.append(i)
    return len(peaks)


def search(group):
    group = group.dropna(subset=["f1", "f2", "s1", "s2"]).sort_values("mjd")
    n = len(group)
    if n < 12:
        return None
    mjd = group.mjd.values
    flux = {b: group[f"f{b}"].values for b in "12"}
    err = {b: group[f"s{b}"].values for b in "12"}
    for k in range(6, n):
        level, sigma, z = {}, {}, {}
        quiet = True
        for b in "12":
            level[b] = np.median(flux[b][:k])
            sigma[b] = max(mad_sigma(flux[b][:k]), np.median(err[b][:k]))
            quiet &= sigma[b] < 0.10 * level[b]
            z[b] = (flux[b][k:] - level[b]) / np.hypot(sigma[b], err[b][k:])
        if not quiet:
            continue
        rise = all(z[b][0] > 5 and flux[b][k] > 1.10 * level[b] for b in "12")
        if not rise:
            continue
        if k + 1 < n and not all(z[b][1] > 3 for b in "12"):
            continue
        excess2 = flux["2"][k:] - level["2"]
        peak = int(np.argmax(excess2))
        excess1 = flux["1"][k + peak] - level["1"]
        return {
            "n_visits": n, "level_w1_mjy": level["1"], "level_w2_mjy": level["2"],
            "base_scatter_w1": sigma["1"] / level["1"], "base_scatter_w2": sigma["2"] / level["2"],
            "last_quiet_mjd": mjd[k - 1], "onset_mjd": mjd[k], "peak_mjd": mjd[k + peak],
            "last_mjd": mjd[-1], "visits_after_onset": n - k,
            "onset_z_w1": z["1"][0], "onset_z_w2": z["2"][0],
            "peak_excess_w1_mjy": excess1, "peak_excess_w2_mjy": excess2[peak],
            "peak_dmag_w1": 2.5 * np.log10(flux["1"][k + peak] / level["1"]),
            "peak_dmag_w2": 2.5 * np.log10(flux["2"][k + peak] / level["2"]),
            "excess_color_w1_w2": (2.5 * np.log10(excess2[peak] / excess1 * 309.54 / 171.787)
                                   if excess1 > 0 else np.nan),
            "last_excess_fraction": excess2[-1] / excess2[peak], "n_peaks": count_peaks(excess2),
            "quiescent_w1_w2": 2.5 * np.log10(level["2"] / level["1"] * 309.54 / 171.787),
        }
    return None


def main():
    table = pd.read_csv(sys.argv[1], dtype={"id": str})
    rows = []
    searched = 0
    for name, group in table.groupby("id", sort=False):
        result = search(group)
        searched += result is not None or group.dropna(subset=["f1", "f2"]).shape[0] >= 12
        if result:
            rows.append({"id": name, **result})
    rows.sort(key=lambda r: -r["onset_mjd"])
    with open(sys.argv[2], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in row.items()})
    print("sources with >=12 two-band visits:", searched, "flares:", len(rows))


if __name__ == "__main__":
    main()

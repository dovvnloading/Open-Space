"""Search the March 2025 OGLE-IV SMC eclipsing binary release for changing eclipse depths."""

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np

EARLY_END = 2457000
LATE_START = 2458300
CENTERS = (-0.03, -0.01, 0.01, 0.03)


def observations(path):
    data = np.loadtxt(path, ndmin=2)
    good = np.all(np.isfinite(data), axis=1) & (data[:, 2] > 0) & (data[:, 2] <= 0.10)
    return data[good].T


def corrected_curve(data, epoch, period):
    t, mag, err = data
    phase = (t - epoch) / period
    phase = (phase + 0.5) % 1.0 - 0.5
    keep = np.zeros(len(t), dtype=bool)
    residual = np.zeros(len(t))
    season = np.floor((t - 2455200) / 365.25)
    quadrature = (abs(phase) > 0.16) & (abs(phase) < 0.34)
    for year in np.unique(season):
        group = season == year
        ref = mag[group & quadrature]
        if len(ref) >= 6:
            residual[group] = mag[group] - np.median(ref)
            keep |= group
    return t[keep], phase[keep], residual[keep], err[keep]


def difference(t, phase, y, err, center):
    offsets = (phase - center + 0.5) % 1 - 0.5
    values, weights = [], []
    for mid in CENTERS:
        window = (offsets >= mid - 0.01) & (offsets < mid + 0.01)
        a, b = window & (t < EARLY_END), window & (t >= LATE_START)
        if min(a.sum(), b.sum()) < 3:
            continue
        values.append(np.median(y[b]) - np.median(y[a]))
        var = 0.0
        for subset in (a, b):
            series = y[subset]
            scatter = 1.4826 * np.median(abs(series - np.median(series)))
            var += 1.57 * max(scatter, np.median(err[subset]), 0.01) ** 2 / len(series)
        weights.append(1 / var)
    if len(values) < 2:
        return None
    return float(np.average(values, weights=weights)), float(1 / np.sqrt(sum(weights)))


def screen(raw, output):
    ids = {
        line[:18]: line[20:23].strip()
        for line in (raw / "ident.dat").read_text().splitlines()
    }
    counts = Counter()
    rows = []
    for line in (raw / "ecl.dat").read_text().splitlines():
        star, _, _, per, ep, primary, _secondary = line.split()
        counts["catalog"] += 1
        if (
            "-" in (per, ep, primary)
            or ids[star] != "NC"
            or not 0.5 <= float(per) <= 12
            or float(primary) < 0.12
        ):
            continue
        counts["catalog_eligible"] += 1
        path = raw / "phot/I" / (star + ".dat")
        if not path.exists():
            continue
        counts["with_I"] += 1
        t, phase, y, err = corrected_curve(observations(path), float(ep), float(per))
        p = difference(t, phase, y, err, 0.0)
        s = difference(t, phase, y, err, 0.5)
        if p is None or s is None:
            continue
        counts["comparable"] += 1
        score = min(abs(p[0] / p[1]), abs(s[0] / s[1])) if p[0] * s[0] > 0 else 0.0
        rows.append(
            (
                star,
                score,
                p[0],
                p[1],
                s[0],
                s[1],
                np.sum(t < EARLY_END),
                np.sum(t >= LATE_START),
            )
        )
    rows.sort(key=lambda row: (-row[1], row[0]))
    with output.open("w", newline="") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(
            [
                "id",
                "ranking_score",
                "primary_delta_mag",
                "primary_nominal_error",
                "secondary_delta_mag",
                "secondary_nominal_error",
                "early_points",
                "late_points",
            ]
        )
        for row in rows:
            w.writerow([row[0], *(round(v, 6) for v in row[1:6]), *row[6:]])
    print(dict(counts))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    screen(args.raw, args.output)

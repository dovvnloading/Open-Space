"""Scan the decade-long KMT difference-flux series at five positions for multiday brightenings."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze import HERE, POSITIONS


def scan(path):
    rows = pd.DataFrame(json.loads(path.read_bytes())["rows"])
    good = rows[(rows.err > 0) & (rows.secz < 2) & rows.see.between(3, 7)].copy()
    good["night"] = np.floor(good.date).astype(int)
    daily = good.groupby(["year", "night", "site"]).flux.agg(["count", "median"])
    daily = daily[daily["count"] >= 3]["median"].unstack()
    peaks = []
    for year, season in daily.groupby(level="year"):
        season = season.droplevel("year")
        # A complete day grid makes the widths independent of weather and season gaps.
        season = season.reindex(range(season.index.min(), season.index.max() + 1))
        times = season.index.to_numpy()
        data = season[["BLG", "SAO", "SSO"]].to_numpy()
        matches = []
        for i in np.where(np.isfinite(data).any(axis=1))[0]:
            on = data[max(0, i-5):i+6]
            baseline = np.vstack([data[max(0, i-40):max(0, i-14)], data[i+15:i+41]])
            if min(np.isfinite(on).sum(axis=0)) < 3 or min(np.isfinite(baseline).sum(axis=0)) < 12:
                continue
            median = np.nanmedian(baseline, axis=0)
            scatter = 1.4826 * np.nanmedian(abs(baseline - median), axis=0)
            if min(scatter) <= 0:
                continue
            z = (np.nanmedian(on, axis=0) - median) / scatter
            if min(z) > 2.5 and max(z) > 5:
                matches.append((int(times[i]), z))
        for cluster in np.split(np.arange(len(matches)), np.where(np.diff([t for t, _ in matches]) > 10)[0]+1):
            if len(cluster):
                t, z = max((matches[i] for i in cluster), key=lambda v: min(v[1]))
                peaks.append((int(year), t, len(cluster), *z))
    return peaks


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--output", type=Path, default=HERE)
    args = p.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    rows = []
    for name, lcid in POSITIONS.items():
        rows.extend((name, *r) for r in scan(args.cache / f"{lcid}.json"))
    pd.DataFrame(rows, columns=["position", "year", "center", "n_passing_centers",
                               "Chile", "South_Africa", "Australia"]).to_csv(
                                   args.output / "decade_scan.csv", index=False, float_format="%.3f")

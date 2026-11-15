"""Recompute late Swift XRT/BAT fluence ratios from the frozen sample.

Requires Python 3.10+ and NumPy. From this directory run:
    python measure.py --download --cache raw --output measured.csv
Cache files are the unmodified public products. SHA-256 checks prevent silent
changes at the live archive from altering this campaign.
"""
import argparse
import csv
import hashlib
import io
import urllib.request
from pathlib import Path

import numpy as np

BASE = "https://www.swift.ac.uk"
HERE = Path(__file__).resolve().parent
BANDS = ((10000, 40000), (10000, 20000), (5000, 30000), (10000, 100000))
FIELDS = ["grb", "target_id", "bat_fluence_erg_cm2", "start_s", "stop_s",
          "bins", "exposure_s", "xrt_fluence_erg_cm2", "xrt_bat_ratio",
          "mean_flux_per_bat_fluence_s_inv"]


def public_url(target, kind):
    if kind == "rate":
        return f"{BASE}/xrt_curves/{target}/PCCURVE.qdp"
    return f"{BASE}/burst_analyser/{target}/xrt/xrt_flux_pc_XRTBAND.qdp"


def product(row, kind, cache, download):
    target = row["target_id"]
    file = cache / kind / f"{target}.qdp"
    if not file.exists() and download:
        file.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(public_url(target, kind), timeout=60) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != row[kind + "_sha256"]:
            raise ValueError(f"Archive has changed: {target} {kind}")
        file.write_bytes(data)
    data = file.read_bytes()
    if hashlib.sha256(data).hexdigest() != row[kind + "_sha256"]:
        raise ValueError(f"SHA-256 mismatch: {target} {kind}")
    columns = 15 if kind == "rate" else 6
    values = np.loadtxt(io.BytesIO(data), comments=("!", "READ", "NO"),
                        usecols=range(columns), ndmin=2)
    return values


def measure(row, rate, flux):
    # Match all three bin edges. Upper limits in separate files are excluded.
    rate_by_bin = {tuple(r[:3]): r for r in rate}
    flux_by_bin = {tuple(f[:3]): f for f in flux}
    if len(rate_by_bin) != len(rate) or len(flux_by_bin) != len(flux):
        raise ValueError(f"Duplicate bins in {row['target_id']}")
    intervals = []
    for key in sorted(rate_by_bin.keys() & flux_by_bin.keys()):
        t, plus, minus = key
        r, f = rate_by_bin[key], flux_by_bin[key]
        if min(r[3], r[12], f[3]) > 0:
            intervals.append((t + minus, t + plus, r[12], f[3]))
    if not intervals:
        return []
    points = np.array(intervals)
    out = []
    bat = float(row["bat_fluence_erg_cm2"])
    for start, stop in BANDS:
        # Complete bins only. Never bridge an orbit gap or assign its duration
        # as exposure. A spanning bin with GTIs uses its actual exposure.
        use = points[(points[:, 0] >= start) & (points[:, 1] <= stop)]
        if not len(use):
            continue
        exposure = sum(use[:, 2])
        fluence = use[:, 2] @ use[:, 3]
        out.append(dict(zip(FIELDS, (row["grb"], row["target_id"], bat, start,
            stop, len(use), exposure, fluence, fluence / bat, fluence / bat / exposure))))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=HERE / "raw")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "measured.csv")
    args = parser.parse_args()
    results = []
    with (HERE / "sample.csv").open() as source:
        for row in csv.DictReader(source):
            if row["status"] != "usable products":
                continue
            rate, flux = [product(row, k, args.cache, args.download) for k in ("rate", "flux")]
            results.extend(measure(row, rate, flux))
    results.sort(key=lambda r: (r["start_s"], r["stop_s"], -r["xrt_bat_ratio"], r["grb"]))
    with args.output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in results:
            writer.writerow({k: f"{v:.10g}" if isinstance(v, float) else v for k,v in row.items()})


if __name__ == "__main__":
    main()

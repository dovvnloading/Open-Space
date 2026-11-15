"""Compare a 17.30 km/s velocity feature in two overlapping 4.6-arcminute beams."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from measure import PAGES, ROOT, design, digitize, fit, serialize

SOURCES = {
    "12.20-01": ("blue", [(2, 648, 12, 167), (2, 650, 12, 167), (2, 648, 12, 167),
                           (0, 714, 14, 158), (0, 713, 20, 172), (0, 713, 20, 172)]),
    "12.18-01": ("red", [(2, 627, 10, 118), (2, 630, 10, 118), (2, 627, 10, 118),
                         (0, 714, 12, 120), (0, 713, 16, 153), (0, 713, 16, 153)]),
}


def analyze(cache, output):
    manifest = json.loads((ROOT / "g12-panels.json").read_text())
    results, curves, references = {}, {}, {}
    for source, (color, ticks) in SOURCES.items():
        panels, controls = [], []
        for i, page in enumerate(PAGES):
            name = source + "_" + page + ".jpg"
            path = cache / name
            if not path.exists():
                with urlopen(manifest[name]["url"], timeout=45) as response:
                    content = response.read()
                if hashlib.sha256(content).hexdigest() != manifest[name]["sha256"]:
                    raise ValueError(f"Live panel changed: {name}")
                path.write_bytes(content)
            if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[name]["sha256"]:
                raise ValueError(f"Unexpected panel contents: {path}")
            panels.append(digitize(path, i, color, ticks))
            # The ~19.8 km/s line appears as green in one beam and yellow in the other.
            controls.append(digitize(path, i, "green" if color == "blue" else "yellow", ticks))
        curves[source], references[source] = panels, controls
        training = np.vstack(panels[2:5])
        period, coeff = fit(training, 175, 215)
        phase_date = np.linspace(59000, 59000 + period, 10000)
        model = design(phase_date, period) @ coeff
        max_mjd = phase_date[np.argmax(model)]
        # Evaluate the frozen model including three panels outside the source-selection interval.
        phase_checks = []
        for data, control in zip(panels, controls):
            t, flux = data.T
            prediction = design(t, period) @ coeff
            delta, near = cKDTree(control[:, :1]).query(data[:, :1])
            matched = (delta < 2.5) & (control[near, 1] > 0.25) & (flux > 0.25)
            ratio = np.log(flux[matched] / control[near[matched], 1])
            phase_checks.append(dict(n=len(t), r=float(np.corrcoef(flux, prediction)[0, 1]),
                                     ratio_bins=len(ratio), log_ratio_r=float(np.corrcoef(ratio, prediction[matched])[0, 1]),
                                     flux_p10_p90=np.quantile(flux, [.1, .9]).tolist()))
        results[source] = dict(period_days=period, maximum_mjd=max_mjd, coefficients=coeff.tolist(),
                               panels=phase_checks,
                               leave_one_middle_panel_out=[fit(np.vstack([panels[j] for j in (2, 3, 4) if j != k]),
                                                                175, 215)[0] for k in (2, 3, 4)])
    with (output / "g12-bins.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["pointing", "feature", "panel", "mjd_approx", "flux_jy_approx"])
        for source in SOURCES:
            for label, collection in (("17.30", curves), ("19.8", references)):
                for page, data in zip(PAGES, collection[source]):
                    writer.writerows((source, label, page, f"{t:.3f}", f"{v:.4f}") for t, v in data)
    serialize(output / "g12-measurements.json", results)
    fig, axes = plt.subplots(2, 1, figsize=(9, 4.6), sharex=True, layout="constrained")
    for ax, (source, panels) in zip(axes, curves.items()):
        summary = results[source]
        for i, data in enumerate(panels):
            color = "#216b81" if i in (2, 3, 4) else "#b15c35"
            ax.scatter(data[:, 0], data[:, 1], c=color, s=7, zorder=3, alpha=.83)
        t = np.arange(56290, 61300, 2)
        ax.plot(t, design(t, summary["period_days"]) @ summary["coefficients"], color="#535872", lw=.9)
        ax.set(ylabel="Flux density (Jy)", title=f"iMet pointing {source}, 17.30 km/s")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Modified Julian Date")
    fig.savefig(output / "g12-cycle.png", dpi=135)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    analyze(args.cache, args.out)

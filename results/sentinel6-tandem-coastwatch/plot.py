"""Plot daily match geometry and the coastal comparison from summary.json."""

import argparse
from datetime import date
import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def plot(output):
    data = json.loads((ROOT / "summary.json").read_bytes())
    rows = [(date.fromisoformat(day), item["separation_km"])
            for day, item in data["daily_counts"].items() if item["coordinate_pairs"]]
    days, separations = zip(*rows, strict=True)
    median = np.array([row["median"] for row in separations])
    low = np.array([row["q25"] for row in separations])
    high = np.array([row["q75"] for row in separations])
    fig, (geometry, coastal) = plt.subplots(2, 1, figsize=(11, 7),
                                            layout="constrained", height_ratios=(1.35, 1))
    fig.set_facecolor("#f7f9fc")
    for axis in (geometry, coastal):
        axis.set_facecolor("#f7f9fc")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="x", alpha=0.2)

    geometry.axvspan(date(2026, 11, 1), date(2026, 11, 16), color="#81c9d1", alpha=0.24, label="November holdout")
    geometry.axhline(2, color="#9c531b", linestyle="--", label="Tighter 2 km cutoff")
    geometry.errorbar(days, median, yerr=[median - low, high - median], fmt="o",
                      markersize=3.5, capsize=0, color="#125a7a", alpha=0.85, label="Daily median and quartiles")
    geometry.xaxis.set_major_locator(mdates.WeekdayLocator(interval=2))
    geometry.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    geometry.set_title("Daily distance between matched Sentinel-6A and 6B samples", loc="left", fontweight="bold")
    geometry.set_ylabel("Surface distance (km)")
    geometry.legend(ncol=3, frameon=False, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, -0.24))
    geometry.set_ylim(0.65, 3.55)

    levels = [
        ("Sep 9–Oct 29, 5 km", "primary", "reference"),
        ("Nov 1–15, 5 km", "primary", "holdout"),
        ("Sep 9–Oct 29, 2 km", "two_km_sensitivity", "reference"),
        ("Nov 1–15, 2 km", "two_km_sensitivity", "holdout"),
    ]
    coastal.axvline(0, color="#6d7880", linewidth=1)
    coastal.axvline(10, color="#9c531b", linestyle="--", linewidth=1)
    coastal.set_yticks(range(len(levels)), [row[0] for row in levels])
    coastal.invert_yaxis()
    for row, (_, experiment, period) in enumerate(levels):
        values = data[experiment]["by_period"][period]
        stat = values["visit_near_minus_offshore_m"]
        value = 1000 * stat["median"]
        left, right = np.array(values["visit_day_bootstrap_95_m"]) * 1000
        color = "#125a7a" if experiment == "primary" else "#43868c"
        coastal.errorbar(value, row, xerr=[[value - left], [right - value]], marker="o",
                         capsize=3, color=color)
        coastal.text(12.1, row, f"{stat['n']} visits", ha="left", va="center", fontsize=9)
    coastal.set_xlim(-8, 17)
    coastal.set_title("Near-shore minus offshore mismatch by track visit", loc="left", fontweight="bold")
    coastal.set_xlabel("Difference in median absolute SLA disagreement (mm)")
    fig.savefig(output, dpi=150, facecolor=fig.get_facecolor())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    plot(parser.parse_args().output)

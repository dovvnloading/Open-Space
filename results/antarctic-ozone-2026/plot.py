"""Plot the location and size of the jointly observed November ozone core."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

from analyze import ROOT, daily_deep_masks, read


def run(output):
    days, complete = daily_deep_masks(read())
    daily_core = complete & np.all(days[2026] < 180, axis=(0, 1))
    categories = complete[:30].astype(int) + daily_core[:30].astype(int)
    summary = json.loads((ROOT / "summary.json").read_bytes())
    fig = plt.figure(figsize=(10.8, 5.5), layout="constrained")
    layout = fig.add_gridspec(1, 2, width_ratios=[1, 1.16])
    polar = fig.add_subplot(layout[0], projection="polar")
    bars = fig.add_subplot(layout[1])

    ring = np.arange(31)
    angle = np.deg2rad(np.arange(-180, 181))
    palette = ListedColormap(("#dce1e5", "#e7f4f6", "#086c84"))
    polar.pcolormesh(angle, ring, categories, cmap=palette, vmin=-.5, vmax=2.5, shading="flat")
    polar.set_theta_zero_location("N")
    polar.set_xticks(np.deg2rad(np.arange(0, 360, 60)),
                     ["0°", "60°E", "120°E", "180°", "120°W", "60°W"])
    polar.set_yticks([10, 20], labels=["80°S", "70°S"])
    polar.set_rlim(0, 30)
    polar.tick_params(labelsize=9)
    polar.set_title("2026 deep depletion on all 15 daily maps", pad=22, weight="bold", fontsize=11)
    polar.legend(handles=[
        Patch(color="#086c84", label="Both instruments below 180 DU on all 15 dates"),
        Patch(color="#e7f4f6", label="Other cells with complete coverage"),
        Patch(color="#dce1e5", label="Incomplete coverage across the three years"),
    ], loc="upper center", bbox_to_anchor=(.5, -.14), frameon=False, fontsize=8)

    years = ["2024", "2025", "2026"]
    measures = [summary["seasons"][year]["fixed_across_years"]["180"] for year in years]
    rows = np.arange(len(years))
    initial = [item["early_area_mkm2"] for item in measures]
    retained = [item["retained_area_mkm2"] for item in measures]
    bars.barh(rows, initial, color="#b7d7df", height=.43, label="Deep on at least four early dates")
    bars.barh(rows, retained, color="#086c84", height=.43, label="Also deep on all five late dates")
    for position, amount in zip(rows, retained):
        bars.text(amount + .07, position + .02, f"{amount:.2f}", va="center", color="#14394a", fontsize=9)
    bars.set_yticks(rows, years)
    bars.invert_yaxis()
    bars.set(title="Early-November area still depleted on November 11–15",
             xlabel="Area (million km²) on the fixed two-window footprint")
    bars.title.set_fontweight("bold")
    bars.title.set_fontsize(11)
    bars.spines[["top", "right", "left"]].set_visible(False)
    bars.grid(axis="x", color="#e5edef")
    bars.set_axisbelow(True)
    bars.legend(loc="upper right", frameon=False, fontsize=8)
    fig.savefig(output, dpi=170, metadata={"Software": "Open-Space ozone analysis"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output)

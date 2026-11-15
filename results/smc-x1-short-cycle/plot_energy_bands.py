"""Plot measured band contrasts and the separate XRT pointings in the four Swift cycles."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
BLUE, RUST = "#245b80", "#a9472b"


def main():
    maxi = json.loads((ROOT / "energy-bands.json").read_text())
    xrt = json.loads((ROOT / "xrt-counts.json").read_text())
    fig, panels = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    ax, bx = panels
    x = np.arange(3)
    for label, key, color, shift, marker in [
        ("Outside eclipse", "out_of_eclipse", BLUE, -.08, "o"),
        ("Eclipse core", "eclipse_core", RUST, .08, "s"),
    ]:
        y = np.array(maxi[key]["middle_minus_edge"])
        lower, upper = np.array(maxi[key]["day_bootstrap_95_contrast"]).T
        ax.errorbar(x + shift, y, yerr=[y - lower, upper - y], color=color, fmt=marker,
                    capsize=3, markersize=5, label=label)
    ax.axhline(0, color="#777777", linewidth=.8)
    ax.set_xticks(x, ["2–4", "4–10", "10–20"])
    ax.set_ylabel("Middle minus outer flux (photons cm⁻² s⁻¹)")
    ax.set_title("MAXI / Gas Slit Camera", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=9)

    for state, color, shift, marker, label in [
        ("middle", BLUE, -.10, "o", "Middle third (4 visits)"),
        ("edge", RUST, .10, "s", "Outer sixths (2 visits)"),
    ]:
        subset = xrt["outside_eclipse"][state]
        visits = [v for v in xrt["visits"] if v["obsid"] in subset["obsids"]]
        for i, visit in enumerate(visits):
            s = np.array(visit["counts"]["wing6"])
            b = np.array(visit["counts"]["background"])
            area = visit["exposure_area_px_s"]
            exposure = area["wing6"] / (np.pi * (20**2 - 6**2))
            net = (s - b * area["wing6"] / area["background"]) / exposure
            bx.scatter(x + shift + (i - (len(visits)-1)/2)*.03, net, marker=marker, s=20,
                       facecolors="none", edgecolors=color, linewidths=.75, alpha=.65)
        bx.scatter(x + shift, subset["net_annular_rates_per_s"]["wing6"],
                   marker=marker, s=42, color=color, label=label)
    bx.set_xticks(x, ["0.3–1", "1–2", "2–10"])
    bx.set_ylabel("Net annular rate (counts s⁻¹)")
    bx.set_title("Swift / X-Ray Telescope", loc="left", fontweight="bold")
    bx.legend(frameon=False, fontsize=9)
    for axis in panels:
        axis.set(xlabel="Energy band (keV)", xlim=(-.45, 2.45))
        axis.spines[["right", "top"]].set_visible(False)
        axis.grid(axis="y", color="#dedede", linewidth=.6)
        axis.set_axisbelow(True)
    fig.savefig(ROOT / "xray-energy-bands.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    main()

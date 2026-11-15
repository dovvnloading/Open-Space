"""Plot the eight unconstrained TESS centroids on a common sky-direction scale."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def main():
    data = json.loads((ROOT / "measurements.json").read_text())
    sectors = list(data)
    fig, (ax, flux) = plt.subplots(1, 2, figsize=(9.2, 4.7), sharey=True,
                                  gridspec_kw={"width_ratios": [1.5, 1]})
    for i, sector in enumerate(sectors):
        row = data[sector]
        baseline = row["toi_to_binary_px"]
        center = row["projected_offset_px"] / baseline
        low, high = np.asarray(row["projected_bootstrap_5_95_px"]) / baseline
        color = "#246486" if sector not in ("12", "65") else "#996041"
        ax.plot([low, high], [i, i], color=color, lw=2, solid_capstyle="round")
        ax.scatter(center, i, color=color, s=26, zorder=3)
        f = row["phase_aperture_electrons_per_second"]
        flux.scatter(f["0"], i, color=color, s=28, zorder=2)
        flux.scatter([f["-0.25"], f["0.25"]], [i-.10, i+.10], marker="|", color="#808995", s=56)
    ax.axvline(0, color="#888d91", ls="--", lw=1)
    ax.axvline(1, color="#246486", ls=":", lw=1.3)
    ax.set_xlabel("Centroid projected along target → binary\n(target = 0, neighboring binary = 1)")
    ax.set_yticks(range(len(sectors)), [f"Sector {s}" for s in sectors])
    ax.invert_yaxis()
    flux.axvline(0, color="#888d91", lw=.8)
    flux.set_xscale("symlog", linthresh=3)
    flux.set_xlim(-30, 260)
    flux.set_xticks([-20, 0, 3, 15, 180], ["−20", "0", "3", "15", "180"])
    flux.set_xlabel("Aperture flux loss (electrons/s)\npoint = eclipse, ticks = controls")
    ax.set_title("Where the light disappears", loc="left", fontsize=12)
    flux.set_title("Eclipse and off-phase controls", loc="left", fontsize=12)
    for a in (ax, flux):
        a.spines[["right", "top"]].set_visible(False)
        a.grid(axis="y", alpha=.15)
        a.tick_params(length=0, pad=5)
    fig.suptitle("TOI-3143.01 across eight TESS sectors", x=.12, ha="left", fontsize=15, fontweight="bold")
    fig.subplots_adjust(wspace=.17, left=.12, right=.97, top=.85, bottom=.19)
    fig.savefig(ROOT / "centroids.png", dpi=165, facecolor="white")


if __name__ == "__main__":
    main()

"""Plot resolved ESPRESSO calcium/hydrogen profiles and the epoch-wise activity test."""
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from replay import EPOCH, ROOT, design, remove

BLUE = "#195d78"
ORANGE = "#b55e38"
GRAY = "#68717a"


def main():
    with (ROOT / "s1d-profiles.csv").open() as source:
        profiles = list(csv.DictReader(source))
    with (ROOT / "s1d-activity.csv").open() as source:
        activity = list(csv.DictReader(source))
    fig = plt.figure(figsize=(11, 6.1), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.2])
    labels = [("CaII_K", "Ca II K", 3933.664), ("CaII_H", "Ca II H", 3968.470),
              ("Halpha", "Hα", 6562.808)]
    for i, (key, title, center) in enumerate(labels):
        ax = fig.add_subplot(grid[0, i])
        selected = [r for r in profiles if r["feature"] == key]
        x = np.array([float(r["rest_air_angstrom"]) - center for r in selected])
        for visit in activity:
            ax.plot(x, [float(r[visit["product"]]) for r in selected],
                    color=BLUE, alpha=0.12, linewidth=0.55)
        ax.plot(x, [float(r["median"]) for r in selected], color=BLUE, linewidth=1.5,
                label="Median of 25 visits")
        ax.set(xlabel="Rest-frame offset (Å)", ylabel="Flux / local continuum" if i == 0 else "",
               title=f"{title}  ·  {center:.3f} Å", xlim=(-2, 2))
        if i == 0:
            ax.legend(loc="best", frameon=False, fontsize=8)

    ax = fig.add_subplot(grid[1, :2])
    t = np.array([float(r["bjd"]) for r in activity]) - EPOCH
    base = np.column_stack([np.ones(len(t)), t - t.mean()])
    models = json.loads((ROOT / "harps_predictions.json").read_text())
    frozen = min(models, key=lambda r: r["rms"])
    full = remove(design(t, frozen["P"]) @ frozen["params"][:8], base)
    ha = remove(np.array([float(r["halpha16"]) for r in activity]), base) * 1000
    years = np.floor((t + EPOCH - 2460175) / 365.25).astype(int)
    for season, color, symbol in zip(np.unique(years), [BLUE, ORANGE, GRAY], ["o", "s", "^"]):
        mask = years == season
        ax.scatter(full[mask], ha[mask], c=color, marker=symbol, s=30,
                   label=f"Season {season + 1} ({mask.sum()})")
    ax.axhline(0, color="#c8ced2", lw=0.7, zorder=0)
    ax.axvline(0, color="#c8ced2", lw=0.7, zorder=0)
    lo, hi = np.min(full), np.max(full)
    slope = np.dot(full, ha) / np.dot(full, full)
    ax.plot([lo, hi], [lo * slope, hi * slope], c="#333f46", lw=1.2)
    ax.set(xlabel="Frozen HARPS waveform, detrended (m/s)",
           ylabel="Hα index residual (× 10⁻³)", title="Hα and the four-component velocity prediction")
    ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)

    ax = fig.add_subplot(grid[1, 2])
    with (ROOT / "s1d-velocities.csv").open() as source:
        velocity = list(csv.DictReader(source))
    for i, band in enumerate(["4200_4600", "4600_5000", "5000_5500"]):
        rvs = remove(np.array([float(r["rv_" + band + "_m_s"]) for r in velocity]), base)
        for val, marker, color, label in [(rvs, "o", GRAY, "Offset + trend"),
                                           (rvs - full, "s", BLUE, "Plus frozen waveform")]:
            ax.scatter(i, np.sqrt(np.mean(val**2)), marker=marker, color=color,
                       label=label if i == 0 else None, s=42)
    ax.set(xticks=range(3), xticklabels=["420–460", "460–500", "500–550"],
           xlabel="Vacuum wavelength (nm)", ylabel="RMS velocity (m/s)", ylim=(0, 5),
           title="Flux-derived velocities")
    ax.legend(frameon=False, fontsize=8, loc="lower center")
    for axis in fig.axes:
        axis.spines[["right", "top"]].set_visible(False)
        axis.tick_params(labelsize=8)
    fig.suptitle("HD 32564   |   ESPRESSO spectroscopy", fontsize=16, weight="bold")
    fig.savefig(ROOT / "s1d-evidence.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

"""Plot the measured conditional MIRI injection limits for GJ 205."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt

d = Path(__file__).resolve().parent
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8), layout="constrained")
for band in ["f1065c", "f1550c"]:
    rows = list(csv.DictReader((d / (band + "_miri_contrast.csv")).open()))
    sep = [float(r["radius_arcsec"]) for r in rows]
    axes[0].plot(
        sep, [float(r["f95_uJy"]) for r in rows], label=band.upper(), marker="o"
    )
    axes[1].plot(
        sep,
        [float(r["contrast95"]) * 1e5 for r in rows],
        label=band.upper(),
        marker="o",
    )
for ax in axes:
    ax.set_xlabel("Separation (arcsec)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.23)
axes[0].set_ylabel("95% sampled recovery flux (µJy)")
axes[1].set_ylabel("95% sampled recovery contrast (×10⁻⁵)")
fig.savefig(d / "miri_contrast.png", dpi=150)

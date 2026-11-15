"""Render the Kepler-18 orientation comparison. Run from this directory."""

from pathlib import Path
import pandas as pd
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

root = Path(__file__).resolve().parent
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.right": False,
        "axes.spines.top": False,
    }
)
fig, (upper, lower) = plt.subplots(
    2, 1, figsize=(8.1, 5.5), height_ratios=[1.1, 1], layout="constrained"
)
angles = pd.read_csv(root / "orientations.csv")
models = [
    ("DR25", "#9b4d23", "DR25"),
    ("Cochran2011", "#5e6475", "Cochran et al. (2011)"),
    ("ShortCadence", "#125b69", "175 isolated Kepler transits"),
]
for key, color, label in models:
    rows = angles[
        (angles.geometry == key)
        & (angles.photometry == "pdc")
        & (angles.cadence == "best")
        & (angles.sign_d == 1)
    ]
    upper.plot(rows.angle_deg, rows.observed_delta_chi2, color=color, label=label, lw=1.9)
upper.axhline(0, color="black", lw=0.6)
upper.set(
    xlim=(-180, 180),
    ylabel="Linearized Δχ²",
    xlabel="Relative sky-plane rotation (degrees)",
    title="Kepler-18 c/d  |  Same-side chord configurations",
)
upper.legend(loc="upper left", frameon=False, fontsize=8)
direct = pd.read_csv(root / "nonlinear.csv").groupby(["col", "ang"], as_index=False).delta.sum()
for key, color, marker, label in [
    ("pdc", "#125b69", "o", "PDCSAP"),
    ("sap", "#ca713c", "s", "SAP"),
]:
    rows = direct[direct.col == key]
    lower.plot(rows.ang, rows.delta, color=color, marker=marker, label=label, lw=1.6, ms=5)
lower.set(
    xlim=(57.5, 64.5),
    xticks=[58, 60, 62, 64],
    ylim=(0, 30),
    ylabel="Refitted Δχ²",
    xlabel="Relative sky-plane rotation (degrees)",
    title="Independent shape fit  |  Direct nonlinear checks",
)
lower.legend(loc="lower left", frameon=False, ncols=2)
fig.savefig(root / "orientation-test.png", dpi=145)

"""Plot the archived astrometry and the backward-integration sensitivity samples."""

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
M = ROOT / "measurements"
if __name__ == "__main__":
    history = pd.read_csv(M / "history.csv")
    samples = pd.read_csv(M / "clones.csv").iloc[1:]
    fig, axes = plt.subplots(1, 2, figsize=(10.7, 4.15), constrained_layout=True)
    a, b = axes
    a.plot(
        history.lookback_years / 1000,
        history.nominal_separation_km,
        lw=0.9,
        color="#20506a",
    )
    a.set_yscale("log")
    a.set(
        xlabel="Thousands of years before June 2026",
        ylabel="Separation (km)",
        title="Nominal orbit history",
    )
    points = b.scatter(
        samples.lookback_years / 1000,
        samples.sampled_min_km,
        c=samples.relative_m_s,
        cmap="viridis",
        s=32,
        edgecolors="white",
        linewidths=0.35,
    )
    b.set(
        xlabel="Thousands of years before June 2026",
        ylabel="Smallest sampled separation (km)",
        title="48 perturbed orbit pairs",
    )
    for axis in axes:
        axis.grid(alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)
    fig.colorbar(points, ax=b, label="Relative speed (m/s)")
    fig.savefig(ROOT / "encounter.png", dpi=150)

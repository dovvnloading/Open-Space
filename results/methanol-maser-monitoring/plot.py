"""Plot digitized, measured bins separately from the fitted template."""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from measure import design

ROOT = Path(__file__).resolve().parent
model = json.loads((ROOT / "measurements.json").read_text())["blue"]
period, maximum = model["period_days"], model["maximum_mjd"]
rows = [r for r in csv.DictReader((ROOT / "digitized.csv").open()) if r["color"] == "blue"]
time = np.array([float(r["mjd_approx"]) for r in rows])
flux = np.array([float(r["flux_jy_approx"]) for r in rows])
train = np.array([r["panel"] in ("3rd", "4th", "5th") for r in rows])
torun = np.genfromtxt(ROOT / "torun/torun.csv", delimiter=",", names=True)
torun = torun[torun["velocity_km_s"] == 16.92]
phase = lambda t: (t - maximum) / period % 1.0

fig, axes = plt.subplots(1, 2, figsize=(9, 3.65), layout="constrained")
grid = np.linspace(0, 1, 300)
ax = axes[0]
ax.scatter(phase(time[train]), flux[train], c="#90a4b3", alpha=.6, s=7, label="Model fitting panels")
ax.scatter(phase(time[~train]), flux[~train], facecolors="none", edgecolors="#137a77",
           linewidths=.7, s=15, label="Earlier and later panels")
ax.plot(grid, design(maximum + grid * period, period) @ model["coefficients"], color="#243a5a", lw=2)
ax.set(title="Ibaraki · ~17.01 km/s", ylabel="Digitized flux density (Jy)")
ax.legend(frameon=False, fontsize=8, loc="upper center")

ax = axes[1]
t, f = torun["mjd"], torun["flux_jy"]
trend = t - np.mean(t)
matrix = np.column_stack((np.ones_like(t), trend, np.cos(2 * np.pi * phase(t))))
mean, slope, amplitude = np.linalg.lstsq(matrix, f, rcond=None)[0]
ax.scatter(phase(t), f - slope * trend, color="#b1572e", s=21, alpha=.85, label="44 distinct epochs")
ax.plot(grid, mean + amplitude * np.cos(2 * np.pi * grid), c="#6b361d", lw=2)
ax.set(title="Toruń · 16.92 km/s", ylabel="Trend-corrected flux density (Jy)")
ax.legend(frameon=False, fontsize=8)
for ax in axes:
    ax.set(xlabel="Phase of the Ibaraki ephemeris", xlim=(-.025, 1.025))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.16)
fig.savefig(ROOT / "phase.png", dpi=145)

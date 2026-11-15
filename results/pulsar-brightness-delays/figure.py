"""Plot the measured slopes and their separate receiver-band test."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).parent
with (root / "visits.csv").open() as file:
    rows = [r for r in csv.DictReader(file) if r["pulsar"] == "J1825-0935"]


def col(k):
    return np.array([float(r[k]) for r in rows])


f = col("flux_mjy")
t = (col("mjd") - col("mjd").mean()) / 365.25
x = np.column_stack((np.ones(len(t)), t, f))
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
for panel, (ax, y, title) in enumerate(zip(axes, [col("dm"), col("relative_dm_4_7") - col("relative_dm_0_3")],
                                         ["Published dispersion measure", "High-band minus low-band estimate"])):
    offset, drift, slope = np.linalg.lstsq(x, y, rcond=None)[0]
    corrected = y - drift * t
    corrected -= np.median(corrected)
    ax.scatter(f, corrected, s=19, color="#006d77", alpha=.8, edgecolors="none")
    fit = offset + slope * f - np.median(y - drift * t)
    ax.plot(f[np.argsort(f)], fit[np.argsort(f)], color="#983d26", lw=1.7)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Flux density (mJy)")
    ax.set_ylabel("Relative DM (pc cm⁻³)")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("PSR B1822−09 · 52 MeerKAT visits", fontweight="bold")
fig.savefig(root / "flux-delay.png", dpi=155)

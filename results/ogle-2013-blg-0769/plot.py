"""Plot the archived OGLE and KMTNet nightly measurements."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

here = Path(__file__).resolve().parent
data = pd.read_csv(here / "nightly.csv")
fig, axes = plt.subplots(1, 3, figsize=(12.8, 3.6), constrained_layout=True)
for ax, year, center in zip(axes[:2], [2013, 2016], [56431, 57536]):
    group = data[data.site_or_extraction == f"OGLE_{year}"]
    ax.scatter(group.hjd_minus_2400000_floor - center, group["median"], s=19, color="#176a85")
    ax.axhline(group.baseline.iloc[0], color="#758487", lw=1, ls="--")
    ax.set(xlim=(-18,18), ylim=(20.9,17.45), xlabel=f"Days from HJD {2400000+center}",
           ylabel="OGLE I magnitude", title=f"{year}   OGLE")
ax = axes[2]
for site, name, color in [("BLG","Chile","#248596"),("SAO","South Africa","#b65424"),("SSO","Australia","#315692")]:
    obs = data[(data.site_or_extraction == f"scene_B_{site}") & (data.center == 60525)]
    ax.scatter(obs.hjd_minus_2400000_floor-60525, obs.contrast, color=color, s=15, label=name)
ax.axhline(0, color="#758487", lw=1, ls="--")
ax.set(xlim=(-18,18), xlabel="Days from HJD 2460525", ylabel="Difference flux / baseline scatter",
       title="2024   KMTNet")
ax.legend(fontsize=8, frameon=False)
for ax in axes:
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", alpha=.13)
fig.savefig(here / "recurrence.png", dpi=155)

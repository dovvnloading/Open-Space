"""Visit timing residuals against each binary clock measured in sectors 5 and 32."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
events = pd.read_csv(ROOT / "timings.csv.gz")
visits = pd.read_csv(ROOT / "visits.csv")
events = events[(events.aperture == "box") & (events.half_window_d == .13)]
visits = visits[(visits.aperture == "box") & (visits.half_window_d == .13)]
colors = ["#075d89", "#097a54", "#b65319"]
fig, axes = plt.subplots(1, 3, figsize=(9, 3.5), sharey=True)
for ax, (label, group), color in zip(axes, visits.groupby("clock"), colors):
    group = group.sort_values("sector")
    t = group.median_btjd.to_numpy()
    delta = group.residual_min.to_numpy()
    slope = (delta[1] - delta[0]) / (t[1] - t[0])
    project = lambda x: delta[0] + slope * (x - t[0])
    for ix, row in enumerate(group.itertuples()):
        e = events[(events.clock == label) & (events.sector == row.sector)]
        residual = e.residual_min - project(e.nominal_btjd)
        jitter = (e.nominal_btjd - e.nominal_btjd.median()) * .005
        ax.scatter(ix+jitter, residual, s=8, c=color, alpha=.25, edgecolors="none")
        ax.plot(ix, row.residual_min - project(row.median_btjd), "D", color=color, ms=6)
    ax.axhline(0, color="#888", lw=.7)
    ax.set_xticks(range(4), ["5", "32", "98", "109"])
    ax.set_xlim(-.45,3.45)
    ax.set_xlabel("TESS sector")
    ax.set_title(label, color=color, fontsize=12)
    ax.spines[["right","top"]].set_visible(False)
axes[0].set_ylabel("Eclipse arrival offset (min)")
fig.suptitle("TIC 43472931  ·  eclipse times relative to the 2018–2020 clocks", fontsize=13)
fig.tight_layout()
fig.savefig(ROOT / "timing.png", dpi=170)

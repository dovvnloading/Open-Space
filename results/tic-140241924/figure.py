"""Draw independent-ephemeris residuals for the two eclipsing binaries."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
seasons = pd.read_csv(ROOT / "seasons.csv").query('aperture == "box" and window == 0.17')
visits = pd.read_csv(ROOT / "visits.csv").query('aperture == "box" and window == 0.17')
curves = pd.read_csv(ROOT / "curvature.csv").query('aperture == "box" and window == 0.17')
fig, ax = plt.subplots(figsize=(7.7, 4.3))
styles = {"A primary": ("#06598b", "o"), "A secondary": ("#0e835f", "s"), "B primary": ("#ac4627", "D")}
for clock, g in seasons.groupby("clock"):
    t = g.mean_btjd.to_numpy()
    y = g.delta.to_numpy()
    slope = (y[-1] - y[0]) / (t[-1] - t[0])
    ref = lambda x: y[0] + slope*(x - t[0])
    col, mark = styles[clock]
    v = visits[visits.clock == clock]
    ax.scatter(v.median_btjd, v.dt - ref(v.median_btjd), s=14, color=col, alpha=.4)
    residual = y - ref(t)
    ax.plot(t, residual, color=col, alpha=.55, lw=1)
    offset = {"A primary": -14, "A secondary": 14, "B primary": 0}[clock]
    ax.scatter(t+offset, residual, marker=mark, s=50, color=col, label=clock)
    ci = curves[(curves.clock == clock) & curves.metric.str.contains("chord")]
    ax.errorbar(t[1:3]+offset, residual[1:3], yerr=np.array([ci.value-ci.low,ci.high-ci.value]), color=col, fmt="none", capsize=3)
ax.axhline(0, color="#777", alpha=.5, lw=.8)
ax.set_xticks([1556,2268,3002,3718], ["2019", "2021", "2023", "2025"])
ax.set_ylabel("Eclipse arrival offset (minutes)")
ax.set_xlabel("TESS observing season")
ax.set_title("TIC 140241924 · offsets from separate 2019–2025 linear clocks")
ax.legend(ncol=3, loc="lower right", fontsize=9)
ax.spines[["right","top"]].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT/"timing.png",dpi=180)

"""Make the research figure from the released numerical tables."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
read = lambda name: np.genfromtxt(root / name, delimiter=",", names=True)
moa = read("moa-nights.csv")
episodes = read("moa-episodes.csv")
day = read("daily.csv")
pix = read("k2-pixels.csv")
fig = plt.figure(figsize=(9.2, 6.5), layout="constrained")
g = fig.add_gridspec(2, 2, width_ratios=(1.7, 1), height_ratios=(.88, 1))
a, b, c = fig.add_subplot(g[0, :]), fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1])
a.scatter(moa["HJD_minus_2450000"], moa["corrected_flux"]/1000, s=5,
          color="#80939e", alpha=.8, linewidths=0, rasterized=True)
for ix, (peak, label) in enumerate(zip(episodes, ["2007", "2008", "2010", "2013", "2014"])):
    t, y = peak[0], peak[1]/1000
    a.scatter(t, y, color="#a64127", s=22, zorder=3)
    a.annotate(label, (t, y), xytext=(0, 6), textcoords="offset points", ha="center", size=9)
a.set(title="Five bright episodes in the Mount John difference photometry",
      xlabel="Heliocentric Julian date − 2,450,000", ylabel="MOA difference flux\n(thousands of archive units)",
      ylim=(-2.0, 10.0))
sel = day["integer_JD_minus_2450000"] <= 7547
q = day[sel]
t = q["median_K2_BJD_minus_2450000"]
k2 = q["K2_weighted_residual_e_per_s"]
norm = np.median(k2[(t>7533)&(t<7539)])
b.plot(t, k2 / norm, "o-", lw=1.6, ms=4, color="#146161", label="Kepler")
b.plot(q["integer_JD_minus_2450000"]+.5, q["KMT_normalized_ground_flux"], "s--",
       lw=1.2, ms=4, color="#a64127", label="KMTNet")
b.axhline(0, color="#b7bec0", lw=.8)
b.set(title="The 2016 eruption from space and ground", xlabel="Julian date − 2,450,000",
      ylabel="Normalized excess flux")
b.legend(frameon=False)
field = np.full((10, 14), np.nan)
for x, y, value in pix:
    if 547 <= x <= 560 and 762 <= y <= 771:
        field[int(y-762),int(x-547)] = value
im = c.imshow(field, origin="lower", extent=(546.5, 560.5, 761.5, 771.5),
              cmap="cividis", vmin=-15, vmax=125, interpolation="nearest")
c.scatter(554.087, 763.630, marker="+", color="#fa441d", s=105, lw=1.5)
c.set(title="Kepler difference image", xlabel="Detector column", ylabel="Detector row")
fig.colorbar(im, ax=c, location="bottom", shrink=.84, label="Residual flux (electrons s⁻¹)")
for ax in (a,b):
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", alpha=.15)
fig.savefig(root/"observations.png", dpi=140)

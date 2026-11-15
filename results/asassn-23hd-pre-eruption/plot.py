"""Draw the image-space cross-identification and the archival light curve."""
import json
from datetime import datetime, timedelta

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from analyze import INPUTS, ROOT, read_image

report = json.loads((ROOT / "measurements.json").read_text())
fig = plt.figure(figsize=(9.0, 6.0), layout="constrained")
grid = fig.add_gridspec(2, 3, height_ratios=[1.05, 1])
panels = [("reference.fits", "quiescent_reference", "2012 reference"),
          ("science.fits", "science", "2014 exposure"),
          ("difference.fits", "difference", "2014 difference")]
for col, (filename, key, title) in enumerate(panels):
    image, date = read_image(filename)
    fit = report["images"][key]
    x, y = fit["gaia_pixel_xy"]
    scale = fit["scale_arcsec_per_pixel"]
    dx, dy = fit["fits"]["7"]["pixel_xy"]
    x0, y0 = int(x) - 19, int(y) - 19
    cut = image[y0:y0+39, x0:x0+39]
    edges = [(x0 - x - .5) * scale, (x0 + 38.5 - x) * scale,
             (y0 - y - .5) * scale, (y0 + 38.5 - y) * scale]
    ax = fig.add_subplot(grid[0, col])
    ax.imshow(cut, origin="lower", cmap="gray_r", extent=edges,
              vmin=np.percentile(cut, 16), vmax=np.percentile(cut, 99.8))
    ax.scatter(0, 0, marker="+", s=180, color="#bd3e26", linewidths=1.5,
               label="Gaia counterpart", zorder=3)
    ax.scatter((dx-x)*scale, (dy-y)*scale, facecolor="none", edgecolor="#176082", s=62,
               linewidths=1.45, label="Fitted source")
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_xticks([-4, 0, 4])
    ax.set_yticks([-4, 0, 4])
    ax.set_xlabel("Image x offset (arcsec)", fontsize=8.5)
    if col == 0:
        ax.set_ylabel("Image y offset (arcsec)", fontsize=8.5)
    ax.tick_params(labelsize=8)
    if col == 2:
        ax.legend(loc="upper left", fontsize=7, framealpha=.93, handlelength=1)
ax = fig.add_subplot(grid[1, :])
data = np.loadtxt(INPUTS / "ogle_photometry.dat")
dates = np.array([datetime(2000,1,1,12) + timedelta(days=jd-2451545) for jd in data[:,0]])
selected = data[:,2] > 0
ax.errorbar(dates[selected], data[selected,1], yerr=data[selected,2], fmt="o",
            color="#235776", ecolor="#708a98", markersize=4, elinewidth=.9, zorder=2,
            label="OGLE I measurement")
ax.plot(dates[~selected], np.full((~selected).sum(), .045), "|", color="#9da4a8",
        markersize=5, transform=ax.get_xaxis_transform(), label="Visit with flagged photometry")
ax.set_ylim(22.05,19.15)
ax.set_ylabel("OGLE I (mag)")
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax.tick_params(labelsize=8.5)
ax.grid(axis="y", alpha=.18)
ax.legend(loc="upper right", fontsize=8, frameon=False)
ax.set_title("Pre-eruption photometry", loc="left", fontsize=11)
fig.suptitle("ASASSN-23hd progenitor field", x=.01, ha="left", fontsize=16)
fig.savefig(ROOT / "field-and-photometry.png", dpi=170)

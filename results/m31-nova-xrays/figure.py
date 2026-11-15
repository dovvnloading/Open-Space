"""Plot archived HRC events at both optical positions using a common count-rate scale."""
import csv
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
import numpy as np
from scipy.ndimage import gaussian_filter

root = Path(__file__).resolve().parent
with gzip.open(root / "photons.csv.gz", "rt") as stream:
    photons = list(csv.DictReader(stream))
with (root / "observations.csv").open() as stream:
    observations = list(csv.DictReader(stream))
eras = [("October 2021", {"2021"}), ("2023–2024", {"2023", "2024"})]
targets = ["M31N 2017-06c", "M31N 2021-06b"]
fig, axes = plt.subplots(2, 2, figsize=(6.4, 6.2), layout="constrained")
edges = np.linspace(-5, 5, 41)
for col, (label, years) in enumerate(eras):
    obs = [r for r in observations if r["date_obs"][:4] in years and "/hrcf" in r["url"]]
    obsids = {r["obsid"] for r in obs}
    exposure_ks = sum(float(r["exposure_s"]) for r in obs) / 1000
    for row, name in enumerate(targets):
        axis = axes[row, col]
        selected = [p for p in photons if p["obsid"] in obsids and p["target"] == name]
        x = [float(p["x_offset_arcsec"]) for p in selected]
        y = [float(p["y_offset_arcsec"]) for p in selected]
        density = np.histogram2d(x, y, bins=(edges, edges))[0].T / (exposure_ks * 0.25**2)
        img = axis.imshow(gaussian_filter(density, 2), origin="lower", extent=(-5,5,-5,5),
                          norm=Normalize(0,0.42), cmap="magma", interpolation="nearest")
        axis.add_patch(Circle((0,0), 2, fill=False, lw=1.0, ec="white", ls=":"))
        axis.set_xticks([-4,-2,0,2,4])
        axis.set_yticks([-4,-2,0,2,4])
        if col == 0:
            axis.set_ylabel(f"{name}\nSky y offset (arcsec)")
        if row == 0:
            axis.set_title(f"{label}\n{exposure_ks:.1f} ks HRC")
        if row == 1:
            axis.set_xlabel("Sky x offset (arcsec)")
fig.colorbar(img, ax=axes, shrink=.75, label="Counts per ks per square arcsec")
fig.savefig(root / "hrc-comparison.png", dpi=150)

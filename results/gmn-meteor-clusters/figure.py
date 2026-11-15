"""Plot the four preferred GMN tracks and their geocentric radiants."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
group = pd.read_csv(root/"geographic/groups.csv").query("multiplicity == 4").iloc[0, 1:]
d = pd.read_csv(root/"tracks.csv").set_index("id").loc[group]
t = pd.to_datetime(d.utc)
colors = ["#176b87", "#a44832", "#5c5393", "#317553"]
fig, (site, sky) = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
for k, (_, row) in enumerate(d.iterrows()):
    tag = "ABCD"[k]
    site.plot([row.lon, row.lon2], [row.lat, row.lat2], color=colors[k], lw=2)
    site.scatter(row.lon, row.lat, color=colors[k], s=40)
    site.annotate(f"{tag}   +{(t.iloc[k]-t.iloc[0]).total_seconds():.2f} s",
                  (row.lon, row.lat), xytext=(6, 7), textcoords="offset points", fontsize=9)
    sky.errorbar(row.ra, row.dec, xerr=row.sigma_ra, yerr=row.sigma_dec, fmt="o",
                 color=colors[k], capsize=2, alpha=.85)
    sky.annotate(tag, (row.ra, row.dec), xytext=(6, 5), textcoords="offset points", fontsize=9)
site.set(xlabel="Longitude (degrees east)", ylabel="Latitude (degrees north)",
         title="Atmospheric entry positions", xlim=(6.8, 15.6), ylim=(51.8, 57.1))
site.set_aspect(1/np.cos(np.deg2rad(54)))
sky.set(xlabel="Right ascension (degrees)", ylabel="Declination (degrees)",
        title="Apparent origin after correcting for Earth’s gravity")
for ax in (site, sky):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=.2)
fig.suptitle("Four meteor trajectories  •  24 August 2026, 00:36 UTC", fontsize=15)
fig.savefig(root/"tracks.png", dpi=150)

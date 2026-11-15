"""Plot time-resolved soft events and the summer Swift PI distribution."""
from datetime import date
from pathlib import Path

import matplotlib.dates as dates
import matplotlib.pyplot as plt
import numpy as np

from xray_bands import read

here = Path(__file__).resolve().parent
rows = [r for r in read("swift-timeline-counts.csv") if r["center"] == "nova" and r["radius_pix"] == "6"]
fig, (lc, sp) = plt.subplots(1, 2, figsize=(10, 3.9), width_ratios=(1.55, 1), layout="constrained")
visits = sorted((r for r in rows if (r["pi_min"], r["pi_max_inclusive"]) == ("30", "99")),
                key=lambda r: r["date_utc"])
t = [date.fromisoformat(r["date_utc"]) for r in visits]
on = np.array([int(r["on"]) for r in visits])
bg = np.array([float(r["alpha"]) * int(r["off"]) for r in visits])
lc.scatter(t, on, color="#17577f", s=19, label="Aperture")
lc.scatter(t, bg, color="#c45c28", marker="+", s=21, label="Scaled annulus")
lc.axvline(date(2021, 6, 12), color="#555555", ls=":", label="Optical eruption")
lc.set(ylabel="0.3–1 keV counts per visit", title="Swift/XRT · M31N 2021-06b sightline",
       xlabel="UTC date (visit exposures vary)")
lc.xaxis.set_major_locator(dates.MonthLocator(interval=3))
lc.xaxis.set_major_formatter(dates.DateFormatter("%b %Y"))
lc.legend(frameon=False, fontsize=8)
bands = [(30, 39), (40, 49), (50, 59), (60, 69), (70, 84), (85, 99), (100, 199), (200, 1000)]
for field, label, offset, color in [("on", "Aperture", -.19, "#17577f"),
                                    ("off", "Scaled annulus", .19, "#c45c28")]:
    values = []
    for lo, hi in bands:
        selected = [r for r in rows if r["period"] == "summer"
                    and (int(r["pi_min"]), int(r["pi_max_inclusive"])) == (lo, hi)]
        values.append(sum(int(r[field]) * (float(r["alpha"]) if field == "off" else 1)
                          for r in selected))
    sp.bar(np.arange(len(bands)) + offset, values, width=.36, color=color, label=label)
sp.set_xticks(range(len(bands)), ["0.3–0.4", "0.4–0.5", "0.5–0.6", "0.6–0.7", "0.7–0.85", "0.85–1", "1–2", "2–10"],
              rotation=45, ha="right")
sp.set(xlabel="Recorded energy band (keV)", ylabel="Counts per band", title="June–September 2021")
for ax in (lc, sp):
    ax.spines[["top", "right"]].set_visible(False)
fig.savefig(here / "swift-history.png", dpi=150)

"""Render the 2032 Phoebe prediction from the published five-minute light curve."""

import csv
from datetime import datetime, timezone

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from search import DIRECTORY

with (DIRECTORY / "phoebe-lightcurve.csv").open() as f:
    rows = list(csv.DictReader(f))
times = [
    datetime.fromisoformat(r["earth_utc"]).replace(tzinfo=timezone.utc) for r in rows
]


def column(name):
    return np.array([float(r[name]) for r in rows])


plt.rcParams.update(
    {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
)
fig, ax = plt.subplots(figsize=(10, 4.4), layout="constrained")
a, b = column("cassini_aldebaran"), column("cassini_gacrux")
ax.fill_between(
    times,
    np.minimum(a, b),
    np.maximum(a, b),
    color="#1a8884",
    alpha=0.24,
    label="Range of two Cassini ring profiles",
)
ax.plot(
    times,
    column("planet_and_rings"),
    color="#006a67",
    lw=2.4,
    label="Voyager ring profile",
)
ax.plot(
    times,
    column("planet_only"),
    color="#b56529",
    lw=1.6,
    ls="--",
    label="Saturn globe alone",
)
ax.set(
    title="Saturn’s rings dim Phoebe on both sides of totality",
    ylabel="Fraction of direct sunlight reaching Phoebe",
    xlabel="Earth reception time (UTC), March 2032",
    ylim=(-0.04, 1.08),
)
ax.set_xlim(times[0], times[-1])
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M", tz=timezone.utc))
ax.grid(axis="y", alpha=0.16)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.48), ncol=3, frameon=False)
fig.savefig(DIRECTORY / "phoebe.png", dpi=160)

"""Plot the eclipse and the quadratic-ephemeris timing comparison."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from measure import DATA, OUT, PERIOD, REF, ROOT, fit, gaia_delay

colors = ["#224d7a", "#277464", "#ab6625", "#855395", "#b3354e"]
sectors = [12, 39, 66, 103, 104]
curves = ROOT / "raw" / "curves.csv"
if curves.exists():
    lc = pd.read_csv(curves)
    rows = []
    for s, a in lc.groupby("sector"):
        ph = ((a.time - REF) / PERIOD) % 1
        # Each sector uses its own common eclipse phase for the shape comparison only.
        from measure import clock

        shift = clock(a.time.to_numpy(), a.flux.to_numpy())[0] / (PERIOD * 86400)
        phase = (ph - shift) % 1
        for i in range(100):
            f = a.flux[(phase >= i / 100) & (phase < (i + 1) / 100)]
            rows.append((s, (i + 0.5) / 100, f.median()))
    pd.DataFrame(rows, columns=["sector", "phase", "flux"]).to_csv(
        OUT / "folds.csv", index=False, float_format="%.9f"
    )
fold = pd.read_csv(OUT / "folds.csv")
times = pd.read_csv(OUT / "timings.csv").dropna(subset=["both"])
t = times.time.to_numpy()
beta, _, _ = fit(times, quadratic=True)
x = t - t.mean()
base = beta[0] + beta[1] * x + beta[2] * (x / 1000) ** 2
residual = times.both.to_numpy() - base
fig, (ax, later) = plt.subplots(2, 1, figsize=(8.6, 6.6), layout="constrained")
for s, c in zip(sectors, colors):
    a = fold[fold.sector == s]
    ax.plot(a.phase, a.flux, color=c, lw=1.5, label=f"S{s}", alpha=0.86)
ax.legend(
    ncol=5, fontsize=8, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.27)
)
ax.set(
    xlabel="Inner orbital phase (minima aligned by sector)",
    ylabel="Normalized TESS flux",
    xlim=(0, 1),
)
ax.set_title(
    "TIC 96100077  |  Two eclipses every 7.93 hours", loc="left", fontweight="bold"
)
outer_zero = (
    Time(2016.0, format="jyear", scale="tcb").tdb.jd
    - 2457000
    + float(pd.read_csv(DATA / "gaia-orbit.csv").iloc[0].t_periastron)
)
period = float(pd.read_csv(DATA / "gaia-orbit.csv").iloc[0].period)
phase = ((t - outer_zero) / period) % 1
grid = np.linspace(0, 1, 401)
later.plot(
    grid,
    beta[-1] * gaia_delay(outer_zero + period * grid),
    color="#243444",
    lw=1.6,
    label="Gaia orbit × fitted amplitude",
)
for s, c in zip(sectors, colors):
    mask = times.sector == s
    later.scatter(
        phase[mask], residual[mask], color=c, s=26, edgecolor="white", linewidth=0.5
    )
later.set(
    xlabel="Gaia outer orbital phase from periastron",
    ylabel="Eclipse timing residual (s)",
    xlim=(0, 1),
)
later.legend(frameon=False, fontsize=9, loc="lower right")
later.set_title(
    "Timing after a fitted quadratic ephemeris", loc="left", fontweight="bold"
)
for axis in (ax, later):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", alpha=0.2)
fig.savefig(ROOT / "timing.png", dpi=175)

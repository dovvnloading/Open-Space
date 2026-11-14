"""Show the conditional velocity forecasts and geometrically usable nights at La Silla."""

from pathlib import Path
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time

from orbit import CLOCK_SCALE, NOMINAL, sample_orbits, velocity

ROOT = Path(__file__).resolve().parent


if __name__ == "__main__":
    dates = pd.date_range("2027-02-01", "2027-10-01", freq="D")
    jd = Time(dates.to_numpy(), scale="utc").tcb.jd
    elements = sample_orbits()
    v = velocity(jd[:, None], elements)
    g = elements[:, 2]
    fig, (ax, alt) = plt.subplots(
        2, 1, figsize=(8.7, 4.8), height_ratios=[3, 1], sharex=True
    )
    fig.subplots_adjust(left=0.10, right=0.975, top=0.88, bottom=0.10, hspace=0.14)
    fig.suptitle(
        "TIC 96100077  |  Where spectroscopy can test the triple",
        x=0.10,
        ha="left",
        fontsize=14,
    )
    colors = ["#145467", "#b45531"]
    for y, label, color in zip(
        [v, g + CLOCK_SCALE * (v - g)],
        ["Gaia RV tracer", "Eclipsing pair’s centre of mass"],
        colors,
    ):
        lo, mid, hi = np.quantile(y, [0.1, 0.5, 0.9], axis=1)
        ax.fill_between(dates, lo, hi, color=color, alpha=0.13, lw=0)
        ax.plot(dates, mid, color=color, label=label, lw=2.2)
    ax.axhline(NOMINAL[2], color="#727878", ls=":", lw=1)
    ax.set_ylabel("Barycentric RV (km/s)")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.text(
        0.99,
        0.04,
        "Bands show the formal 10th–90th percentiles",
        ha="right",
        color="#525b61",
        fontsize=8,
        transform=ax.transAxes,
    )
    nights = pd.read_csv(ROOT / "nightly-visibility.csv")
    alt.fill_between(
        pd.to_datetime(nights.evening_utc_date),
        nights.visible_hours,
        color="#80959e",
        alpha=0.63,
    )
    alt.set_ylabel("Dark hours\nabove 35°", fontsize=9)
    windows = pd.read_csv(ROOT / "quadrature-windows.csv")
    for role, d in windows.groupby("role"):
        t = pd.Timestamp(d.week_start.iloc[0])
        for a in (ax, alt):
            a.axvspan(
                t, t + pd.Timedelta(days=7), color="#33464b", alpha=0.10, ec="none"
            )
    alt.set_xlim(dates[0], dates[-1])
    alt.xaxis.set_major_locator(mdates.MonthLocator())
    alt.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    alt.set_xlabel("2027 · La Silla, Chile")
    for a in (ax, alt):
        a.spines[["right", "top"]].set_visible(False)
        a.tick_params(labelsize=9)
    fig.savefig(ROOT / "velocity-forecast.png", dpi=150)

"""Regenerate frequency.png from the two published timing tables."""
import csv
from datetime import date, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.ticker import ScalarFormatter

from replay import ROOT, load


def day(mjd):
    return dates.date2num(date(1858,11,17)) + mjd


rows = load("gbm.csv")
with (ROOT / "nustar.csv").open() as stream:
    nustar = list(csv.DictReader(stream))
fig,axes = plt.subplots(1,3,figsize=(12.5,3.85),width_ratios=[1.05,1.1,1.15])
for ax,title,limits,scale in zip(axes,
        ["18 years of monitoring","Spin-down weakens","Return to spin-down in 2026"],
        [(54600,61430),(60240,61390),(61064,61370)], [1e3,1e6,1e6]):
    subset = [r for r in rows if limits[0] <= r["barytime"] <= limits[1]]
    ts = [day(r["barytime"]) for r in subset]
    ax.errorbar(ts,[r["frequency"]*scale for r in subset],
                yerr=[r["frequency_err"]*scale for r in subset],fmt=".",ms=4,
                elinewidth=.8,color="#234c70",label="Fermi/GBM")
    ax.set_xlim(*map(day,limits))
    ax.xaxis_date()
    ax.set_title(title,fontsize=10,loc="left",weight="bold")
    ax.grid(axis="y",alpha=.22)
    ax.spines[["top","right"]].set_visible(False)
    fmt=ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(fmt)
    ax.set_ylabel("Pulse frequency (mHz)" if scale==1e3 else "Pulse frequency (µHz)")
axes[0].xaxis.set_major_locator(dates.YearLocator(5))
axes[0].xaxis.set_major_formatter(dates.DateFormatter("%Y"))
axes[1].xaxis.set_major_locator(dates.MonthLocator(bymonth=[1,7]))
axes[1].xaxis.set_major_formatter(dates.DateFormatter("%b\n%Y"))
axes[2].xaxis.set_major_locator(dates.MonthLocator(interval=2))
axes[2].xaxis.set_major_formatter(dates.DateFormatter("%b"))
axes[2].set_xlabel("2026")
axes[1].scatter([day(float(r["mjd_tt"])) for r in nustar],
                [float(r["frequency_uHz"]) for r in nustar],marker="D",s=20,
                color="#a73712",label="NuSTAR A/B",zorder=4)
axes[1].legend(frameon=False,fontsize=8)
axes[2].axvspan(day(61111),day(61176),alpha=.13,color="#e78e3f")
axes[2].axvspan(day(61296),day(61320),alpha=.12,color="#278a68")
fig.tight_layout()
fig.savefig(ROOT/"frequency.png",dpi=180)

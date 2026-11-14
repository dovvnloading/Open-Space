"""Plot the two independent measurements of the 4U 1538−522 timing change."""
import csv
from datetime import date
import json
from pathlib import Path

import matplotlib.dates as dates
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

root=Path(__file__).resolve().parent
with (root/"gbm.csv").open() as stream:
    rows=list(csv.DictReader(stream))
t=np.array([float(r["barytime"]) for r in rows])
f=np.array([float(r["frequency"])*1e6 for r in rows])
e=np.array([float(r["frequency_err"])*1e6 for r in rows])
m=json.loads((root/"maxi_results.json").read_text())["61132-61325"]
dayzero=dates.date2num(date(1858,11,17))
fig,axes=plt.subplots(1,2,figsize=(9.5,3.5),width_ratios=[1.15,1])
for ax,(lo,hi),title in zip(axes,[(54640,61430),(60830,61380)],
                             ["Fermi pulse-frequency history","Two instruments measure the 2026 decline"]):
    i=(t>=lo)&(t<=hi)
    ax.errorbar(t[i]+dayzero,f[i],yerr=e[i],fmt=".",color="#244b6a",elinewidth=.85,
                ms=4,label="Fermi/GBM")
    ax.set_xlim(lo+dayzero,hi+dayzero)
    ax.set_title(title,loc="left",weight="bold",fontsize=10)
    ax.set_ylabel("Pulse frequency (µHz)")
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.grid(axis="y",alpha=.2)
    ax.spines[["right","top"]].set_visible(False)
    ax.xaxis_date()
axes[0].xaxis.set_major_locator(dates.YearLocator(4))
axes[0].xaxis.set_major_formatter(dates.DateFormatter("%Y"))
axes[1].xaxis.set_major_locator(dates.MonthLocator(interval=3))
axes[1].xaxis.set_major_formatter(dates.DateFormatter("%b\n%Y"))
endpoints=np.array([61132,61325])
params=m["best_fundamental"]
axes[1].plot(endpoints+dayzero,params["frequency_uHz"]+(endpoints-m["epoch_mjd"])*params["slope_uHz_per_year"]/365.25,
             color="#c44e16",lw=2,label="MAXI/GSC pulse track")
axes[1].legend(frameon=False,fontsize=8)
fig.tight_layout()
fig.savefig(root/"spin.png",dpi=175)

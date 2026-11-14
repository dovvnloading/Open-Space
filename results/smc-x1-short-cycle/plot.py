"""Plot the independent measurements and the five minima used for the four-cycle count."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as md
import matplotlib.pyplot as plt
import numpy as np
from replay import ROOT, WINDOWS, data, fit, design, select

series = data()
results = json.loads((ROOT / "measurements.json").read_text())
origin = md.datestr2num("1858-11-17")
fig, (top, bottom) = plt.subplots(2, 1, figsize=(9.2, 5.8), gridspec_kw={"height_ratios": [1.45, 1]})
for name, mission, color, marker in [("smcbat", "Swift/BAT", "#273d75", "o"),
                                     ("smcmaxi", "MAXI/GSC", "#b44c21", "s")]:
    a = select(series[name], *WINDOWS["fast"])
    # Scale each instrument separately only for the overlay.
    center = np.median(a[:, 1])
    scale = np.percentile(a[:, 1], 90) - np.percentile(a[:, 1], 10)
    for i in range(0, 165, 4):
        rows = select(a, 60695+i, 60699+i)
        if len(rows):
            w = 1/rows[:, 2]**2
            top.errorbar(np.average(rows[:, 0], weights=w)+origin, (np.average(rows[:, 1], weights=w)-center)/scale,
                         1/np.sqrt(w.sum())/scale, fmt=marker, color=color, ms=3, alpha=.72, lw=.8)
    p = results["timing"][name]["fast"]["period_d"]
    _, b = fit(a, p)
    t = np.arange(60695,60860,.2)
    top.plot(t+origin, (design(t,p) @ b - center)/scale, color=color, lw=1.5, label=mission)
    for j, epoch in enumerate(WINDOWS):
        item = results["timing"][name][epoch]
        if "period_d" in item:
            lo,hi=WINDOWS[epoch]
            bottom.plot([lo+origin,hi+origin],[item["period_d"]]*2,color=color,lw=1,alpha=.55)
            bottom.plot((lo+hi)/2+origin,item["period_d"],marker=marker,color=color,ms=6)
for t in results["minima"]["Swift"]["mjd"]:
    top.axvline(t+origin, ls=":", lw=.8, color="#898989", zorder=0)
top.text(0,1.04,"Four cycles, January–June 2025",transform=top.transAxes,weight="bold")
bottom.text(0,1.07,"Cycle length before, during and after the short-period episode",
            transform=bottom.transAxes,weight="bold")
top.set_ylabel("Relative X-ray flux")
bottom.set_ylabel("Cycle length (days)")
top.legend(ncol=2,loc="upper right",frameon=False)
bottom.set_ylim(35,63)
top.xaxis.set_major_formatter(md.DateFormatter("%b %d"))
bottom.xaxis.set_major_locator(md.MonthLocator(interval=4))
bottom.xaxis.set_major_formatter(md.DateFormatter("%b %Y"))
for ax in (top,bottom):
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y",alpha=.17)
fig.subplots_adjust(left=.1,right=.98,bottom=.09,top=.93,hspace=.53)
fig.savefig(ROOT / "cycle-timing.png",dpi=160,metadata={"Software":"Open-Space"})

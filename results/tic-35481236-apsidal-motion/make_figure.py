"""Plot archived eclipse phases and same-cadence TESS shape measurements."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

root = Path(__file__).resolve().parent
def read(name):
    with (root/name).open() as stream:
        return list(csv.DictReader(stream))
def year(jd):
    return Time(jd, format="jd").decimalyear

sectors = read("sectors.csv")
timings = [r for r in sectors if r["product"]=="SPOC PDCSAP_FLUX" or r["sector"]=="33"]
spoc = [r for r in timings if r["product"]=="SPOC PDCSAP_FLUX"]
ground = [r for r in read("ground.csv") if r["sample"]!="all"]
fig,axes = plt.subplots(1,2,figsize=(9.1,3.65),gridspec_kw={"width_ratios":[1.4,1]})
blue,orange = "#0e5870","#a2472d"
for row in ground:
    y = float(row["phase"])
    err = [[y-float(row["bootstrap_low"])],[float(row["bootstrap_high"])-y]]
    axes[0].errorbar(year(float(row["median_hjd"])),y,yerr=err,fmt="s",color=orange,capsize=3,
                     label="ASAS-3 half sample" if row["sample"]=="early" else None)
axes[0].plot([year(float(r["median_btjd"])+2457000) for r in timings],[float(r["phase"]) for r in timings],
             "o",color=blue,label="TESS sector")
axes[0].set_ylabel("Secondary eclipse phase")
axes[0].legend(frameon=False,fontsize=8)
for kind, color, marker in [("primary",blue,"o"),("secondary",orange,"s")]:
    axes[1].plot([year(float(r["median_btjd"])+2457000) for r in spoc],
                 [float(r[kind+"_fwhm_min"]) for r in spoc],marker+"-",linewidth=1.25,color=color,label=kind.capitalize())
axes[1].set_ylabel("Eclipse FWHM (minutes)")
axes[1].legend(frameon=False,fontsize=8)
for ax in axes:
    ax.set_xlabel("Year")
    ax.grid(axis="y",alpha=.2)
    ax.spines[["top","right"]].set_visible(False)
axes[0].set_title("When the secondary eclipse occurs",loc="left",fontsize=10)
axes[1].set_title("Widths at 2-minute cadence",loc="left",fontsize=10)
fig.tight_layout()
fig.savefig(root/"timing.png",dpi=160)

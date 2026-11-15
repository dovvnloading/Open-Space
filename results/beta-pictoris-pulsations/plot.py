"""Plot the measured amplitude ratio against observation time."""
import csv
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
ORIGIN = datetime(2014, 12, 8, 12)

def main():
    with (ROOT/"measurements.csv").open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["flux"] == "PDCSAP_FLUX"]
    start = np.array([float(r["btjd_start"]) for r in rows])
    end = np.array([float(r["btjd_end"]) for r in rows])
    time = mdates.date2num([ORIGIN + timedelta(days=float(v)) for v in (start+end)/2])
    width = (end-start)/2
    low, high, ratio = ([float(r[k]) for r in rows]
                        for k in ("low_ppm", "high_ppm", "doublet_band_ratio"))
    fig, axes = plt.subplots(2,1,figsize=(8.6,5.6),height_ratios=(1.4,1),
                             sharex=True,layout="constrained")
    colors = ["#00656e","#a43c25"]
    for ax in axes:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.18)
    for value,color,label,marker in [(low,colors[0],"47.28 cycles/day group","o"),
                                      (high,colors[1],"47.44 cycles/day peak","s")]:
        axes[0].errorbar(time,value,xerr=width,fmt=marker,color=color,
                         markersize=4.8,elinewidth=1.2,label=label)
    axes[0].set(ylabel="Semi-amplitude (ppm)",ylim=(0,1200))
    axes[0].set_title("β Pictoris changed its dominant pulsation in the 47 cycles/day band",
                      loc="left",fontsize=13)
    axes[0].legend(frameon=False,loc="upper right",fontsize=9)
    axes[1].axhline(1,lw=1,ls="--",color="0.35")
    axes[1].errorbar(time,ratio,xerr=width,fmt="o",color="#2c3548",markersize=4.8)
    axes[1].set(ylabel="47.28 group / 47.44 peak\nRMS amplitude ratio",ylim=(0,2.65),
                xlabel="Observation date · horizontal bars show each sector’s time span")
    for row,when,value in zip(rows[-3:],time[-3:],ratio[-3:]):
        axes[1].annotate("S"+row["sector"],xy=(when,value),xytext=(-7,7),
                         textcoords="offset points",fontsize=8,ha="right")
    axes[1].xaxis.set_major_locator(mdates.YearLocator())
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.savefig(ROOT/"amplitude-reversal.png",dpi=145,bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()

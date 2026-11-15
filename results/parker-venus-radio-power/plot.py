"""Draw both dipoles against the Venus sunlight boundary and report solar-center contacts."""
import csv
import json
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from analyze import ANTENNAS, HERE, load_minutes


def load_sightlines():
    dates = {}
    with (HERE/"sightlines.csv").open() as f:
        for row in csv.DictReader(f):
            dates.setdefault(row["utc"][:10],[]).append(row)
    return dates


def solar_crossings(rows):
    times = [datetime.fromisoformat(r["utc"]) for r in rows]
    solar = [float(r["solar_center_clearance_deg"]) for r in rows]
    crossings = []
    for i,(a,b) in enumerate(zip(solar,solar[1:])):
        if (a < 0) != (b < 0):
            crossings.append(times[i] + (times[i+1]-times[i])*(-a/(b-a)))
    return crossings


def main():
    all_rows = load_sightlines()
    series = load_minutes(HERE/"minutes.csv.gz")
    manifest = json.loads((HERE/"manifest.json").read_text())
    with (HERE/"solar-contacts.csv").open("w",newline="") as fp:
        out = csv.writer(fp)
        out.writerow(["date","apparent_radius_max_15s_deg","maximum_grid_utc","sun_center_ingress_utc","sun_center_egress_utc"])
        for day,rows in all_rows.items():
            peak = max(rows,key=lambda r:float(r["venus_radius_deg"]))
            contacts = [t.isoformat(timespec="seconds") for t in solar_crossings(rows)]
            out.writerow([day,peak["venus_radius_deg"],peak["utc"],*(contacts or ["",""])])

    fig,axes=plt.subplots(2,1,figsize=(8.6,5.8),layout="constrained")
    for ax,day in zip(axes,["2024-11-06","2020-07-11"]):
        center = next(m["maximum_radius_utc_minute"] for m in manifest if m["date"]==day)
        midnight=datetime.fromisoformat(day)
        mins=np.arange(center-24,center+24)
        times=[midnight+timedelta(minutes=float(i)+.5) for i in mins]
        for ant,color,label in zip(ANTENNAS,["#006b75","#933e63"],["V1–V2","V3–V4"]):
            data=series[(day,ant,40)][:,2]
            offsets=np.arange(1440)+.5-center
            pre=np.nanmedian(data[(offsets>=-10)&(offsets< -4)])
            post=np.nanmedian(data[(offsets>4)&(offsets<=10)])
            power=data[mins]-(pre+post)/2
            ax.plot(times,power,marker=".",ms=3.5,lw=1.3,color=color,label=label)
        a,b=solar_crossings(all_rows[day])
        ax.axvspan(a,b,facecolor="#ded8c6",alpha=.55,label="Solar center behind Venus",zorder=0)
        ax.axhline(0,color="#8e9294",lw=.65,ls="--")
        ax.set_xlim(times[0],times[-1]);ax.set_title(day,loc="left",weight="bold",fontsize=11)
        ax.set_ylabel("Relative power (dB)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"));ax.grid(axis="y",alpha=.17)
        for side in ["top","right"]:ax.spines[side].set_visible(False)
    axes[0].legend(frameon=False,ncols=3,fontsize=8,loc="upper right");axes[-1].set_xlabel("Spacecraft reception time (UTC)")
    fig.suptitle("Parker radio power at 4–8 MHz during two Venus flybys",fontsize=15,weight="bold")
    fig.savefig(HERE/"radio-power.png",dpi=135)
    plt.close(fig)


if __name__=="__main__":
    main()

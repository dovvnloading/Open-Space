"""Stack the fully paired non-SGR 1935+2154 starts on the published clock.

From the catalog-clock directory run:
  python plot.py /path/to/clock-cache --output .
Requires NumPy, Astropy and Matplotlib. The cache is populated by measure.py.
"""
import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

SHIFT = 69.184


def tt_epoch(header):
    if header["TIMESYS"].strip() != "TT":
        raise ValueError("Unexpected FITS time system")
    return Time(header["MJDREFI"], header["MJDREFF"], format="mjd", scale="tt")


def intervals(starts, stops):
    """Union GTIs, handling abutting or overlapping intervals once."""
    merged = []
    for start, stop in sorted(zip(starts, stops)):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], stop)
        else:
            merged.append([start, stop])
    return merged


def bin_data(t, good_times, center, edges):
    """Count only events in overlapped GTIs and integrate exact live wall time."""
    counts = np.zeros(len(edges)-1, dtype=int)
    exposure = np.zeros(len(counts))
    lower, upper = center + edges[[0,-1]]
    for start, stop in good_times:
        lo, hi = max(lower, start), min(upper, stop)
        if hi <= lo:
            continue
        i, j = np.searchsorted(t, [lo, hi])
        counts += np.histogram(t[i:j] - center, bins=edges)[0]
        exposure += np.maximum(0., np.minimum(edges[1:],hi-center) -
                              np.maximum(edges[:-1],lo-center))
    return counts, exposure


def extract(catalog, cache):
    with (catalog / "windows.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["source"] != "SGR 1935+2154" and
                r["narrow_literal"] != "" and r["narrow_plus69"] != ""]
    by_obs = defaultdict(list)
    for r in rows:
        by_obs[r["obsid"]].append(r)
    with (catalog / "inputs.csv").open() as f:
        hashes = {r["obsid"]:r["events_sha256"] for r in csv.DictReader(f)}
    # Put both tested hypotheses at the center of their respective zoom bins.
    grids = {"overview": np.arange(-20.25,100.26,0.5),
             "published": np.arange(-0.505,1.5051,0.01),
             "corrected": SHIFT+np.arange(-0.505,1.5051,0.01)}
    values = {key:(np.zeros(len(edges)-1,int), np.zeros(len(edges)-1))
              for key,edges in grids.items()}
    for obs, selected in by_obs.items():
        path = cache/f"ni{obs}_0mpu7_cl.evt.gz"
        with path.open("rb") as f:
            if hashlib.file_digest(f,"sha256").hexdigest() != hashes[obs]:
                raise ValueError(f"Input hash differs for {obs}")
        with fits.open(path) as h:
            ev = h["EVENTS"]
            keep = (ev.data["PI"] >= 50) & (ev.data["PI"] <= 800)
            t = np.sort(ev.data["TIME"][keep] + ev.header.get("TIMEZERO",0.))
            gti = h["GTI"]
            offset = gti.header.get("TIMEZERO",0.) + (tt_epoch(gti.header)-tt_epoch(ev.header)).sec
            good = intervals(gti.data["START"]+offset, gti.data["STOP"]+offset)
            for r in selected:
                for key,edges in grids.items():
                    c, e = bin_data(t, good, float(r["predicted_met_tt_s"]), edges)
                    counts, exposure = values[key]
                    counts += c
                    exposure += e
    return grids, values, len(rows), len(by_obs), len({r["source"] for r in rows})


def run(catalog, cache, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grids, values, nr, no, ns = extract(catalog, cache)
    output.mkdir(parents=True,exist_ok=True)
    profiles = {}
    with (output/"profile.csv").open("w",newline="") as f:
        writer = csv.writer(f,lineterminator="\n")
        writer.writerow(["view","offset_lo_s","offset_hi_s","photon_counts","summed_gti_seconds","rate_per_s"])
        for key,(counts,expo) in values.items():
            edges = grids[key]
            rate = np.divide(counts,expo,out=np.full_like(expo,np.nan),where=expo>0)
            mid = (edges[1:]+edges[:-1])/2
            profiles[key] = (mid,rate)
            peak = int(np.nanargmax(rate))
            print(f"{key} maximum at {mid[peak]:.3f} s: {counts[peak]} counts / "
                  f"{expo[peak]:.4f} summed GTI seconds = {rate[peak]:.2f} counts/s")
            for lo, hi, n, e, v in zip(edges[:-1],edges[1:],counts,expo,rate):
                writer.writerow([key,f"{lo:.4f}",f"{hi:.4f}",n,f"{e:.6f}",f"{v:.6f}"])

    # The catalog number of starts and exposure vary independently of the ordinate.
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10.5,
                         "axes.spines.right":False,"axes.spines.top":False})
    fig, axes = plt.subplots(1,2,figsize=(11,4),gridspec_kw={"width_ratios":[1.35,1]})
    fig.subplots_adjust(left=.075,right=.98,top=.77,bottom=.23,wspace=.33)
    fig.suptitle("NICER photons align with the shifted catalog times",
                 x=.075,y=.965,ha="left",fontsize=17,fontweight="bold",color="#152b3d")
    fig.text(.075,.852,f"{nr} catalog starts  ·  {no} observations  ·  {ns} sources  ·  0.5–8 keV",
             color="#456373",fontsize=10.5)
    blue,orange="#006c83","#b94713"
    a,b=axes
    a.stairs(profiles["overview"][1],grids["overview"],baseline=None,color="#233b52",lw=1.25)
    a.axvline(0,color=blue,ls="--",lw=1.2)
    a.axvline(SHIFT,color=orange,ls="--",lw=1.2)
    a.set_xlim(-20,100)
    a.set_title("Full interval · 0.5 s bins",loc="left",fontsize=11.5,fontweight="bold")
    a.set_ylabel("NICER count rate (counts/s)")
    a.set_xlabel("Seconds from published start")
    a.set_ylim(bottom=0)
    for ax in axes:
        ax.grid(axis="y",color="#c9d4db",lw=.7,alpha=.65)
        ax.set_axisbelow(True)
    for key,shift,color,label in [("published",0,blue,"Published"),("corrected",SHIFT,orange,"+69.184 s")]:
        b.stairs(profiles[key][1],grids[key]-shift,baseline=None,color=color,lw=1.3,label=label)
    b.axvline(0,ls="--",lw=.9,color="#53636a")
    b.set_xlim(-.40,1.00)
    b.set_ylim(bottom=0)
    b.set_title("Aligned starts · 10 ms bins",loc="left",fontsize=11.5,fontweight="bold")
    b.set_xlabel("Seconds from each proposed start")
    b.set_ylabel("NICER count rate (counts/s)")
    b.legend(frameon=False,fontsize=9.7,loc="upper right")
    fig.text(.075,.038,"Each bin divides summed photon counts by summed valid observing time. A photon can contribute to overlapping burst windows.",
             fontsize=9.0,color="#425561")
    fig.savefig(output/"time-profile.png",dpi=180,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache",type=Path)
    p.add_argument("--catalog",type=Path,default=Path(__file__).resolve().parent)
    p.add_argument("--output",type=Path,default=Path("."))
    args=p.parse_args()
    run(args.catalog,args.cache,args.output)

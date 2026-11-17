"""Plot one held-out target and the signed data/template measurements for every trigger."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from desi_shard import targets
from lithium_profile import LI, profile


def plot(cache, ledger, output):
    triggers = {}
    for p in ledger.glob("scan-*.jsonl"):
        for line in p.open():
            row = json.loads(line)
            m = row["meta"]
            if any(abs(c["rest_A"]-LI)<1 for arm in row["arms"].values() for c in arm["candidates"]):
                triggers[str(m["targetid"])] = m
    # First held-out target in ascending TARGETID order.
    tid = min(t for t,m in triggers.items() if m["healpix"]==49134)
    meta, arms = next((m,a) for m,a in targets(cache, "main", "bright",49134) if str(m["targetid"])==tid)
    fig, (left,right) = plt.subplots(1,2,figsize=(9.6,4.1),layout="constrained",
                                    gridspec_kw={"width_ratios":[1.25,1]})
    colors = ["#24698c","#af5227","#3a7666"]
    for i,e in enumerate(arms["r"]):
        w = e.wave/e.doppler
        use = (w>LI-8)&(w<LI+8)
        meas = profile(e,LI,15)
        off = 0.27*i
        left.plot(w[use],e.flux[use]/meas["observed_cont17"]+off,marker=".",ms=3,
                  lw=1.1,color=colors[i],label=f"EXPID {e.expid}")
        left.plot(w[use],e.template[use]/meas["template_cont17"]+off,"--",lw=1.25,color=colors[i])
    left.axvline(LI, color="#6c7476",lw=.7,alpha=.65)
    left.set(xlim=(LI-8,LI+8),xlabel="Stellar-rest vacuum wavelength (Å)",
             ylabel="Local continuum-normalized flux + offset",title="A held-out star on three nights")
    left.legend(fontsize=8,loc="lower left",frameon=False)
    rows=[r for p in (output.parent/"measurements").glob("lithium-*.csv")
          for r in csv.DictReader(p.open()) if r["targetid"] in triggers and r["rest_A"]==str(LI)
          and r["half_window_px"]=="15" and r["observed_snr"]]
    for region,color,label in [(False,"#24698c","Initial fields"),(True,"#af5227","Held-out field")]:
        data = [r for r in rows if (triggers[r["targetid"]]["healpix"]==49134)==region]
        right.scatter([float(r["template_ew_A"]) for r in data],[float(r["observed_ew_A"]) for r in data],
                      s=27, color=color, alpha=.8,edgecolors="white",lw=.45,label=label)
    right.axhline(0,color="#555",lw=.8)
    right.set(xlabel="Model line strength (Å, absorption < 0)",ylabel="Observed line strength (Å)",
              title="One point per exposure")
    right.legend(frameon=False,fontsize=8)
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.title.set_fontsize(12)
    fig.savefig(output,dpi=155)
    print(tid,len(rows))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache",type=Path)
    p.add_argument("ledger",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    plot(args.cache,args.ledger,args.output)

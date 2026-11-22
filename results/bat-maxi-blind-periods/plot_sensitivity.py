"""Plot measured joint recovery at the tested periods and amplitudes."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"]="open-space-bat-maxi-recovery"
import matplotlib.pyplot as plt
import numpy as np


def plot(folder,out):
    panels=[
        ("Sinusoidal modulation",folder/"sensitivity.csv","group","relative_semiamplitude","fraction",
         "source_bootstrap_lo","source_bootstrap_hi","Semiamplitude / band RMS"),
        ("Periodic outbursts (10% duty cycle)",folder/"review"/"pulses"/"pulse-sensitivity.csv","instrument",
         "peak_to_baseline_rms","recovery_fraction","source_bootstrap_025","source_bootstrap_975","Peak above baseline / band RMS")]
    fig,axes=plt.subplots(1,2,figsize=(9.4,4),sharey=True,layout="constrained")
    colors=["#2166ac","#b24c00","#268054","#8148a1"]
    markers=["o","s","^","D"]
    for ax,(title,path,group,akey,fkey,lkey,hkey,xlabel) in zip(axes,panels):
        with path.open() as handle:
            rows=[r for r in csv.DictReader(handle) if r[group]=="joint"]
        for i,period in enumerate(sorted({float(r["period_days"]) for r in rows})):
            vals=sorted((r for r in rows if float(r["period_days"])==period),key=lambda r:float(r[akey]))
            x=np.array([float(r[akey]) for r in vals])
            y=100*np.array([float(r[fkey]) for r in vals])
            lo=100*np.array([float(r[lkey]) for r in vals])
            hi=100*np.array([float(r[hkey]) for r in vals])
            ax.errorbar(x*1.026**(i-1.5),y,yerr=[np.maximum(0,y-lo),np.maximum(0,hi-y)],
                        color=colors[i],fmt=markers[i]+"-",linewidth=1.25,markersize=4,capsize=2,label=f"{period:g} d")
        ax.set_xscale("log")
        ax.set_xticks(x,labels=[f"{v:g}" for v in x])
        ax.minorticks_off()
        ax.set(title=title,xlabel=xlabel,ylim=(-3,104))
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.2)
    axes[0].set_ylabel("Recovered in both monitors (%)")
    axes[1].legend(loc="lower right",ncols=2,fontsize=8,framealpha=.9)
    fig.savefig(out,metadata={"Date":None})
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    plot(a.evidence,a.output)

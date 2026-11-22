"""Fold both instruments at one measured ephemeris for direct comparison of the halves."""
import argparse
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "open-space-bat-maxi"
import matplotlib.pyplot as plt
import numpy as np
from field_checks import load_archives, read_rows
from search import bounded_uncertainties


def draw(source, period, curves, path):
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.0), sharex=True, layout="constrained")
    for ax,band,unit in zip(axes,("bat","maxi"),("counts cm⁻² s⁻¹","photons cm⁻² s⁻¹")):
        t,y,err=curves[band]
        t=np.asarray(t)
        design_error,_=bounded_uncertainties(t,err)
        half=(t.min()+t.max())/2
        phase=((t-55000.)/period)%1.
        for color,label,offset,sel in (("#17659B","Earlier half",-.006,t<=half),("#BB4B24","Later half",.006,t>half)):
            xs,ys,es=[],[],[]
            for i in range(16):
                use=sel & (phase>=i/16) & (phase<(i+1)/16)
                if np.any(use):
                    w=design_error[use]**-2
                    xs.append((i+.5)/16+offset)
                    ys.append(np.average(y[use],weights=w))
                    es.append(np.sqrt(np.sum((w/w.sum())**2 * err[use]**2)))
            ax.errorbar(xs,ys,yerr=es,fmt="o-",markersize=3.2,capsize=1.3,lw=.9,alpha=.9,label=label,color=color)
        ax.set_ylabel(band.upper()+"\n"+unit)
        ax.grid(axis="y",alpha=.18)
        ax.spines[["right","top"]].set_visible(False)
    axes[0].legend(ncol=2,frameon=False,loc="best")
    axes[1].set(xlabel="Phase at the BAT period, measured from MJD 55000",xlim=(0,1))
    fig.suptitle(f"{source}    {period:.8g} days",fontsize=13,weight="bold")
    fig.savefig(path,metadata={"Date":None})
    plt.close(fig)


def main(folder, out):
    out.mkdir(exist_ok=True,parents=True)
    candidates=list(read_rows(folder/"survivors.jsonl"))
    ids={r["source_id"] for r in candidates}
    curves=load_archives(sorted(folder.glob("curves-*.tar.gz")),{sid+"."+b for sid in ids for b in ("bat","maxi")})
    index=["# Phase folds of the independent numerical survivors","",
           "Both instruments use the refined BAT period and MJD 55000 as phase zero. Each panel divides its own observations at the time midpoint. Bin means use the bounded analysis weights. Error bars propagate the original quoted uncertainties through those weights; the period decision also uses a correlated-noise model. The two panels have independent rate scales.","",
           "| Source | BAT period (days) | Fold |","|---|---:|---|"]
    for i,r in enumerate(candidates):
        sid=r["source_id"]
        name=f'{i+1:03d}-{re.sub(r"[^A-Za-z0-9.+_-]","_",sid)}.svg'
        draw(r["source"],r["period_days"],{b:curves[sid+"."+b] for b in ("bat","maxi")},out/name)
        index.append(f'| {r["source"]} | {r["period_days"]:.9g} | [Both instruments]({name}) |')
    (out/"index.md").write_text("\n".join(index)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    main(args.evidence,args.out)

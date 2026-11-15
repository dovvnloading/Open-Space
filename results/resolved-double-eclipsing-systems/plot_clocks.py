"""Draw red-band, time-split phase curves for the three candidates."""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
epochs = pd.read_csv(root / "ztf_epochs.csv.gz", dtype={"source_id":str,"oid":str})
validation = json.loads((root / "validation.json").read_text())
systems = [r for r in csv.DictReader((root / "four_recoveries.csv").open()) if r["previous_fezenko"]=="False"]
fig, axes = plt.subplots(3, 2, figsize=(9.7,7.7),layout="constrained")
for row,system in enumerate(systems):
    for col in range(2):
        ax=axes[row,col]
        sid=system[f"source_id{col+1}"]
        fit=validation[sid]["zr"]
        sel=epochs[(epochs.source_id==sid)&(epochs.oid==fit["oid"])&(epochs.catflags==0)
                   &(epochs.magerr>0)&(epochs.magerr<.1)].sort_values("hjd").drop_duplicates("hjd")
        t=sel.hjd.to_numpy()
        y=sel.mag.to_numpy()
        keep=abs(y-np.median(y))<max(1.5,8*np.median(abs(y-np.median(y))))
        t,y=t[keep],y[keep]
        phase=((t-t[0])/fit["own_period"])%1
        early=t<=(t[0]+t[-1])/2
        for mask,color,label,marker in [(early,"#2e668a","Earlier half","o"),(~early,"#bb492d","Later half","^")]:
            bins=np.linspace(0,1,26)
            m=(bins[1:]+bins[:-1])/2
            val=[np.median(y[mask&(phase>=lo)&(phase<hi)]) if np.any(mask&(phase>=lo)&(phase<hi)) else np.nan
                 for lo,hi in zip(bins[:-1],bins[1:])]
            ax.plot(m,np.array(val)-np.median(y[early]),marker=marker,lw=.9,ms=3.4,color=color,label=label)
        ax.invert_yaxis()
        ax.grid(alpha=.18)
        ax.set_title(f"Gaia {sid}\n{fit['own_period']:.6f} d",fontsize=10)
        if col==0:ax.set_ylabel("Relative r magnitude")
        if row==2:ax.set_xlabel("Orbital phase")
axes[0,0].legend(loc="best",fontsize=8)
fig.suptitle("Resolved ZTF eclipse clocks in three candidate quadruples",fontsize=14)
fig.savefig(root/"folded_clocks.png",dpi=145)

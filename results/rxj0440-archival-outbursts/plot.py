"""Render the two contemporary light curves around the 2010 pulsed flare."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
from datetime import datetime, timedelta
from replay import load, ROOT

fig, axes = plt.subplots(2, 1, figsize=(9, 5.4), sharex=True)
fig.subplots_adjust(top=.83, bottom=.13, left=.14, right=.97, hspace=.09)
fig.suptitle("RX J0440.9+4431 brightened before its reported first outburst",
             x=.14, y=.97, ha="left", fontsize=14, weight="bold")
fig.text(.14,.885,"The shaded interval marks the 204.903-second pulse detection by Fermi/GBM.",
         ha="left", fontsize=10, color="#374151")
date = lambda t: datetime(1858,11,17) + timedelta(days=float(t))
for ax,name,color,label in zip(axes,["bat","maxi"],["#276b93","#9b4535"],
                               ["Swift/BAT, 15–50 keV\ncounts cm⁻² s⁻¹","MAXI/GSC, 2–20 keV\nphotons cm⁻² s⁻¹"]):
    d=load(name)
    if name=="bat":
        d=d[d[:,3]==0]
    d=d[(d[:,0]>=55216)&(d[:,0]<55262)&(d[:,2]>0)]
    ax.axvspan(date(55236),date(55238),alpha=.28,color="#d3a72b",lw=0)
    ax.axhline(0, color="#88929d",lw=.7)
    ax.errorbar([date(t) for t in d[:,0]],d[:,1],d[:,2],fmt="o",ms=3.6,
                elinewidth=1,color=color,capsize=2)
    ax.set_ylabel(label,fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=9)
axes[1].xaxis.set_major_formatter(DateFormatter("%b %d"))
axes[1].set_xlabel("2010 UTC calendar date (archive day labels)")
fig.savefig(ROOT/"february.png",dpi=140,metadata={"Software":"Open-Space"})

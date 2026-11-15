"""Plot the released, calibrated measurements. Run in this result directory."""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p=Path(__file__).resolve().parent
p1=pd.read_csv(p/"ps1-measurements.csv").query("source == 'target' and band == 'i'").copy()
p1["day"]=np.floor(p1.mjd)
early=p1.groupby("day").dmag.median()
z=pd.read_csv(p/"ztf-forced-blocks.csv").query("band == 'ZTF_g'")
r=json.loads((p/"measurement.json").read_text())
colors=["#007E80","#8C4991"]
fig,(ax,az)=plt.subplots(2,1,figsize=(8.1,5.5),gridspec_kw={"height_ratios":[.8,1]},layout="constrained")
ax.plot(early.index,early,"o",c=colors[0],ms=6)
ax.axvspan(56280,56380,color=colors[0],alpha=.10)
ax.axvspan(56640,56745,color="gray",alpha=.13)
ax.set(ylabel="PS1 i (mag)",xlabel="Modified Julian Date",title="WD J0923+7326 · Pan-STARRS1 nightly photometry")
az.plot(z.t,z.m,".",ms=5.5,color=colors[1])
for k in range(6):
    t=r["anchor_mjd"][0]+k*r["period_days"]
    az.axvspan(t-56,t+56,alpha=.1,color=colors[1])
az.set(ylabel="ZTF g (mag)",xlabel="Modified Julian Date",title="ZTF 14-day medians · shaded windows use the 2018–19 timing")
for a in (ax,az):
    a.invert_yaxis()
    a.grid(axis="y",alpha=.18)
    a.spines[["right","top"]].set_visible(False)
fig.savefig(p/"light-curves.png",dpi=150)

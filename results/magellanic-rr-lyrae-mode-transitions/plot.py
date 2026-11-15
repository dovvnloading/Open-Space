"""Plot same-band mode amplitudes and fundamental periods, with observing spans."""
import json
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).parent
data = json.loads((root/"measurements.json").read_text())
checks = json.loads((root/"independent-checks.json").read_text())
day0 = datetime(1995,10,9,12)
plt.rcParams.update({"svg.fonttype":"none","svg.hashsalt":"rrlyr09022","font.size":10,
                     "axes.spines.top":False,"axes.spines.right":False})
fig, (top, bottom) = plt.subplots(2,1,figsize=(8.6,6.4),sharex=True,layout="constrained")


def bar(ax,window,center,low,high,label,color,marker):
    t1,t2=window["first"],window["last"]
    ax.errorbar(day0+timedelta(days=(t1+t2)/2),center,xerr=timedelta(days=(t2-t1)/2),
                yerr=[[max(0,center-low)],[max(0,high-center)]],color=color,marker=marker,
                capsize=2,label=label)


for key,w in data["windows"].items():
    if "I_" not in key:
        continue
    for i,(color,marker,label) in enumerate([("#19557d","o","0.360-day overtone"),
                                              ("#b0441e","s","0.485-day fundamental")]):
        bar(top,w,w["semiamp"][i],*w["night_bootstrap_95"][i],
            label if key.startswith("3") else None,color,marker)
top.legend(frameon=False,loc="upper left")
top.set(ylabel="OGLE I semi-amplitude (mag)",ylim=(0,.355),title="OGLE-LMC-RRLYR-09022")
for key,w in checks["windows"].items():
    low,mid,hi=[(x-.48525)*86400 for x in w["fundamental_period_night_interval"]]
    label={"I_5200_5750":"OGLE I","V_6485_6850":"OGLE V"}.get(key)
    band=key[0]
    bar(bottom,data["windows"]["4"+key],mid,low,hi,label,
        "#8a572d" if band=="V" else "#303e58","D" if band=="V" else "o")
g=data["gaia"]
center=(g["catalog_period"]-.48525)*86400
delta=1.96*g["catalog_period_error"]*86400
bar(bottom,g,center,center-delta,center+delta,"Gaia DR3 G","#397358","*")
bottom.set(ylabel="Period offset from 0.48525 d (s)",xlabel="Observation year")
bottom.legend(frameon=False,loc="upper left",ncol=3)
for ax in (top,bottom):
    ax.grid(axis="y",alpha=.18)
fig.savefig(root/"mode-history.svg",metadata={"Date":None})

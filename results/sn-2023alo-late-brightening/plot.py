"""Plot direct photometry and difference-flux measurements from the same exposures."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from analyze import INPUTS, OBJECTS, ROOT, WINDOWS

colors = {"zg":"#177658","zr":"#98452d"}
data = pd.read_csv(INPUTS / "ztf_field.csv.gz")
data = data[(data.catflags == 0) & (data.magerr < .15)]
alerts = pd.read_json(INPUTS / "ztf_alerts.json.gz")
alerts = alerts[(alerts.isdiffpos > 0) & (alerts.rb > .55)]
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9,5.8), layout="constrained")
for band, oid in OBJECTS.items():
    direct = data[data.oid == oid].copy()
    direct["mag"] -= direct[direct.mjd.between(*WINDOWS["pre"])].mag.median()
    weekly = direct.groupby(np.floor(direct.mjd / 7)).agg(mjd=("mjd","median"),mag=("mag","median"))
    axes[0].scatter(weekly.mjd,weekly.mag,s=12,alpha=.75,color=colors[band],label=band[-1])
    a = alerts[alerts.fid == (1 if band == "zg" else 2)]
    flux = 10**((23.9-a.magpsf)/2.5)
    axes[1].errorbar(a.mjd,flux, yerr=flux*a.sigmapsf*np.log(10)/2.5, fmt=".",
                     color=colors[band],alpha=.75,markersize=4,elinewidth=.65)
for ax in axes:
    for k in ["trough","late"]:
        ax.axvspan(*WINDOWS[k],color="#82939d",alpha=.15)
    ax.grid(axis="y",alpha=.18)
axes[0].invert_yaxis()
axes[0].legend(title="ZTF filter",loc="upper left",fontsize=8)
axes[0].set_ylabel("Magnitude relative to\npre-2023 median")
axes[0].set_title("Direct photometry, including the host galaxy",loc="left",fontsize=11)
axes[1].set_ylabel("Positive difference flux (µJy)")
axes[1].set_xlabel("Modified Julian date")
axes[1].set_title("Detections in ZTF difference images",loc="left",fontsize=11)
axes[1].set_xlim(58900,60900)
fig.suptitle("Late brightening at the position of SN 2023alo",ha="left",x=.01,fontsize=15)
fig.savefig(ROOT/"light-curve.png",dpi=160)

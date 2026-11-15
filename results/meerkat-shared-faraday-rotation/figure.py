"""Regenerate the comparison figure from the measured tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
target = np.genfromtxt(root / "J1057-5226.csv", delimiter=",", names=True)
dates = np.genfromtxt(root / "dates.csv", delimiter=",", names=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top":False, "axes.spines.right":False})
fig, (left, right) = plt.subplots(1, 2, figsize=(9.4,3.7), gridspec_kw={"width_ratios":[1.35,1]})
t = (target["mjd"] - 58849) / 365.25 + 2020
y = target["corrected_rm"]
err = target["formal_error"]
peer = target["other_pulsars_offset"]
zero = np.average(y, weights=err**-2)
series = [(y-zero, "#245784", "Published ionosphere correction"), (y-peer-zero, "#ae5628", "After subtracting nearby pulsars")]
for v,color,label in series:
    left.scatter(t, v, color=color, s=10, alpha=.47)
    coef = np.polyfit(t, v, 1, w=1/err)
    xx = np.array([t.min(),t.max()])
    left.plot(xx, np.polyval(coef, xx), color=color, lw=2, label=label)
left.set(xlabel="Year",ylabel="Relative rotation measure (rad m⁻²)",title="PSR J1057−5226 · 53 matched visits")
left.legend(fontsize=8,frameon=False,loc="upper left")
right.axline((0,0),slope=1,lw=1,ls="--",color="#acacac")
right.scatter(dates["resid_A"],dates["resid_B"],s=14,color="#245784",alpha=.7)
right.set(xlabel="Pulsar set A (rad m⁻²)",ylabel="Pulsar set B (rad m⁻²)",title="Daily residuals · 124 dates",aspect="equal")
lower=min(dates["resid_A"].min(),dates["resid_B"].min())-.08
upper=max(dates["resid_A"].max(),dates["resid_B"].max())+.08
right.set(xlim=(lower,upper),ylim=(lower,upper))
right.text(.05,.94,"r = 0.730",transform=right.transAxes,va="top")
fig.tight_layout()
fig.savefig(root / "shared-rotation.png",dpi=160)

"""Render the observed velocity bounds and conditional injection efficiencies."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
obs = pd.read_csv(HERE / "search/kinematics.csv", dtype={"source_id": str})
sim = pd.read_csv(HERE / "search/injection.csv")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"svg.fonttype":"none",
                     "svg.hashsalt":"gaia-dr3-supernova-survivors"})
fig, (a,b) = plt.subplots(1,2,figsize=(10,4.6),layout="constrained",
                         gridspec_kw={"width_ratios":[1.06,1]})
for axis in (a,b):
    axis.spines[["top","right"]].set_visible(False)
    axis.set_axisbelow(True)
for j,row in obs.iterrows():
    col = "#265a87" if pd.notna(row["name"]) else "#787674"
    a.plot([row.vt_lower95,row.corrected_vt],[j,j],color=col,lw=2)
    a.scatter(row.corrected_vt,j,s=21,color=col)
    a.plot(row.vt_lower95,j,marker="|",markersize=8,color=col)
labels=[str(n) if pd.notna(n) else "…"+s[-6:] for n,s in zip(obs["name"],obs.source_id)]
a.set_yticks(range(len(obs)),labels)
a.invert_yaxis()
a.axvline(1000,color="#b55725",ls="--",lw=1)
a.grid(axis="x",alpha=.17)
a.set_title("Observed Gaia sources",loc="left",fontweight="bold")
a.set_xlabel("Heliocentric transverse speed (km/s)")
a.text(0,1.067,"Tick: 95% lower bound   •   Dot: corrected estimate",
       transform=a.transAxes,fontsize=8,color="#42474b")
colours=["#173f58","#007e85","#c77321","#74688b"]
chosen=sim[sim.true_3d_kms.eq(1500)]
for (mag,s),colour in zip(chosen.groupby("g_lo"),colours):
    b.errorbar(s.distance_kpc,100*s.strong_fraction,yerr=100*s.strong_template_se,
               label=f"{mag} ≤ G < {s.g_hi.iloc[0]:g}",color=colour,marker="o",ms=4,capsize=2)
b.set_title("Recovery of simulated survivors",loc="left",fontweight="bold")
b.text(0,1.067,"1,500 km/s in 3D, isotropic directions",transform=b.transAxes,fontsize=8,color="#42474b")
b.set_ylim(-2,70);b.set_xscale("log",base=2);b.set_xticks([1,2,4,8],["1","2","4","8"])
b.set_xlabel("Input distance (kpc)");b.set_ylabel("Recovery after all cuts (%)")
b.grid(axis="y",alpha=.2);b.legend(frameon=False,fontsize=8)
fig.savefig(HERE/"sensitivity.svg",metadata={"Date":None},bbox_inches="tight")

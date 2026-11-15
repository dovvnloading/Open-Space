"""Draw the held-out ESPRESSO velocity comparison and temporal residuals."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from replay import EPOCH, design

ROOT = Path(__file__).resolve().parent
with (ROOT/"espresso.csv").open() as handle:
    rows = list(csv.DictReader(handle))
bjd = np.array([float(row["bjd"]) for row in rows])
v = np.array([float(row["rv"])*1000 for row in rows])
dates = [datetime.fromtimestamp((jd-2440587.5)*86400, tz=timezone.utc) for jd in bjd]
models = json.loads((ROOT/"harps_predictions.json").read_text())
model = min(models, key=lambda m:m["rms"])
predicted = design(bjd-EPOCH, model["P"]) @ model["params"][:8]
nuisance = np.column_stack([np.ones(len(bjd)), (bjd-bjd.mean())/1000])
correction = nuisance @ np.linalg.lstsq(nuisance, v-predicted, rcond=None)[0]
measured = v-correction
residual = measured-predicted
years = np.floor((bjd-2460175)/365.25).astype(int)
colors = ["#17608A", "#AA5422", "#386C55"]
symbols = ["o", "s", "D"]

plt.rcParams.update({"font.size":10, "axes.spines.top":False,
                     "axes.spines.right":False, "axes.edgecolor":"#84919A",
                     "axes.titlelocation":"left","savefig.facecolor":"white"})
fig = plt.figure(figsize=(9.0, 4.6), layout="constrained")
grid = fig.add_gridspec(1,2,width_ratios=[1,1.13])
ax, bx = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
for i in range(3):
    use = years==i
    ax.scatter(predicted[use], measured[use], color=colors[i], marker=symbols[i],
               s=39, edgecolor="white", linewidth=.6, zorder=3)
    bx.scatter(np.array(dates)[use], residual[use], color=colors[i], marker=symbols[i],
               s=39, edgecolor="white", linewidth=.6,zorder=3)
span = [-12, 12]
ax.plot(span,span,"--",color="#7E8B94",lw=1, zorder=1)
ax.set(xlim=span,ylim=span,xlabel="HARPS prediction (m/s)",ylabel="ESPRESSO velocity (m/s)",
       title="Predicted and measured velocity",aspect="equal")
ax.set_xticks([-10,-5,0,5,10]);ax.set_yticks([-10,-5,0,5,10])
ax.text(.03,.96,"Dashed line = exact agreement",transform=ax.transAxes,va="top",
        color="#4C5A63",fontsize=8)
bx.axhline(0,color="#7E8B94",ls="--",lw=1,zorder=1)
bx.set(ylabel="Measured minus predicted (m/s)",xlabel="ESPRESSO observation",
       title="Residuals over two years",ylim=(-10,10))
bx.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
bx.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
for panel in [ax,bx]:
    panel.grid(color="#DCE3E7",alpha=.8,linewidth=.6,zorder=0)
fig.suptitle("HD 32564  |  Independent test of the HARPS waveform",
             x=.015, ha="left",fontsize=15,fontweight="bold",color="#233642")
filename = Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/"espresso_validation.png"
fig.savefig(filename,dpi=155,metadata={"Software":"Open-Space"})

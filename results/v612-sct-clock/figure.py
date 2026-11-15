"""Plot the independently measured phase and pixel localization."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from replay import DATA, EPOCH, ROOT, main

def draw():
    report, rows, t, residual, f, coef = main()
    epoch=report["ztf"]["max_brightness_bjd_near_epoch"]
    sectors=json.loads((DATA/"tess-manifest.json").read_text())
    fig,ax=plt.subplots(2,2,figsize=(9.3,6.8),layout="constrained")
    fig.suptitle("V612 Scuti  ·  a 3.8295-hour photometric modulation",fontsize=15)
    for axis,sector,color in zip(ax.T,sectors,["#136f7b","#a1472e"]):
        s=sector["sector"];d=pd.read_csv(DATA/f"tess-{s}.csv.gz")
        phase=((d.btjd+2457000-epoch)*f+.5)%1-.5
        categories=pd.cut(phase,np.linspace(-.5,.5,19),labels=False)
        series=d.corrected_e_per_s.groupby(categories)
        x=(np.arange(18)+.5)/18-.5
        axis[0].plot(x,series.mean(),marker="o",color=color,label="TESS")
        amp=report["tess"][s]["semiamplitude"]
        tx=np.linspace(-.5,.5,501)
        axis[0].plot(tx,amp*np.cos(2*np.pi*tx),c="0.23",ls="--",label="ZTF phase, TESS amplitude")
        axis[0].set(xlabel="Cycles from ZTF brightness maximum",ylabel="Flux residual (electrons/s)",title=f"TESS sector {s}")
        axis[0].legend(frameon=False,fontsize=8);axis[0].axhline(0,color="0.8",lw=.5)
        im=axis[1].imshow(sector["pixel_semiamplitude_e_per_s"],origin="lower",cmap="magma",interpolation="nearest")
        axis[1].scatter(*sector["target_xy"],facecolors="none",edgecolors="cyan",s=90,linewidth=1.4,label="Gaia position")
        axis[1].set(xlabel="Cutout pixel x",ylabel="Cutout pixel y",xticks=range(5),yticks=range(5));axis[1].legend(frameon=False,labelcolor="white",fontsize=8,loc="upper left")
        fig.colorbar(im,ax=axis[1],label="Semiamplitude (electrons/s)",shrink=.8)
    fig.savefig(ROOT/"phase-and-pixels.png",dpi=155)
if __name__=="__main__":
    draw()

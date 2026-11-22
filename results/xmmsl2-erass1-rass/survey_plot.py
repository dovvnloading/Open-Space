"""Compare the measured soft fluxes and show where a map sets the faint epoch."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def make_plot(out):
    cat=pd.read_csv(out/"catalog_sources.csv.gz",low_memory=False).set_index("detection_uid")
    pairs=pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    fig,axes=plt.subplots(1,2,figsize=(9.5,4.8),sharex=True,sharey=True,layout="constrained")
    ticks=[-15,-14,-13,-12,-11,-10,-9,-8]
    collections=[]
    for ax,m,label in zip(axes,["ero","rass"],["eROSITA eRASS1","ROSAT RASS"]):
        valid=(cat.xmm_quality & cat[m+"_quality"] & ~cat[m+"_ambiguous"]
               & cat.xmm_flux.gt(0) & cat[m+"_flux"].gt(0))
        x=np.log10(cat.loc[valid,"xmm_flux"])
        y=np.log10(cat.loc[valid,m+"_flux"])
        collections.append(ax.hexbin(x,y,gridsize=50,mincnt=1,extent=(-15,-8,-15,-8),
                                     linewidths=0,cmap="Blues",norm=LogNorm()))
        p=pairs.loc[pairs.photometric_lead & pairs.high.isin([m,"xmm"]) & pairs.low.isin([m,"xmm"])]
        det=p.loc[p.low_kind.eq("detection")].drop_duplicates("detection_uid")
        if len(det):
            c=cat.loc[det.detection_uid]
            ax.scatter(np.log10(c.xmm_flux),np.log10(c[m+"_flux"]),s=15,
                       facecolors="none",edgecolors="#b3421d",linewidths=.9,zorder=3)
        ul=p.loc[p.low_kind.eq("upper_limit") & p.low.eq(m)]
        if len(ul):
            ax.scatter(np.log10(ul.high_flux),np.log10(ul.low_flux_or_limit),s=18,
                       marker="v",color="#8050a3",alpha=.85,zorder=4)
        for offset,style in [(0,"-"),(-1,"--"),(1,"--")]:
            ax.plot(ticks,np.array(ticks)+offset,ls=style,color="#707070",lw=.8,zorder=1)
        ax.set(title=label,xlim=(-15,-8),ylim=(-15,-8),xticks=ticks[::2],yticks=ticks[::2],
               xlabel="XMM slew flux")
        ax.set_aspect("equal",adjustable="box")
        ax.set_xticklabels([rf"$10^{{{v}}}$" for v in ticks[::2]])
        ax.set_yticklabels([rf"$10^{{{v}}}$" for v in ticks[::2]])
    axes[0].set_ylabel("Comparison survey flux or upper limit")
    maximum=max(int(c.get_array().max()) for c in collections)
    for c in collections:
        c.set_clim(1,max(2,maximum))
    fig.colorbar(collections[0],ax=axes,shrink=.73,label="Clean slew detections per bin")
    handles=[
        Line2D([],[],ls="",marker="o",mec="#b3421d",mfc="none",label="Two detections passing the tenfold cut"),
        Line2D([],[],ls="",marker="v",color="#8050a3",label="Map upper limit passing the tenfold cut")]
    fig.legend(handles=handles,loc="outside lower center",ncol=2,frameon=False,fontsize=9)
    fig.suptitle("Soft X-ray brightness across survey epochs",fontsize=15,weight="bold")
    fig.savefig(out/"survey-comparison.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    make_plot(p.parse_args().out)

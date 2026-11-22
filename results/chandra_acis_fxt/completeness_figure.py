"""Plot the measured primary injection recovery grid from the final survey report."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COUNTS=(8, 16, 32, 64)
DURATIONS=(60, 180, 600, 1800, 3600)


def table(path):
    recovered=np.zeros((len(COUNTS),len(DURATIONS)),dtype=np.int64)
    attempted=np.zeros_like(recovered)
    seen=set()
    with path.open(encoding="utf-8",newline="") as fh:
        for row in csv.DictReader(fh):
            n,t=int(row["input_counts"]),int(row["duration_s"])
            ident=n,t,float(row["psf_sigma_sky_px"])
            if n not in COUNTS or t not in DURATIONS or ident in seen:
                raise ValueError(f"Unexpected injection grid cell {ident}")
            seen.add(ident)
            i,j=COUNTS.index(n),DURATIONS.index(t)
            k,m,planned=(int(row[x]) for x in ("recovered","attempted","scheduled"))
            if not 0<=k<=m<=planned:
                raise ValueError("Invalid injection recovery counts")
            recovered[i,j]+=k
            attempted[i,j]+=m
    return recovered,attempted


def render(source, output, expected_attempted, expected_recovered):
    good,total=table(source)
    if good.sum()!=expected_recovered or total.sum()!=expected_attempted:
        raise ValueError("The plotted injections do not match the survey totals")
    fraction=np.divide(good,total,out=np.full(good.shape,np.nan),where=total>0)
    colors=matplotlib.colormaps["viridis"].copy()
    colors.set_bad("#eeeeee")
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":10,"savefig.facecolor":"white"}):
        fig,ax=plt.subplots(figsize=(9,4.6),layout="constrained")
        ax.imshow(np.ma.masked_invalid(fraction),vmin=0,vmax=1,cmap=colors,aspect="auto",origin="lower")
        ax.set_title("Recovery of synthetic flares in the dated Chandra fields",loc="left",fontweight="bold",pad=14)
        ax.set_xticks(range(len(DURATIONS)),[str(t//60) for t in DURATIONS])
        ax.set_yticks(range(len(COUNTS)),[str(n) for n in COUNTS])
        ax.set_xlabel("Injected flare duration (minutes)",labelpad=8)
        ax.set_ylabel("Injected source photons",labelpad=8)
        ax.tick_params(length=0,pad=7)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for i in range(len(COUNTS)):
            for j in range(len(DURATIONS)):
                p=fraction[i,j]
                if np.isnan(p):
                    label,color="No attempts","#374151"
                else:
                    label=f"{p:.0%}\n{good[i,j]:,} / {total[i,j]:,}"
                    rgb=np.asarray(colors(p)[:3])
                    linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
                    color="black" if linear @ np.array([.2126,.7152,.0722])>=.179 else "white"
                ax.text(j,i,label,ha="center",va="center",color=color,fontsize=10)
        output=Path(output)
        temp=output.with_suffix(".tmp")
        fig.savefig(temp,format="png",dpi=145,metadata={"Software":"Open-Space Chandra ACIS survey"})
        plt.close(fig)
        temp.replace(output)

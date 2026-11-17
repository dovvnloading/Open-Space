"""Plot the published IRAS endpoint and the July-compatible moving-track predictions."""
import csv
from pathlib import Path
import re
import tarfile

from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scan_tracks import ROOT, filtered


def main():
    directory=ROOT/"scanpi"
    with (directory/"track_measurements.csv").open() as file:
        metadata={int(row["scan"]):row for row in csv.DictReader(file) if row["band_um"]=="60"}
    table=Table.read(directory/"b3_summary.tbl",format="ascii.ipac")
    fig,axes=plt.subplots(1,2,figsize=(10,3.65),sharey=True,layout="constrained")
    with tarfile.open(directory/"profiles.tar.gz") as profiles:
        for ax,ids,title in zip(axes,[range(4,10),range(10,13)],["July 9–10, 1983","August 28 and September 1, 1983"]):
            responses=[]
            weights=[]
            bounds=[]
            for i in ids:
                row=table[table["scan"]==i][0]
                path="b3/"+Path(re.search(r'href="([^"]+)"',row["table"]).group(1)).name
                raw=profiles.extractfile(path).read().decode().splitlines()
                x,f=np.array([[float(v) for v in s.split()[:2]] for s in raw if s.startswith(" ")]).T
                signal=filtered(x,f,1.5)
                meta=metadata[i]
                responses.append(signal)
                weights.append(1/max(float(meta["control_robust_scatter_Jy"]),.04)**2)
                bounds.extend(float(meta[k]) for k in ("july_gated_along_min_arcmin","july_gated_along_max_arcmin"))
                ax.plot(x,signal,color="#8b949b",lw=.6,alpha=.6)
            ax.axvspan(min(bounds),max(bounds),color="#d79722",alpha=.3,label="July-compatible moving positions")
            ax.axvline(0,color="#233d56",ls="--",lw=1,label="Published IRAS position")
            ax.plot(x,np.average(responses,axis=0,weights=weights),color="#006f82",lw=2.1,label="Weighted detector mean")
            ax.axhline(0,color="#252525",lw=.5,alpha=.55)
            ax.set_title(title,fontsize=12)
            ax.set_xlabel("Along-scan offset from IRAS position (arcmin)",fontsize=9)
            ax.set_xlim(-8,5)
            ax.set_xticks([-8,-6,-4,-2,0,2,4])
            ax.set_ylim(-.32,.5)
            ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("60 µm local template amplitude (Jy)")
    axes[0].legend(loc="lower left",frameon=False,fontsize=7.6)
    fig.savefig(directory/"scan_profiles.png",dpi=160)
    plt.close(fig)


if __name__=="__main__":main()

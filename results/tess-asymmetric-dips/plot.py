"""Plot measured Sector 105 eclipse profiles and injection recovery."""
import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    data = json.loads(gzip.decompress((here/"confirmation/event-profiles.json.gz").read_bytes()))
    event = next(r for r in data if int(r["tic"])==316937197)
    t = np.array([b["btjd"] for b in event["bins"]])
    hours = 24*(t-event["peak_btjd"])
    def series(k):
        y=np.array([b[k] if b[k] is not None else np.nan for b in event["bins"]])
        return 100*(y/np.nanmedian(y[(abs(hours)>12)&(abs(hours)<36)])-1)
    rate=json.loads((here/"confirmation/report.json").read_text())["limits"]
    heat=np.array([100*r["recovered"]/r["trials"] for r in rate]).reshape(3,3)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","svg.hashsalt":"tess105"})
    fig=plt.figure(figsize=(11,4.4),layout="constrained")
    gs=fig.add_gridspec(2,2,width_ratios=(1.6,1))
    ax=fig.add_subplot(gs[0,0])
    bx=fig.add_subplot(gs[1,0],sharex=ax)
    for panel,k,c,lab in [(ax,"DET_FLUX","#22537b","QLP corrected flux"),(bx,"pix3","#975126","3 × 3 pixel aperture")]:
        use=abs(hours)<15
        panel.plot(hours[use],series(k)[use],color=c,marker=".",lw=1.1,ms=3)
        panel.axhline(0,color="#777777",lw=.5)
        panel.set_ylabel("Flux change (%)")
        panel.text(.02,.08,lab,transform=panel.transAxes,color=c,fontsize=9)
    ax.tick_params(labelbottom=False)
    ax.set_title("TIC 316937197 • Sector 105",loc="left",weight="bold")
    bx.set_xlabel("Hours from the QLP template peak (BTJD 4217.632)")
    cx=fig.add_subplot(gs[:,1])
    im=cx.imshow(heat,vmin=0,vmax=75,cmap="Blues",aspect="auto")
    cx.set_xticks(range(3),["6 h","12 h","24 h"])
    cx.set_yticks(range(3),["1%","2%","5%"])
    cx.set_xlabel("Injected duration")
    cx.set_ylabel("Injected peak depth")
    cx.set_title("Recovery in 256 held-out stars",loc="left",weight="bold")
    for i in range(3):
        for j in range(3):
            cx.text(j,i,f"{heat[i,j]:.1f}%",ha="center",va="center",
                    color="white" if heat[i,j]>44 else "#183044",fontsize=13,weight="bold")
    fig.savefig(args.out,metadata={"Date":None})
    plt.close(fig)


if __name__=="__main__":
    main()

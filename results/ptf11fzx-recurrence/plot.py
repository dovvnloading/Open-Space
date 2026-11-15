"""Plot the g-band outbursts from frozen direct-image photometry (matplotlib)."""
from collections import defaultdict
from statistics import median
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze import ROOT, read


def main():
    direct = [r for r in read("ztf_direct.csv") if r["filtercode"] == "zg"
              and int(r["catflags"]) == 0 and float(r["magerr"]) <= 0.25]
    ptf = [r for r in read("ptf.csv") if r["fid"] == "1" and r["goodflag"] == "1" and r["qaflags"] == "0"]
    data = [(float(r["mjd"]), float(r["mag"]), float(r["magerr"])) for r in direct]
    data += [(float(r["obsmjd"]), float(r["mag_autocorr"]), float(r["magerr_auto"])) for r in ptf]
    nights = defaultdict(list)
    for t,m,e in data:
        nights[int(t)].append((t,m,e))
    points = [tuple(median(c) for c in zip(*rows)) for _, rows in sorted(nights.items())]
    fig, (overview, detail) = plt.subplots(1,2,figsize=(10,4.9),sharey=True,gridspec_kw={"width_ratios":[1.05,1]})
    series = [(56949,56956,"2014 · PTF","#8d52a1","s"),(58426,58435,"2018 · ZTF","#33658a","o"),
              (59182,59219,"2020 · ZTF","#bf4a20","D"),(60141,60149,"2023 · ZTF","#26816f","^"),
              (60511,60518,"2024 · ZTF","#766338","v")]
    for start,end,label,color,marker in series:
        group = [p for p in points if start <= p[0] < end]
        ts,ms,errors=zip(*group)
        t0=ts[0]
        for ax in [overview,detail]:
            if ax is overview and start not in (58426,59182):
                continue
            ax.errorbar([t-t0 for t in ts],ms,yerr=errors,fmt=marker,ms=4.5,elinewidth=.8,
                        capsize=2,alpha=.93,color=color,label=label)
    for ax,xend,title in [(overview,37,"The 2018 and 2020 declines"),(detail,8.6,"First eight days across outbursts")]:
        ax.set_xlim(-.8,xend)
        ax.set_xlabel("Days since the first plotted night")
        ax.set_title(title,loc="left",fontsize=11,pad=10,weight="bold")
        ax.tick_params(length=3)
        ax.grid(axis="y",color="#ddd",lw=.6)
        for side in ["right","top"]:ax.spines[side].set_visible(False)
        ax.legend(loc="lower left",frameon=False,fontsize=8,ncols=2)
    overview.set_ylabel("g magnitude (brighter upward)")
    overview.set_ylim(21.85,18.05)
    fig.suptitle("Optical outbursts of PTF11fzx",x=.075,y=.96,ha="left",fontsize=17,weight="bold")
    fig.subplots_adjust(left=.075,right=.99,top=.79,bottom=.16,wspace=.14)
    fig.savefig(ROOT / "outbursts.png",dpi=160,facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()

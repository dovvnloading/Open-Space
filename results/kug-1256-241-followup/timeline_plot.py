"""Plot deposited epochs with separate magnitude systems, apertures and X-ray products."""
import csv
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as dates
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "outputs"
with (OUT / "timeline.csv").open() as f:
    rows = list(csv.DictReader(f))


def date(text):
    return datetime.fromisoformat(text) if text else None


def draw():
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True,
                             gridspec_kw={"height_ratios":[1.35,1,.95,.75]})
    ax, opt, ir, rug = axes
    ax.set_title("KUG 1256+241", loc="left", weight="bold", fontsize=16)
    ax.set_ylabel("Observed 0.2–2 keV\n(erg cm⁻² s⁻¹)")
    ax.set_yscale("log")
    all_labels=set()
    colors = {"ROSAT":"#8041a1","eROSITA":"#c14542","XMM":"#b26b1a","Chandra":"#127c78","Swift":"#245c9c"}
    for row in rows:
        start, end = date(row["utc"]), date(row["end_utc"])
        if not start: continue
        a, b, unit = row["archive"], row["band"], row["unit"]
        if b in ("0.2–2.0 keV","0.2-2.0 keV") and row["value"]:
            group = next((s for s in colors if s.lower() in a.lower()), "Other")
            c=colors.get(group,"#444")
            y = float(row["value"])
            x = start+(end-start)/2 if end and (end-start).days < 366 else start
            if end and (end-start).days > 365:
                # Catalogue averages over several years can conceal the visit-to-visit evolution.
                continue
            label = group if group not in all_labels else None
            all_labels.add(group)
            if end and end > start:
                ax.hlines(y,start,end,color=c,lw=1.2)
            ax.scatter(x,y,marker="v" if "limit" in row["kind"].lower() else "o",
                       s=25,color=c,alpha=.88,label=label,zorder=3)
            if row["error"] and float(row["error"]) < y:
                ax.errorbar(x,y,yerr=float(row["error"]),fmt="none",ecolor=c,lw=.9)
        elif unit.startswith("mag") and row["value"] and row["kind"]=="detection":
            if a in ("AllWISE","NEOWISE"):
                group="WISE "+b; panel=ir
            elif a in ("ZTF","CRTS","Pan-STARRS1 DR2") and b in ("g","r","zg","zr","CV"):
                group=a.split()[0]+" "+b; panel=opt
            else:
                continue
            palette={"WISE W1":"#146b61","WISE W2":"#bd6335","CRTS CV":"#777777",
                     "Pan-STARRS1 g":"#73782e","Pan-STARRS1 r":"#925a4c","ZTF zg":"#358641","ZTF zr":"#ba594a",
                     "ZTF g":"#358641","ZTF r":"#ba594a"}
            panel.plot(start,float(row["value"]),".",ms=2,alpha=.6,color=palette.get(group,"#556382"),
                       label=group if group not in all_labels else None)
            all_labels.add(group)
        if row["kind"]=="spectrum":
            rug.plot(start,2,"|",ms=14,color="#314d68")
        elif ("GHz" in b or "VLASS" in a or "FIRST" in a) and row["kind"] in ("detection","upper limit"):
            rug.plot(start,0,"|",ms=13,color="#945e32")
        elif any(s in a.lower() for s in ("uvot","galex")) and row["kind"]=="detection":
            rug.plot(start,1,"|",ms=13,color="#675a8d")
    for panel,title in ((opt,"Optical magnitude"),(ir,"WISE magnitude\n(Vega)")):
        panel.invert_yaxis()
        panel.set_ylabel(title)
    rug.set_yticks([0,1,2],["Radio","Ultraviolet","Spectrum"])
    rug.set_ylim(-.6,2.65)
    rug.set_xlabel("Observing year (UTC)")
    for panel in axes:
        panel.grid(axis="x",alpha=.18)
        panel.spines[["right","top"]].set_visible(False)
        handles, labels=panel.get_legend_handles_labels()
        if labels:
            panel.legend(ncol=3,loc="best",fontsize=7,frameon=False)
    rug.xaxis.set_major_locator(dates.YearLocator(5))
    rug.xaxis.set_major_formatter(dates.DateFormatter("%Y"))
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(OUT / "timeline.png",dpi=155,bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    draw()

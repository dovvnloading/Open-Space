"""Recompute plotted completeness limits directly from deposited injection trials."""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
NIRCAM = BASE / "nircam-hd53143-hd38206/outputs"
DISKS = BASE / "hd21997-hd202628/measurements"


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def curve(path, ratio_key, minimum, selection=None, sustained=False):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows(path):
        if selection and not selection(row):
            continue
        grouped[float(row["radius_arcsec"])][float(row[ratio_key])].append(int(row["recovered"]))
    output = []
    for radius, levels in sorted(grouped.items()):
        fractions = {c: sum(z)/len(z) for c,z in levels.items()}
        for c in sorted(levels):
            fraction = min(fractions[v] for v in levels if v >= c) if sustained else fractions[c]
            if fraction >= minimum:
                output.append((radius,c,sum(levels[c]),len(levels[c])))
                break
        else:
            output.append((radius,None,0,0))
    return output


def main():
    fig, axs = plt.subplots(2,3,figsize=(11,6.5),layout="constrained")
    colors = ["#006D77", "#BF5B17", "#60509A"]
    summary = []
    def draw(ax, field, band, label, data, completeness, color):
        valid = [v for v in data if v[1] is not None]
        ax.plot([v[0] for v in valid],[v[1] for v in valid], marker="o",lw=1.4,
                markersize=4,color=color,label=label)
        for radius,ratio,n,total in data:
            summary.append(dict(field=field,band=band,series=label,radius_arcsec=radius,
                                recovery_fraction_required=completeness,limiting_tested_flux_ratio=ratio,
                                recovered_at_limit=n if ratio else None,tried_at_limit=total if ratio else None))
    ax=axs[0,0]
    for color,band in zip(colors,["F164N","F277W","F323N"]):
        data=curve(BASE/"lawd66/joint_injections.csv","ratio",.90,lambda r:r["band"]==band)
        levels=rows(BASE/"lawd66/joint_contrast.csv")
        for rad,ratio,n,total in data:
            if ratio:
                match=next(r for r in levels if r["band"]==band and float(r["radius_arcsec"])==rad and float(r["ratio"])==ratio)
                assert (int(match["recovered"]),int(match["trials"]))==(n,total)
        draw(ax,"LAWD 66",band,band,data,.90,color)
    ax.set_title("LAWD 66   ·   paired epochs",loc="left")
    groups=[("HD 53143",axs[0,1],[(NIRCAM/"hd53143_o001","Roll 001"),(NIRCAM/"hd53143_o008","Roll 008")]),
            ("HD 38206",axs[0,2],[(NIRCAM/"hd38206_o141","Jan 2025"),(NIRCAM/"hd38206_o150","Oct 2025")]),
            ("HD 21997",axs[1,0],[(DISKS/"HD21997","Feb 2024")]),
            ("HD 202628",axs[1,1],[(DISKS/"HD202628","Jun 2024")])]
    for field,ax,series in groups:
        for color,(prefix,label) in zip(colors,series):
            data=curve(prefix.with_name(prefix.name+"_injections.csv"),"contrast",.95,sustained=True)
            expected=rows(prefix.with_name(prefix.name+"_contrast.csv"))
            assert len(data)==len(expected)
            for (radius,ratio,_,_),e in zip(data,expected):
                assert radius==float(e["radius_arcsec"])
                assert (ratio is None and not e["contrast95"]) or math.isclose(ratio,float(e["contrast95"]),rel_tol=1e-12)
            draw(ax,field,"F444W",label,data,.95,color)
        ax.set_title(field+"   ·   F444W",loc="left")
    ax=axs[1,2]
    for color,(obs,label) in zip(colors,[("o003","Roll 003"),("o004","Roll 004"),("c1005","Combined")]):
        folder=BASE/"hd53143-miri/measurements"/("jw08142-"+obs)
        data=curve(folder/"injection_trials.csv","contrast",.95)
        expected=rows(folder/"contrast_recovery.csv")
        for rad,ratio,n,total in data:
            match=next(r for r in expected if float(r["radius_arcsec"])==rad and float(r["contrast"])==ratio)
            assert (int(match["recovered"]),int(match["tried"]))==(n,total)
        draw(ax,"HD 53143","F2550W",label,data,.95,color)
    ax.set_title("HD 53143   ·   F2550W",loc="left")
    for ax in axs.flat:
        ax.set_yscale("log")
        ax.grid(axis="y",alpha=.17,which="both")
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(frameon=False,fontsize=8)
        ax.set_xlabel("Separation (arcsec)")
    for ax in axs[:,0]: ax.set_ylabel("Companion / host flux")
    fig.savefig(BASE/"contrast_curves.png",dpi=150,metadata={"Software":"Open-Space"})
    with (BASE/"contrast_summary.csv").open("w") as stream:
        writer=csv.DictWriter(stream,list(summary[0]));writer.writeheader();writer.writerows(summary)


if __name__=="__main__":
    main()

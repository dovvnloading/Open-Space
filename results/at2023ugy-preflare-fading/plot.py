"""Draw the direct NEOWISE and Gaia measurements and same-field comparison."""
import argparse
import csv
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from measure import EPOCH

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    with (args.out / "visits.csv").open() as f:
        visits = list(csv.DictReader(f))
    with (args.out / "field.csv").open() as f:
        field = list(csv.DictReader(f))
    with (args.cache / "gaia.csv").open() as f:
        nights = defaultdict(list)
        for row in list(csv.reader(f))[2:]:
            try:
                day = int(float(row[1]) - 2400000.5)
                mag = float(row[2])
                nights[day].append(mag)
            except (IndexError, ValueError):
                continue
    colors = ["#0c7190", "#ba5727", "#704F88"]
    plt.rcParams.update({"font.size": 10.5, "svg.fonttype": "none", "svg.hashsalt": "at2023ugy",
                         "axes.spines.right": False, "axes.spines.top": False,
                         "axes.labelcolor": "#2b353c", "text.color": "#25333a"})
    fig = plt.figure(figsize=(10.8, 5.6))
    layout = fig.add_gridspec(2, 2, width_ratios=[1.85, 1], hspace=.16, wspace=.34)
    ax = fig.add_subplot(layout[0, 0])
    gx = fig.add_subplot(layout[1, 0], sharex=ax)
    fx = fig.add_subplot(layout[:, 1])
    for k, band in enumerate(["w1", "w2"]):
        series = [(EPOCH + timedelta(days=float(v["mjd"])), float(v[band+"mpro"])) for v in visits]
        ax.plot(*zip(*series), "o-", color=colors[k], label=band.upper(), markersize=4, lw=1.15)
    gx.scatter([EPOCH + timedelta(days=t) for t in sorted(nights)],
               [median(nights[t]) for t in sorted(nights)], s=14, color=colors[2])
    for axis, band in [(ax,"NEOWISE magnitude"), (gx, "Gaia G magnitude")]:
        axis.set_ylabel(band)
        axis.invert_yaxis()
        axis.grid(axis="y", color="#e3e6e7", linewidth=.6)
        for start, end in [(2014,2016), (2021,2023)]:
            axis.axvspan(date(start,1,1), date(end,1,1), color="#7aaa9d", lw=0, alpha=.16, zorder=-1)
        axis.axvline(date(2023,9,22), color="#383d40", ls="--", lw=1)
    ax.tick_params(labelbottom=False)
    ax.legend(frameon=False, loc="upper left", ncol=2)
    gx.xaxis.set_major_locator(mdates.YearLocator(2))
    gx.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    gx.set_xlabel("Calendar year")
    ax.set_title("Brightness at the transient position", loc="left", fontsize=13, fontweight="bold", pad=10)
    gx.annotate("2023 optical detection", xy=(date(2023,9,22),.57), xycoords=("data","axes fraction"),
                xytext=(-8,8), textcoords="offset points", ha="right", fontsize=9)
    ax.set_xlim(date(2013,10,1),date(2024,12,31))
    fx.axvline(0,lw=.8,color="#c8cdcf")
    fx.axhline(0,lw=.8,color="#c8cdcf")
    for r in field:
        target = r["target"] == "True"
        fx.scatter(float(r["w1_fading"]),float(r["w2_fading"]),color=colors[0] if target else "#777e83",
                   s=85 if target else 33,marker="D" if target else "o",zorder=2)
    fx.annotate("AT 2023ugy host", (.299,.3945), xytext=(-6,-19), textcoords="offset points", ha="right",
                color=colors[0], fontsize=10, fontweight="bold")
    fx.text(.14,-.10,"8 field controls",ha="center",color="#626970",fontsize=10)
    fx.set(xlabel="W1 fading (mag)",ylabel="W2 fading (mag)",xlim=(-.065,.35),ylim=(-.21,.46))
    fx.xaxis.set_major_locator(MultipleLocator(.1))
    fx.yaxis.set_major_locator(MultipleLocator(.1))
    fx.set_title("2014–15 to 2021–22", loc="left",fontsize=13,fontweight="bold",pad=10)
    fig.subplots_adjust(left=.09, bottom=.18, right=.97, top=.87)
    fig.text(.09,.045,"Shading marks the comparison windows. Each infrared point is a visit median. Each Gaia point is a daily median.",
             fontsize=9,color="#555f64")
    fig.savefig(args.out / "fading.svg", format="svg", metadata={"Date":None})
    plt.close(fig)

if __name__ == "__main__":
    main()

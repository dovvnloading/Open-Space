"""Plot the sky location of every published branch-b candidate and Ret II target."""
import argparse
from math import cos, radians
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

from reproduce import get_sources


def make_figure(cache, output):
    get_sources(cache)
    atlas = []
    for row in (cache/"ibata.txt").read_text().splitlines():
        if row[:19].strip().isdigit() and row[95:97].strip() == "14":
            atlas.append((float(row[20:30]), float(row[31:41]), int(row[92:94]) > 0))
    ji = []
    for row in (cache/"ji.txt").read_text().splitlines():
        if not row[:3].strip().isdigit() or row[15:16] != "-":
            continue
        ra = 15*(int(row[4:6])+int(row[7:9])/60+float(row[10:14])/3600)
        dec = -(int(row[16:18])+int(row[19:21])/60+float(row[22:27])/3600)
        ji.append((ra, dec, row[71] == "1", row[73] == "1"))
    def points(ax, rows, **opts):
        ax.scatter([r[0] for r in rows], [r[1] for r in rows], **opts)

    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.labelcolor": "#26332f",
                         "xtick.color": "#34453f", "ytick.color": "#34453f"})
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
    fig.subplots_adjust(left=.075, right=.98, bottom=.22, top=.85, wspace=.3)
    line_color = "#cd462c"
    all_rv = [star for star in atlas if star[2]]
    points(a, [star for star in atlas if not star[2]], s=18, color="#688d9b", label="No atlas RV (58)")
    points(a, all_rv, s=26, color=line_color, marker="o", label="Measured atlas RV (12)")
    a.set(xlim=(59, 41.7), ylim=(-55.4, -48.6), xlabel="Right ascension (deg)", ylabel="Declination (deg)")
    a.set_aspect(1/cos(radians(-52)))
    a.set_title("NGC 1261b candidate catalog", loc="left", fontweight="bold")
    a.add_patch(Rectangle((53.55, -54.29), .72, .48, edgecolor="#3c4a43", facecolor="none",
                          linewidth=.9, linestyle="--"))
    a.annotate("Reticulum II", (53.9, -53.8), (51.2, -52.9), fontsize=9,
               arrowprops={"arrowstyle":"-", "color":"#55665a", "lw":.75})
    a.legend(loc="upper left", bbox_to_anchor=(0,-.22), frameon=False, ncol=2, fontsize=8,
             handletextpad=.2, columnspacing=.6, borderaxespad=0)
    points(b, [star for star in ji if star[2]], s=27, color="#167c6f", label="Ji firm member")
    points(b, [star for star in ji if star[3]], s=46, marker="^", facecolors="white",
           edgecolors="#167c6f", linewidth=.95, label="Ji candidate member")
    points(b, all_rv, s=89, facecolors="none", edgecolors=line_color, linewidth=1.15, label="Atlas RV star")
    b.set(xlim=(54.27, 53.55), ylim=(-54.29,-53.81), xlabel="Right ascension (deg)", ylabel="Declination (deg)")
    b.set_aspect(1/cos(radians(-54.05)))
    b.set_title("Reticulum II field", loc="left", fontweight="bold")
    b.xaxis.set_major_locator(MultipleLocator(.2))
    b.yaxis.set_major_locator(MultipleLocator(.1))
    b.legend(loc="upper left", bbox_to_anchor=(-.02,-.22), frameon=False, ncol=2, fontsize=8,
             handletextpad=.2, columnspacing=.6, borderaxespad=0)
    for ax in (a,b):
        ax.tick_params(length=3)
        ax.spines[["left","bottom"]].set_color("#8d9b99")
    fig.savefig(output, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    make_figure(args.cache, args.output)

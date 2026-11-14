"""Plot actual phase-binned OGLE-IV I photometry for the four eclipse-depth leads."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import itertools

import matplotlib.pyplot as plt
import numpy as np
from followup import phase
from matplotlib.lines import Line2D


def plot(raw, output):
    cat = {l.split()[0]: l.split() for l in (raw / "ecl.dat").read_text().splitlines()}
    fig, ax = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
    fig.patch.set_facecolor("white")
    colors = ["#006762", "#B34122"]
    for i, number in enumerate(["07730", "06703", "04349", "09564"]):
        axis = ax.ravel()[i]
        star = "OGLE-SMC-ECL-" + number
        per, ep = map(float, cat[star][3:5])
        t, m, e = np.loadtxt(raw / "phot/I" / (star + ".dat")).T
        p = phase(t, ep, per)
        good = (e > 0) & (e <= 0.1)
        for j, subset in enumerate([t < 2457000, t >= 2459700]):
            take = good & subset
            ref = np.median(m[take & (abs(p) > 0.16) & (abs(p) < 0.34)])
            ph = (p + 0.20) % 1 - 0.20
            bins = np.linspace(-0.20, 0.80, 51)
            x = (bins[1:] + bins[:-1]) / 2
            y = []
            for low, high in itertools.pairwise(bins):
                select = take & (ph >= low) & (ph < high)
                y.append(np.median(m[select]) - ref if select.sum() >= 3 else np.nan)
            axis.plot(x, y, marker=".", markersize=3, lw=1.5, color=colors[j], zorder=3)
        axis.axhline(0, lw=0.75, color="#9eacb2")
        axis.set(
            xlim=(-0.20, 0.8),
            ylim=(0.56, -0.065),
            xticks=[0, 0.25, 0.5, 0.75],
            yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5],
        )
        axis.set_title(
            star + "  ·  " + f"{per:.4f} d", loc="left", fontsize=12, weight="semibold"
        )
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#e2e8e8", linewidth=0.6)
        if i // 2 == 1:
            axis.set_xlabel("Orbital phase")
        if i % 2 == 0:
            axis.set_ylabel("Dimming in I band (mag)")
    fig.legend(
        [Line2D([0], [0], color=c, lw=2) for c in colors],
        ["2010–2014", "2022–2024"],
        loc="outside upper center",
        ncols=2,
        frameon=False,
    )
    fig.savefig(output, dpi=145)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    a = parser.parse_args()
    plot(a.raw, a.output)

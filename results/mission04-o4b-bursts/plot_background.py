"""Plot the observed distribution of whole-campaign maxima. No tail extrapolation."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, default=Path(__file__).parent / "results")
    args = p.parse_args()
    survey = json.loads((args.results / "survey.json").read_text())
    background = np.sort([row[-1] for row in survey["campaign_slides"] if row[0]])
    loudest = survey["foreground_loudest"]
    n = len(background)
    # Q(x) = count(maximum >= x) / n, with exact inclusion of ties at x.
    x = np.r_[2., np.nextafter(background, np.inf)]
    fraction = np.arange(n, -1, -1) / n
    fig, ax = plt.subplots(figsize=(8, 4.7), layout="constrained")
    ax.step(np.r_[x[:-1], background[-1]], np.r_[fraction[:-1], 1/n],
            where="post", color="#244e6b", linewidth=2)
    ax.set(xscale="log", yscale="log", xlim=(2.25, 280), ylim=(.005, 1.24),
           xlabel="Maximum coherent score across the ten survey intervals",
           ylabel="Fraction of shifted campaigns reaching score")
    ax.set_title("O4b burst search  |  128 time-shifted campaigns", loc="left", weight="bold")
    ax.set_xticks([3, 5, 10, 20, 50, 100, 200], labels=["3", "5", "10", "20", "50", "100", "200"])
    ax.set_yticks([1, .5, .2, .1, .05, .02, .01], labels=["100%", "50%", "20%", "10%", "5%", "2%", "1%"])
    ax.grid(axis="y", color="#d6dadd", linewidth=.6)
    ax.spines[["right", "top"]].set_visible(False)
    y = np.sum(background >= loudest) / n
    ax.scatter(loudest, y, color="#b44621", s=45, zorder=5)
    ax.annotate(f"Observed maximum {loudest:.2f}\n{round(y*n)}/{n} shifted campaigns reached it",
                xy=(loudest, y), xytext=(6.0, .49), color="#883512",
                arrowprops={"arrowstyle":"-", "color":"#b44621"}, fontsize=10)
    ax.axvline(20, color="#6a6d6b", linestyle=":", linewidth=1.3)
    ax.text(22, .10, "Fixed rate-limit\nthreshold 20", color="#555956", fontsize=10)
    ax.scatter(background[-1], 1/n, color="#244e6b", s=30, zorder=4)
    ax.annotate(f"Loudest time shift {background[-1]:.0f}", (background[-1], 1/n),
                xytext=(145, .016), ha="right", fontsize=9, color="#244e6b",
                arrowprops={"arrowstyle":"-", "color":"#244e6b"})
    control = next(r for r in (json.loads(f.read_text()) for f in (args.results / "epochs").glob("*.json"))
                   if r["role"] == "control")
    ctrl = control["follow_up"][0]["score"]
    ax.axvline(ctrl, color="#43784f", linestyle="--", linewidth=1.2)
    ax.text(ctrl*.95, .43, f"GW250114  {ctrl:.0f}\nHeld-out control", ha="right",
            color="#35623f", fontsize=10)
    fig.savefig(args.results / "background.png", dpi=170, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()

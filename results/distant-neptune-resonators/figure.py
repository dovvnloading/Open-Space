"""Render the two resonant trajectories and their finite-horizon orbit support."""

import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from current import TARGETS, setup
from dynamics import WORK, integrate, angle

if __name__ == "__main__":
    fig, ax = plt.subplots(
        2, 2, figsize=(9.4, 5.7), gridspec_kw={"width_ratios": [1.7, 1]}
    )
    colors = ["#15596a", "#9b5025"]
    survival = list(csv.DictReader((WORK / "duration.csv").open()))
    for i, (target, color) in enumerate(zip(TARGETS, colors)):
        name, _, p, q = target
        for sign in (-1, 1):
            sim = setup(target, 0, sign * 0.5)
            t, track = integrate(sim, sign * 400000, 200)
            a = angle(track[:, 0], p, q)
            ax[i, 0].plot(t / 1e6, np.rad2deg(a), color=color, linewidth=0.8)
        ax[i, 0].set(
            xlim=(-0.4, 0.4),
            ylim=(0, 360),
            yticks=[0, 90, 180, 270, 360],
            ylabel="Critical angle (deg)",
            title=f"{name.replace('-', ' ')}  ·  {p}:{q} with Neptune",
        )
        records = [r for r in survival if r["object"] == name]
        x = [float(r["half_span_yr"]) / 1e6 for r in records]
        y = [100 * int(r["both_directions_bounded"]) / int(r["draws"]) for r in records]
        ax[i, 1].plot(x, y, "o-", ms=4, color=color, linewidth=1.6)
        ax[i, 1].set(
            xscale="log",
            xlim=(0.04, 12),
            ylim=(-4, 104),
            xticks=[0.1, 1, 10],
            xticklabels=["0.1", "1", "10"],
            ylabel="Bounded samples (%)",
        )
        for a in ax[i]:
            a.spines[["top", "right"]].set_visible(False)
            a.grid(axis="y", alpha=0.22)
    ax[0, 1].set_title("64 current-orbit draws per object")
    ax[1, 0].set_xlabel("Time from 2026 orbit epoch (Myr)")
    ax[1, 1].set_xlabel("Integration each side of epoch (Myr)")
    fig.tight_layout(pad=1.4)
    fig.savefig(Path(__file__).resolve().parent / "resonances.png", dpi=150)

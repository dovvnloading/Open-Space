"""Render amplitude measurements and the reconstruction trained before 2010."""

import csv
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from modulation import ROOT, design


def main():
    rows = [
        r
        for r in csv.DictReader((ROOT / "segments.csv").open())
        if r["star"].endswith("2332")
    ]
    report = json.loads((ROOT / "validation_2332.json").read_text())
    split = report["split"]
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": "#223340",
            "xtick.color": "#344752",
            "ytick.color": "#344752",
        }
    )
    fig, (ax, bx) = plt.subplots(
        2,
        1,
        figsize=(9.5, 6.5),
        sharex=True,
        layout="constrained",
        height_ratios=[1.2, 1],
    )
    dates = lambda t: [datetime(1995, 10, 9, 12) + timedelta(days=float(d)) for d in t]
    t = np.arange(280, 7600, 4)
    ax.plot(
        dates(t),
        1000 * (design(t, split["period"]) @ split["model_beta"]),
        color="#294858",
        lw=1.4,
        label="Cycle fitted to pre-2010 OGLE",
    )
    ax.axvspan(dates([4955])[0], dates([7594])[0], alpha=0.10, color="#5684bf")

    def series(label, idx):
        group = [r for r in rows if r["survey_passband"] == label]
        return dates([r["mean_HJD_minus_2450000"] for r in group]), [
            float(r[idx]) * 1000 for r in group
        ]

    for lab, c, marker, title in [
        ("OGLE II/III I", "#14686b", "o", "OGLE before 2010"),
        ("OGLE IV I", "#ad4827", "s", "OGLE from 2010"),
    ]:
        ax.scatter(
            *series(lab, "semiamp_2O_mag"),
            c=c,
            marker=marker,
            s=32,
            label=title,
            zorder=3,
        )
    for lab, c, m in [("EROS B", "#405ea6", "D"), ("EROS R", "#a65733", "o")]:
        x, y = series(lab, "semiamp_2O_mag")
        bx.plot(x, y, "-", color=c, lw=0.8, alpha=0.55)
        bx.scatter(x, y, color=c, marker=m, s=33, label=lab)
    bx.plot(
        dates(t[t < 2700]),
        1000 * (design(t[t < 2700], split["period"]) @ split["model_beta"]),
        color="#294858",
        lw=1.2,
        ls="--",
        label="OGLE model in I band",
    )
    for panel in [ax, bx]:
        panel.grid(axis="y", alpha=0.22)
        panel.set_ylim(-3, 100)
        panel.set_ylabel("Second-overtone\nsemi-amplitude (mmag)")
        panel.legend(loc="upper right", ncol=2, frameon=False, fontsize=9)
    ax.set_title("OGLE-SMC-CEP-2332", loc="left", weight="bold", fontsize=16)
    bx.set_xlabel("Observing year")
    fig.savefig(ROOT / "amplitude-history.png", dpi=145)


if __name__ == "__main__":
    main()

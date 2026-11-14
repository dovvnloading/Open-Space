"""Plot measured event depths and differences using the saved CSVs."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def rows(name):
    with (ROOT / name).open() as table:
        return list(csv.DictReader(table))


def main():
    plt.rcParams.update(
        {"font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "tess-sector108"}
    )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.1, 4.1),
        gridspec_kw={"width_ratios": [1.35, 1]},
        layout="constrained",
    )
    ax, control = axes
    observed = rows("followup_events.csv")
    styles = [
        ("SAP", "#216d8b", "o", "SAP"),
        ("core_3x3", "#2e8860", "s", "3 × 3 pixels"),
        ("PDCSAP", "#b64c36", "D", "PDCSAP"),
    ]
    for idx, (kind, color, symbol, label) in enumerate(styles):
        values = [
            r
            for r in observed
            if int(r["sector"]) == 108 and r["extraction"] == kind and r["depth_ppm"]
        ]
        ax.scatter(
            [int(r["cycle"]) + (idx - 1) * 0.13 for r in values],
            [float(r["depth_ppm"]) for r in values],
            color=color,
            marker=symbol,
            s=44,
            label=label,
            zorder=3,
        )
    ax.axhline(1779, color="#6d7882", lw=1.1, ls="--", label="Published depth")
    ax.axhline(0, color="#555", lw=0.8)
    ax.set(
        xticks=range(191, 197),
        xlabel="Predicted eclipse cycle",
        ylabel="Measured depth (ppm)",
        title="Four covered predicted events",
    )
    ax.legend(fontsize=8, frameon=False, ncol=2, loc="lower left")
    ax.set_ylim(-4200, 3900)
    star = rows("detector_events.csv")
    for n in [195, 196]:
        group = {}
        for r in star:
            if int(r["cycle"]) == n and r["depth_ppm"]:
                group.setdefault(r["tic"], {})[r["extraction"]] = float(r["depth_ppm"])
        controls = [v["PDCSAP"] - v["SAP"] for v in group.values() if len(v) == 2]
        control.scatter(
            n + np.linspace(-0.11, 0.11, len(controls)),
            sorted(controls),
            s=18,
            color="#87959f",
            label="Comparison stars" if n == 195 else None,
        )
        target = {
            r["extraction"]: float(r["depth_ppm"])
            for r in observed
            if int(r["sector"]) == 108 and int(r["cycle"]) == n and r["depth_ppm"]
        }
        control.scatter(
            n,
            target["PDCSAP_common"] - target["SAP_common"],
            s=56,
            color="#b64c36",
            marker="D",
            label="TIC 257433294" if n == 195 else None,
            zorder=3,
        )
    control.axhline(0, color="#555", lw=0.8)
    control.set(
        xticks=[195, 196],
        xlim=(194.6, 196.4),
        ylim=(-6300, 1600),
        xlabel="Predicted eclipse cycle",
        ylabel="PDCSAP minus SAP depth (ppm)",
        title="Identical-cadence comparison",
    )
    control.legend(frameon=False, fontsize=8, loc="lower left")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.16)
    fig.suptitle("TIC 257433294 · TESS Sector 108", fontsize=14, fontweight="bold")
    fig.savefig(ROOT / "sector108.svg", metadata={"Date": None})


if __name__ == "__main__":
    main()

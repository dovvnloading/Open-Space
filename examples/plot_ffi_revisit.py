"""Plot measured pixel amplitudes and independent period searches from the TESS revisit."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


def plot(packet, output):
    sectors = [s["product"]["sector"] for s in packet["sectors"]]
    rows = [
        next(c for c in s["cases"] if c["role"] == "unresolved_source")
        for s in packet["sectors"]
    ]
    maps = [np.array(c["fixed_frequency_signal_map_e_per_s"]) for c in rows]
    limit = max(np.abs(im).max() for im in maps)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(11.8, 6.4), facecolor="white")
    fig.text(
        0.055,
        0.951,
        "A repeating signal in the Large Magellanic Cloud",
        fontsize=22,
        weight="bold",
        color="#173547",
    )
    fig.text(
        0.055,
        0.907,
        "TESS full-frame images  •  ICRS 77.05502°, −69.23710°  •  individual star unresolved",
        fontsize=12,
        color="#375963",
    )
    fig.text(
        0.055,
        0.837,
        "Signed flux semiamplitude at the fixed 1.11221-day period",
        fontsize=13,
        weight="bold",
        color="#173547",
    )
    gs = fig.add_gridspec(
        1,
        len(sectors) + 1,
        left=0.055,
        right=0.92,
        bottom=0.545,
        top=0.785,
        wspace=0.26,
        width_ratios=[1] * len(sectors) + [0.09],
    )
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    for col, (sector, im) in enumerate(zip(sectors, maps)):
        ax = fig.add_subplot(gs[0, col])
        img = ax.imshow(
            im, origin="lower", interpolation="nearest", norm=norm, cmap="RdBu_r"
        )
        ax.set_title(f"Sector {sector}", fontsize=10)
        ax.plot([2], [2], marker="+", color="#111111", ms=11, mew=1)
        ax.set_xticks([])
        ax.set_yticks([])
    bar = fig.colorbar(img, cax=fig.add_subplot(gs[0, -1]))
    bar.set_label("electrons / s", fontsize=9)
    fig.text(
        0.055,
        0.504,
        "Each stamp is 5 × 5 detector pixels. The cross marks the pixel nearest the original sky position. Detector orientation varies by sector.",
        fontsize=9,
        color="#415663",
    )
    gs2 = fig.add_gridspec(
        1, 2, left=0.083, right=0.93, bottom=0.135, top=0.388, wspace=0.32
    )
    ax = fig.add_subplot(gs2[0, 0])
    pos = np.arange(len(sectors))
    periods = [c["nearest_blind_peak"]["period_days"] for c in rows]
    offsets = [c["nearest_blind_offset_arcsec"] for c in rows]
    colors = ["#bf5c32" if s in (1, 2) else "#077b78" for s in sectors]
    ax.axhline(rows[0]["fixed_period_days"], color="#a9b8ba", lw=1, ls="--", zorder=0)
    ax.scatter(pos, periods, s=52, color=colors, zorder=3)
    ax.set(ylabel="Period (days)", ylim=(1.086, 1.129), yticks=[1.09, 1.10, 1.11, 1.12])
    ax.set_title(
        "Period found independently in each sector",
        loc="left",
        fontsize=11,
        weight="bold",
    )
    bx = fig.add_subplot(gs2[0, 1])
    bx.bar(pos, offsets, color=colors, width=0.43)
    bx.set(ylabel="Offset (arcseconds)", ylim=(0, 16), yticks=[0, 4, 8, 12, 16])
    bx.set_title(
        "Position relative to the sector 1 detection",
        loc="left",
        fontsize=11,
        weight="bold",
    )
    for axis in (ax, bx):
        axis.set_xticks(pos, labels=[str(v) for v in sectors])
        axis.set_xlabel("TESS sector")
        axis.grid(axis="y", alpha=0.14)
        axis.spines["left"].set_color("#9cabab")
        axis.spines["bottom"].set_color("#9cabab")
    fig.text(
        0.055,
        0.033,
        "Orange marks the original two sectors. Green marks the three subsequent checks. Peak positions and periods come from blind spatial searches.",
        fontsize=9,
        color="#415663",
    )
    fig.savefig(output, dpi=160, metadata={"Date": None})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(json.loads(args.input.read_bytes()), args.output)


if __name__ == "__main__":
    main()

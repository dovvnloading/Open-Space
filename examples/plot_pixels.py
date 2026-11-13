"""Render measured TESS postage stamps and their dimming-source localization."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurements", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--tic", type=int, help="Plot one target in a multi-target report")
    args = parser.parse_args()
    loaded = json.loads(args.measurements.read_text())
    records = loaded if isinstance(loaded, list) else [loaded]
    records = [r for r in records if r.get("images") and (args.tic is None or r["tic_id"] == args.tic)]
    if not records:
        parser.error("No usable difference images")
    plt.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "tess-difference-v1",
                         "font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(len(records), 3, figsize=(10.4, 2.9*len(records)+1.65),
                             squeeze=False, layout="constrained")
    title = ("TIC " + str(records[0]["tic_id"])) if len({r["tic_id"] for r in records}) == 1 else "TESS targets"
    fig.suptitle(f"{title} · Locating the dimming source", fontsize=18, fontweight="bold", x=0.055, ha="left")
    for i, (row, record) in enumerate(zip(axes, records)):
        loc = record["localization"]
        for j, (ax, key) in enumerate(zip(row, ("scene", "difference", "score"))):
            data = np.asarray(record["images"][key], dtype=float)
            if key == "scene":
                finite = data[np.isfinite(data)]
                vmin, vmax = np.percentile(finite, [4, 100])
                norm = colors.AsinhNorm(linear_width=max(1, (vmax-vmin)/70), vmin=vmin, vmax=vmax)
                cmap = "magma"
            else:
                vmax = max(abs(np.nanmin(data)), abs(np.nanmax(data)), 1.)
                norm = colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
                cmap = "RdBu_r"
            img = ax.imshow(data, origin="lower", interpolation="nearest", cmap=cmap, norm=norm)
            if i == 0:
                ax.set_title(("Out of event", "Out minus in", "Difference / error scale")[j], fontsize=12)
            mask = np.asarray(record["science_aperture"])
            ax.contour(mask, levels=[0.5], colors=["#a3a4a6"], linewidths=.8, linestyles="--")
            ax.scatter(*record["target_xy"], marker="+", s=90, c="#60e5d8", linewidths=1.9)
            if loc["centroid_xy"] is not None:
                ax.scatter(*loc["centroid_xy"], marker="o", s=95, facecolors="none",
                           edgecolors="#282027" if j else "white", linewidths=1.5)
            ax.tick_params(length=2, labelsize=8)
            ax.set_xlabel("Stamp column (pixel)")
            if j == 0:
                ax.set_ylabel(f"Sector {record['sector']} · CCD {record['ccd']}\nStamp row (pixel)")
            cb = fig.colorbar(img, ax=ax, orientation="horizontal", pad=.025, shrink=.78, aspect=26)
            cb.set_label("electrons / s" if j < 2 else "pointwise score", fontsize=8)
            cb.ax.tick_params(length=2, labelsize=7)
    legend = [Line2D([], [], marker="+", color="#008f82", lw=0, markersize=9, label="TIC sky position"),
              Line2D([], [], marker="o", markerfacecolor="none", color="#282027", lw=0, markersize=8,
                     label="Difference-image centroid"),
              Line2D([], [], ls="--", color="#777778", label="SPOC aperture")]
    fig.legend(handles=legend, loc="outside lower center", ncol=3, frameon=False)
    fig.savefig(args.output, dpi=145, metadata={"Date": None} if args.output.suffix.lower()==".svg" else None)
    plt.close(fig)


if __name__ == "__main__":
    main()

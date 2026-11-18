"""Render the optical shoreline checks and the verified calibration difference."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from screen import BOUNDS, VENT, run


def render(data, output):
    findings, plots = run(data, "zero")
    scenes = (
        ("l9-20261111", "11 Nov · Landsat 9"),
        ("s2c-20261111", "11 Nov · Sentinel-2C"),
        ("s2b-20261116", "16 Nov · Sentinel-2B"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(10.2, 4), layout="constrained")
    figure.suptitle("Okinahama after correcting the Sentinel-2 reflectance offset", fontsize=14, x=0.02, ha="left")
    extent = (BOUNDS[0], BOUNDS[2], BOUNDS[1], BOUNDS[3])
    for ax, (name, title) in zip(axes, scenes):
        frame, grid, masks, _ = plots[name]
        values = np.stack([frame.values[key] for key in ("red", "green", "blue")], axis=2)
        color = np.clip(np.nan_to_num(values) * 3.2, 0, 1) ** 0.65
        ax.imshow(color, extent=extent)
        if masks["new"].any():
            ax.contour(grid.x, grid.y, masks["new"], levels=[0.5], colors=["#38ffd2"], linewidths=1.15)
        ax.add_patch(Circle(VENT, 450, fill=False, color="#fbed61", linestyle="--", linewidth=1.05))
        ax.set_xlim(VENT[0] - 1100, VENT[0] + 1100)
        ax.set_ylim(VENT[1] - 1100, VENT[1] + 1100)
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[:].set_visible(False)
    figure.supxlabel("Dashed circle · 450 m around approximate vent. Cyan · nominal new-ground proxy.", fontsize=9)
    figure.savefig(output, dpi=155)
    plt.close(figure)
    return findings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == "__main__":
    main()

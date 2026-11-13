"""Plot per-window coverage, sign-matched detections and model-subtracted outcomes."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("holdout", type=Path)
    p.add_argument("template", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    data = json.loads(args.holdout.read_bytes())["summary"]
    templates = {r["id"]: r for r in json.loads(args.template.read_bytes())["summary"]}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "open-space-holdout-v1"})
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    fig.subplots_adjust(left=0.275, right=0.928, top=0.80, bottom=0.16)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#d9e2e7", linewidth=0.8)
    labels = []
    for i, row in enumerate(data):
        num = len(data) - i
        hit, coverage = row["primary"]["hits"], row["primary"]["covered"]
        ax.barh(num, coverage, height=0.62, color="#d6e3e9")
        ax.barh(num, hit, height=0.62, color="#205b78")
        remaining = templates[row["id"]]["residual_detections"]
        ax.scatter(remaining, num, facecolor="#e17c41", edgecolor="#934015", s=31, zorder=5)
        ax.text(coverage + 0.45, num, f"{hit}/{coverage}", va="center", fontsize=9, color="#1d343f")
        group = row["id"] if row["id"] != "CaK" else "Ca K"
        direction = "+" if row["direction"] == "emission" else "−"
        labels.append(f"{group}    {row['vacuum_rest_angstrom']:,.2f} Å  ({direction})")
    ax.set_yticks(list(range(len(data), 0, -1)), labels)
    ax.set_xlim(-0.65, 27.7)
    ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
    ax.tick_params(axis="both", length=0, pad=7, labelcolor="#203541")
    ax.set_xlabel("Number of held-out galaxies", labelpad=9)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(0.054, 0.95, "Frozen predictions in 24 new galaxy spectra", color="#183440",
             fontsize=19, fontweight="bold", va="top")
    fig.text(0.055, 0.876, "Three SDSS plate visits  •  Rest-frame vacuum wavelengths  •  (+) emission, (−) absorption",
             color="#435762", fontsize=10, va="top")
    handles = [Line2D([], [], lw=7, color="#d6e3e9", label="Searchable"),
               Line2D([], [], lw=7, color="#205b78", label="Detected"),
               Line2D([], [], linestyle="", marker="o", color="#e17c41", label="Survived model subtraction")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.05, 0.013),
               ncol=3, frameon=False, fontsize=9.5)
    fig.savefig(args.output, metadata={"Date": None} if args.output.suffix == ".svg" else None,
                dpi=150, facecolor="white")


if __name__ == "__main__":
    main()

"""Plot JEDI witness counts around the frozen encounter times."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def plot(source, output):
    passes = json.loads(source.read_text(encoding="utf-8"))["passes"]
    titles = [("83", "May 1 · Closest approach to Thebe"),
              ("82", "March 29 · Closest approach to Jupiter")]
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 5.4), sharex=True, layout="constrained")
    for axis, (orbit, title) in zip(axes, titles):
        data = passes[orbit]["jedi_primary"]["bins"]
        time = np.array([row["offset_s"] / 60 for row in data])
        rate = np.array([row["cps"][:2] if row["cps"] else [np.nan, np.nan] for row in data])
        rate[rate <= 0] = np.nan
        axis.axvspan(-1, 1, color="#D8E2ED", label="Test window")
        axis.axvspan(-20, -8, color="#EAECEF", label="Baseline")
        axis.axvspan(8, 20, color="#EAECEF")
        axis.axvline(0, color="#526376", linewidth=0.8, linestyle="--")
        axis.plot(time, rate[:, 0], color="#365DA7", linewidth=1.35, label="JEDI T1")
        axis.plot(time, rate[:, 1], color="#B64727", linewidth=1.35, label="JEDI T3")
        axis.set(title=title, ylabel="Counts per second", xlim=(-45, 45), yscale="log")
        axis.grid(axis="y", linewidth=0.5, alpha=0.25)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=4, fontsize=8, loc="lower left")
    axes[1].set_xlabel("Minutes from closest approach")
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=ROOT / "summary.json")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    plot(arguments.summary, arguments.output)

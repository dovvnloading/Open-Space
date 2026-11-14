"""Plot the independent 2024 test from the released event measurements."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    folder = Path(__file__).resolve().parent
    events = np.genfromtxt(folder / "events.csv", names=True, delimiter=",")
    later = events[events["sector"] >= 80]
    epoch = later["orbit"]
    design = np.column_stack((np.ones(len(epoch)), epoch - epoch.mean()))

    def residual(column):
        y = later[column]
        return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]

    predicted = residual("predicted_bias_min")
    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8.8, 4.4), layout="constrained")
    colors = ["#186a73", "#b24a34"]
    for ax, label, col in zip(axes, ("Original flux", "After binary subtraction"),
                              ("raw_min", "corrected_min")):
        measured = residual(col)
        ax.axhline(0, lw=0.9, color="#a5a5a5")
        ax.plot([-19, 19], [-19, 19], ls="--", color="#9b9b9b", lw=1)
        for sec, color in zip((80, 81), colors):
            indices = later["sector"] == sec
            ax.scatter(predicted[indices], measured[indices], color=color, label=f"Sector {sec}", s=40)
        ax.set_title(label, loc="left", fontsize=12)
        ax.set(xlim=(-19, 19), ylim=(-19, 19), aspect="equal", xticks=[-15, -10, -5, 0, 5, 10, 15])
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, loc="upper left", fontsize=9)
    axes[0].set_ylabel("Measured timing residual (minutes)")
    fig.supxlabel("Timing bias predicted from the Gaia binary waveform (minutes)", fontsize=11)
    fig.suptitle("KELT-8 b  |  Twelve TESS transits in 2024", fontsize=14)
    fig.savefig(folder / "later-transits.png", dpi=150)


if __name__ == "__main__":
    main()

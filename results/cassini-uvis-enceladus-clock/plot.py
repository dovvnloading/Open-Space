"""Plot the difference caused by mixing the Cassini UVIS archive clocks."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze import HERE, BASELINE, window


def run(output):
    time, fuv, legacy, mixed, sclk = np.loadtxt(HERE / "curve.csv", delimiter=",", skiprows=1).T
    selected = window(time, (535, 645))
    ref = window(time, BASELINE)
    seconds = time[selected] - 540
    reference = np.mean(sclk[ref])
    fuv_norm = fuv / np.nanmean(fuv[ref])
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8.6, 5.6), layout="constrained",
                             gridspec_kw={"height_ratios": (1.65, 1)})
    top, bottom = axes
    top.set_title("Mixing archive clocks changes which UVIS counts share an exposure", loc="left", fontsize=14)
    top.plot(seconds, 100 * fuv_norm[selected], ".-", color="#295585", alpha=.7, label="FUV spectrum sum")
    top.plot(seconds, 100 * sclk[selected] / reference, color="#08746b", label="HSP with a common clock")
    top.plot(seconds, 100 * mixed[selected] / reference, color="#af6d28", linestyle=":", label="HSP with mixed archive clocks")
    top.set_ylabel("Counts (% of reference)")
    top.legend(loc="lower right", facecolor="white", framealpha=0.97, fontsize=8)
    diff = sclk[selected] - mixed[selected]
    bottom.plot(seconds, diff, ".-", color="#96562e", linewidth=1.1)
    bottom.axhline(0, color="#c5cdd3")
    peak = np.argmax(np.abs(diff) * (seconds >= 10))
    bottom.scatter([seconds[peak]], [diff[peak]], color="#96562e", s=27)
    bottom.annotate(f"{diff[peak]:.0f} counts", (seconds[peak], diff[peak]),
                    xytext=(13, 0), textcoords="offset points", fontsize=9)
    bottom.set_ylabel("Common minus mixed\nHSP counts")
    bottom.set_xlabel("Seconds since 14:39:00 UTC, FUV exposure midpoints")
    for ax in axes:
        ax.grid(axis="y", color="#e4e9ec")
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output, dpi=145)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    run(parser.parse_args().output)

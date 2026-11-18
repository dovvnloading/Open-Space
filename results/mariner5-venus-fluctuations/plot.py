"""Plot the known Venus crossings and the matched cruise-window comparison."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from analyze import ROOT, encounter


BLUE = "#175d84"
TEAL = "#10775f"
INK = "#172f3d"
GRAY = "#657b86"


def plot(output):
    frames, _ = encounter()
    result = json.loads((ROOT / "summary.json").read_bytes())
    fig, (inbound, outbound, controls) = plt.subplots(1, 3, figsize=(12.9, 4.3), layout="constrained")
    for ax, event, metric, index, color, title in (
        (inbound, "inbound", "magnitude", 0, BLUE, "Inbound field magnitude"),
        (outbound, "outbound", "vector_variability", 1, TEAL, "Outbound vector variability"),
    ):
        center = -156 if event == "inbound" else 18
        peak = result["events"][event]["primary"]["maxima"][metric]
        x_peak = peak["split_s"] / 60
        x = [frame.second / 60 for frame in frames if abs(frame.second / 60 - center) <= 19]
        y = [frame.values[index] for frame in frames if abs(frame.second / 60 - center) <= 19]
        ax.scatter(x, y, color=color, s=18, alpha=.72)
        ax.axvline(center, color=GRAY, linewidth=1.15, linestyle=":", label="Historical estimate")
        ax.axvline(x_peak, color=INK, linewidth=1.2, linestyle="--", label="Peak split")
        for key, start, stop in (("median_before_nt", -8, -2), ("median_after_nt", 2, 8)):
            ax.axvspan(x_peak + start, x_peak + stop, color=color, alpha=.1, linewidth=0)
            ax.hlines(peak[key], x_peak + start, x_peak + stop, color=color, linewidth=3)
        ylabel = "Field magnitude (nT)" if index == 0 else "Within-frame deviation (nT)"
        ax.set(title=title, xlabel="Minutes from closest approach", ylabel=ylabel)
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        ax.legend(loc="upper left" if index == 0 else "upper right", frameon=False, fontsize=8)

    metric = "vector_variability"
    peaks = sorted(item["scan"]["maxima"][metric]["abs_log_ratio"]
                   for item in result["cruise_anchors"] if item["scan"]["maxima"][metric])
    peak = result["events"]["outbound"]["primary"]["maxima"][metric]["abs_log_ratio"]
    p95 = result["cruise_score_context"][metric]["p95_nearest_rank"]
    controls.step([0, *peaks, 1.9], [0, *range(1, len(peaks) + 1), len(peaks)],
                  where="post", color=GRAY, linewidth=1.7)
    controls.axvline(p95, color=GRAY, linestyle=":", linewidth=1.1)
    controls.axvline(peak, color=TEAL, linewidth=2)
    controls.text(p95 + .035, 73, "Cruise 95th\npercentile", va="top", ha="left", color=GRAY, fontsize=8)
    controls.text(peak - .035, 49, "Venus\noutbound", va="top", ha="right", color=TEAL, fontsize=8)
    controls.set(title="Cruise-window comparison", xlabel="Peak absolute log ratio", ylabel="Cumulative cruise windows")
    for ax in (inbound, outbound, controls):
        ax.title.set_color(INK)
        ax.title.set_fontweight("bold")
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e5ebee", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#a5b3b9")
        ax.tick_params(labelsize=9, colors=INK)
    fig.savefig(output, dpi=170, metadata={"Software": "Open-Space Mariner 5 analysis"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    plot(parser.parse_args().output)

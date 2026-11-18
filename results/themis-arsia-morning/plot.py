"""Plot the unshifted profiles and the registration stress test."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def render(source, output):
    pair = json.loads(source.read_text(encoding="utf-8"))["pairs"]["primary"]
    main = pair["maps"]["main_720m_linear"]
    blocks = main["blocks"]
    latitudes = [row["latitude_deg"] for row in blocks]
    contrast = np.array([row["median_change"] or [np.nan, np.nan] for row in blocks])
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.2), layout="constrained",
                             gridspec_kw={"height_ratios": [2, 1]})
    profile, shifts = axes
    profile.axvspan(-10, -8.5, color="#E6EDF5", label="Search interval")
    profile.axvspan(-10.55, -10.05, color="#EEEEEE", label="Background strips")
    profile.axvspan(-8.45, -7.95, color="#EEEEEE")
    profile.axhline(10, color="#526376", linestyle="--", linewidth=1, label="Threshold")
    profile.axhline(0, color="#A3A8A9", linewidth=.6)
    for band, color, values in zip(("0.42 µm", "0.75 µm"), ("#365DA7", "#B64727"), contrast.T):
        profile.plot(latitudes, values * 100, color=color, linewidth=1.6, label=band)
    profile.set(xlabel="Planetocentric latitude (°)", ylabel="Relative change (%)", xlim=(-10.6, -7.9),
                title="January brightness after estimating the broad change")
    profile.legend(ncol=3, fontsize=8, loc="upper left", frameon=False)
    labels = ["Unshifted", "East", "West", "North", "South"]
    runs = [main["longest_coincident_run_blocks"]]
    for label in labels[1:]:
        runs.append(pair["registration_sensitivity"][label.lower()]["longest_coincident_run_blocks"])
    shifts.bar(labels, runs, color=["#B64727" if count >= 6 else "#527399" for count in runs])
    shifts.axhline(6, linestyle="--", color="#526376", linewidth=1)
    shifts.set(ylim=(0, 8), ylabel="Consecutive blocks", title="Shifting January by one 720 m cell")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=.18)
    fig.savefig(output, dpi=165)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=ROOT / "summary.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.summary, args.output)

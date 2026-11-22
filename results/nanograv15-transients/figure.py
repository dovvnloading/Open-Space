"""Render the survey-wide calibration figure from deposited workflow evidence."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter


TEAL = "#176b77"
ORANGE = "#b54b12"
INK = "#182838"
MUTED = "#556676"


def records(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def panel(ax, heading, detail):
    ax.set_title(heading, loc="left", color=INK, fontsize=13, fontweight="bold", pad=32)
    ax.text(0, 1.025, detail, transform=ax.transAxes, ha="left", va="bottom",
            color=MUTED, fontsize=9)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#86959e")
    ax.grid(axis="y", alpha=.15)
    ax.set_axisbelow(True)


def sign_panel(ax, evidence):
    path = evidence / "coherence"
    summary = json.loads((path / "summary.json").read_text())
    maxima = [float(row["full_max_log10_score"]) for row in records(path / "null-maxima.csv")]
    draws = summary["draws"]
    if len(maxima) != draws:
        raise ValueError("The sign trial list differs from the recorded draw count")
    tail = summary["observed_full_max_tail"]
    observed = summary["observed_full_max_log10_score"]
    panel(ax, "A   Whole-pulsar sign reversals",
          f"{draws:,} draws   ·   {summary['n_bank']:,} templates per draw")
    ax.hist(maxima, bins=36, color=TEAL, alpha=.83, edgecolor="white", linewidth=.4)
    ax.axvline(observed, color=ORANGE, linewidth=2, zorder=4,
               label=f"Observed maximum  {observed:.2f}")
    ax.legend(loc="upper right", framealpha=.94, edgecolor="none", fontsize=9)
    ax.set_xlabel(r"Full-bank maximum of $-\log_{10}(p_{\rm local})$")
    ax.set_ylabel("Sign draws per bin")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    return (f"Survey-wide sign-test p = {tail['p_sign_max_bank']:.4g}. "
            f"{tail['exceedances']:,} draws met or exceeded the observed maximum.")


def memory_panel(ax, evidence):
    path = evidence / "array"
    rows = [r for r in records(path / "array_limits.csv")
            if r["family"] == "memory" and r["shape"] == "ramp" and r["units"] == "strain"]
    rows.sort(key=lambda r: float(r["amplitude"]))
    if not rows or float(rows[0]["amplitude"]) != 0 or len(rows) < 2:
        raise ValueError("The memory injection ladder or its zero-amplitude control is missing")
    counts = {int(r["n_trials"]) for r in rows}
    if len(counts) != 1:
        raise ValueError("Expected equal injection counts at each amplitude")
    n = counts.pop()
    control, positive = rows[0], rows[1:]
    x = [float(row["amplitude"]) for row in positive]
    y = [float(row["recovery_fraction"]) for row in positive]
    low = [float(row["lower_95"]) for row in positive]
    sensitivity = json.loads((path / "array_summary.json").read_text())["sensitivity"]
    entry = next(r for r in sensitivity if r["family"] == "memory" and r["shape"] == "ramp")
    threshold = entry["first_tested_amplitude_90pct"]
    panel(ax, "B   Memory-strain injection recovery",
          f"{n} trials per amplitude   ·   eligible cells, on- and off-grid epochs")
    ax.set_xscale("log")
    ax.axhline(.9, color=MUTED, linestyle=(0, (3, 3)), linewidth=1)
    ax.text(.02, .9, "90%", transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            color=MUTED, fontsize=8)
    ax.vlines(x, low, y, color=TEAL, linewidth=1, alpha=.65)
    ax.scatter(x, low, marker="_", color=TEAL, s=26, zorder=3, label="Simultaneous 95% lower bound")
    ax.scatter(x, y, color=TEAL, s=17, zorder=4, label="Recovered fraction")
    if threshold is not None:
        i = x.index(float(threshold))
        ax.scatter(x[i], low[i], facecolor="white", edgecolor=ORANGE, s=46,
                   linewidth=1.5, zorder=5, label="First bound reaching 90%")
    ax.set_ylim(-.035, 1.04)
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.yaxis.set_major_formatter(PercentFormatter(1.))
    ax.set_xlabel("Injected Earth-term memory strain (dimensionless)")
    ax.set_ylabel("Recovered injections")
    ax.legend(loc="center left", framealpha=.95, edgecolor="none", fontsize=8)
    first = (f"First strain reaching the 90% recovery criterion was {threshold:.3g}."
             if threshold is not None else "No tested amplitude met the criterion.")
    return f"{first} Zero-strain control yielded {control['n_recovered']}/{n} recovered."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    args = parser.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelpad": 8})
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.8, 5.1))
    fig.subplots_adjust(left=.079, right=.984, bottom=.28, top=.84, wspace=.32)
    sign_note = sign_panel(left, args.evidence)
    memory_note = memory_panel(right, args.evidence)
    fig.text(.079, .104, sign_note, color=INK, fontsize=9)
    fig.text(.079, .061, memory_note, color=INK, fontsize=9)
    fig.savefig(args.evidence / "search-calibration.png", dpi=170, facecolor="white",
                metadata={"Title": "NANOGrav transient search calibration"})
    plt.close(fig)


if __name__ == "__main__":
    main()

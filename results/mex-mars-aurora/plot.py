"""Plot the three near-arrival limb profiles and native calibration coverage."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze import ROOT, SKY, TARGET, groups, inside, load, model


def plot(output):
    pairs = load()
    fig, (limb, audit) = plt.subplots(1, 2, figsize=(10.2, 4.25), layout="constrained",
                                   gridspec_kw={"width_ratios": [1.65, 1]})
    colors = ("#235b83", "#687985", "#a34e31")
    descriptions = (
        ("4432-tcl2", "Mar 17 · 04:29 · New Norcia"),
        ("4433-dsn", "Mar 17 · 11:14 · Madrid"),
        ("4434-tcl2", "Mar 18 · 01:27 · New Norcia"),
    )
    for (key, title), color in zip(descriptions, colors):
        ingress = next(block for block in groups(pairs[key])
                       if block[0]["direction"] == "ingress" and len(inside(block, SKY)) >= 15)
        origin, offset, slope = model(inside(ingress, SKY), "d_hz", "linear")
        plot_rows = inside(ingress, (110, 680))
        x = [row["altitude"] for row in plot_rows]
        y = [row["d_hz"] - offset - slope * (row["second"] - origin) for row in plot_rows]
        limb.plot(x, y, color=color, linewidth=1.0, alpha=.85, label=title)
    limb.axvspan(*TARGET, color="#e8cda8", alpha=.45, linewidth=0)
    limb.axvspan(*SKY, color="#b0c5da", alpha=.25, linewidth=0)
    limb.axhline(0, color="#637885", linewidth=.7, linestyle=":")
    limb.set(title="Differential Doppler during ingress", xlabel="Spherical tangent height (km)",
             ylabel="Background-adjusted S − (3/11)X (Hz)")
    limb.legend(loc="upper right", frameon=True, fontsize=8)
    limb.set_xlim(110, 680)

    data = json.loads((ROOT / "radio-reduction.json").read_bytes())["groups"]
    volumes = ["4433", "4436", "4443", "4445", "4451", "4452"]
    values = []
    for key in volumes:
        rows = data[f"{key}-dsn"]["s"]
        values.append(100 * (1 - rows["missing_reported_residual_rows"] / rows["native_rows"]))
    y = list(range(len(volumes)))
    audit.barh(y, [100] * len(y), color="#e7edef")
    audit.barh(y, values, color="#235b83")
    audit.set_yticks(y, labels=volumes)
    audit.invert_yaxis()
    audit.set(title="Restored DSN S-band tables", xlabel="Rows with a calibrated residual (%)", ylabel="PDS volume")
    for i, pct in enumerate(values):
        text = f"{pct:.2f}%" if 0 < pct < 100 else f"{pct:.0f}%"
        audit.text(max(pct - 2, 2), i, text, ha="right" if pct > 10 else "left",
                   va="center", color="white" if pct > 10 else "#172e3b", fontsize=8)
    for ax in (limb, audit):
        ax.set_axisbelow(True)
        ax.title.set_color("#172e3b")
        ax.title.set_fontweight("bold")
        ax.spines[["top", "right"]].set_visible(False)
    limb.grid(axis="y", color="#e7edef", linewidth=.8)
    audit.grid(axis="x", color="#e7edef", linewidth=.8)
    fig.savefig(output, dpi=180, metadata={"Software": "Open-Space MaRS analysis"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.output)

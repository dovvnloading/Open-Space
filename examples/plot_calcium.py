"""Draw the three measured SDSS galaxy coadds used in the first investigation."""

import argparse
import math
from pathlib import Path

from examples.calcium_doublet import K
from open_space.dataio import load_observation


def plot(cohort, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "svg.fonttype": "none",
        "svg.hashsalt": "open-space-calcium-k-v1",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": "#394b57",
        "ytick.color": "#394b57",
    })
    ink, flux_color, band = "#152e40", "#1d5964", "#d4e8e8"
    flag, reference = "#ad480d", "#805775"
    figure, axes = plt.subplots(3, 1, figsize=(8.4, 6.25), sharex=True, sharey=True)
    figure.subplots_adjust(left=0.115, right=0.965, top=0.77, bottom=0.11, hspace=0.17)
    figure.text(0.115, 0.944, "Absorption around Ca II K", fontsize=19, color=ink,
                weight="semibold")
    figure.text(0.115, 0.897, "Three SDSS DR17 galaxy coadds   ·   Plate 6497   ·   MJD 56329",
                fontsize=10, color="#405969")
    legend = [
        Line2D([], [], color=flux_color, lw=1.4, marker=".", markersize=4, label="Measured flux"),
        Patch(facecolor=band, label="±1σ flux error"),
        Line2D([], [], color=flag, marker="o", markerfacecolor="white", lw=0,
               markersize=4, label="OR flag set"),
        Line2D([], [], color=reference, lw=1.2, linestyle="--", label="Ca II K  3934.775 Å"),
    ]
    figure.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.108, 0.863),
                  frameon=False, ncol=4, borderaxespad=0, columnspacing=1.8,
                  handlelength=1.8, fontsize=8.4)

    for axis, fiber in zip(axes, (88, 96, 100)):
        manifest, rows = load_observation(cohort / f"sdss-6497-56329-{fiber:04d}.json")
        z = float(manifest["redshift"])
        waves, fluxes, errors, flagged = [], [], [], []
        for row in rows:
            wave = float(row["wavelength"]) / (1 + z)
            if not 3913 <= wave <= 3957:
                continue
            flux, weight = float(row["flux"]), float(row["ivar"])
            valid = int(row["mask"]) == 0 and math.isfinite(flux) and math.isfinite(weight) and weight > 0
            waves.append(wave)
            fluxes.append(flux if valid else math.nan)
            errors.append(1 / math.sqrt(weight) if valid else math.nan)
            if valid and int(row["or_mask"]) != 0:
                flagged.append((wave, flux))
        low = [f - e for f, e in zip(fluxes, errors)]
        high = [f + e for f, e in zip(fluxes, errors)]
        axis.fill_between(waves, low, high, color=band, lw=0, zorder=1)
        axis.axvline(K, color=reference, linestyle="--", linewidth=1.0, zorder=2)
        axis.plot(waves, fluxes, color=flux_color, linewidth=1.25, marker=".", markersize=3.2, zorder=3)
        if flagged:
            x, y = zip(*flagged)
            axis.scatter(x, y, edgecolors=flag, facecolors="white", s=13, linewidths=0.9, zorder=4)
        axis.text(0.015, 0.9, f"Fiber {fiber:04d}     z = {z:.6f}",
                  transform=axis.transAxes, fontsize=8.5, color=ink,
                  bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.1})
        axis.set_xlim(3914, 3956)
        axis.set_ylim(4, 25.5)
        axis.set_yticks([5, 10, 15, 20])
        axis.set_xticks([3915, 3920, 3925, 3930, 3935, 3940, 3945, 3950, 3955])
        axis.grid(axis="y", color="#e2e8e9", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.tick_params(length=3, width=0.7)
        for side in ("left", "bottom"):
            axis.spines[side].set_color("#97a9ac")
            axis.spines[side].set_linewidth(0.7)
    axes[-1].set_xlabel("Rest wavelength in vacuum (Å)", color=ink, labelpad=7)
    figure.text(0.035, 0.44, "Observed flux density  (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)",
                rotation=90, va="center", fontsize=10, color=ink)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"Creator": "Open-Space", "Title": "Absorption around Ca II K"}
    if output.suffix.lower() == ".svg":
        metadata["Date"] = None
    figure.savefig(output, dpi=170, metadata=metadata)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort", type=Path, help="directory of fetched investigation spectra")
    parser.add_argument("--output", type=Path, default=Path("docs/figures/calcium-k.svg"))
    args = parser.parse_args()
    plot(args.cohort, args.output)


if __name__ == "__main__":
    main()

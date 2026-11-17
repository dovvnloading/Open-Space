"""Plot empirical F444W limits from the five published curves."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt


fields = {
    "MCC 124": "mcc124",
    "G 99-49": "g99_49",
    "L 788-34": "l788_34",
    "Ross 619": "ross619",
    "V374 Peg": "v374peg",
}
directory = Path(__file__).resolve().parent
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.5), layout="constrained")
for label, stem in fields.items():
    with (directory / (stem + "_contrast.csv")).open() as stream:
        rows = list(csv.DictReader(stream))
    distance = [float(r["radius_arcsec"]) for r in rows]
    for ax, key in zip(axes, ["f95_uJy", "contrast95"]):
        ax.plot(
            distance,
            [float(r[key]) if r[key] else float("nan") for r in rows],
            marker="o",
            ms=4,
            label=label,
        )
for ax in axes:
    ax.set_yscale("log")
    ax.set_xlabel("Projected separation (arcsec)")
    ax.grid(alpha=0.2, which="both")
axes[0].set_ylabel("95% sampled recovery flux (µJy)")
axes[0].set_yticks([0.5, 1, 2, 4, 6], labels=["0.5", "1", "2", "4", "6"])
axes[1].set_yticks(
    [v * 1e-6 for v in (1.5, 2, 3, 5, 10, 20)],
    labels=["1.5", "2", "3", "5", "10", "20"],
)
axes[1].set_ylabel("95% sampled recovery contrast (×10⁻⁶)")
axes[0].legend(loc="upper right", frameon=False, fontsize=8, ncol=2)
fig.savefig(directory / "nircam_contrast.png", dpi=160)

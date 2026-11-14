"""Compare the fixed signed timing waveform with phase shifts on independent surveys."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def series(ax, x, y, color, label=None):
    x = (x + 0.5) % 1 - 0.5
    idx = np.argsort(x)
    ax.plot(x[idx], y[idx], color=color, lw=2, label=label)
    ax.scatter([0], y[np.argmin(abs(x))], color=color, s=35, zorder=5)


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharex=True)
    fig.subplots_adjust(bottom=0.23, top=0.79, left=0.085, right=0.972, wspace=0.32)
    fig.suptitle(
        "Do the independent light curves favor Gaia’s orbital phase?",
        fontsize=15,
        x=0.055,
        ha="left",
        y=0.97,
    )
    kelt = pd.read_csv(ROOT / "kelt/phase-controls.csv").query(
        "degree == 2 and reduction == 'tfa'"
    )
    result = json.loads((ROOT / "kelt/fits.json").read_text())["tfa"][1]
    series(
        axes[0],
        kelt["shift"].to_numpy(),
        result["zero_rms_ppt"] ** 2 - kelt.held_out_mse_ppt2.to_numpy(),
        "#166864",
    )
    ground = pd.read_csv(ROOT.parent / "ground-check/ground-phase-controls.csv")
    for survey, color in (("ASAS-3", "#4066a2"), ("ASAS-SN", "#ac5320")):
        d = ground[ground.survey == survey]
        series(
            axes[1],
            d.phase_shift.to_numpy(),
            d.frozen_improvement.to_numpy(),
            color,
            survey,
        )
    axes[1].legend(frameon=False, loc="lower right")
    for ax, title, ylabel in zip(
        axes,
        ("KELT-South · held-out nights", "ASAS · survey fits"),
        (
            "Reduction in mean squared error (ppt²)",
            "Reduction in weighted squared error",
        ),
    ):
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlabel("Shift from Gaia phase (outer cycles)", fontsize=9)
        ax.axhline(0, color="#686d72", lw=0.8)
        ax.axvline(0, ls=":", color="#747a80", lw=0.8)
        ax.set_xlim(-0.5, 0.5)
        ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
        ax.spines[["right", "top"]].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.text(
        0.085,
        0.04,
        "Dots mark the published Gaia phase. Positive values favor the timing prediction over zero delay.",
        fontsize=9,
        color="#41474b",
    )
    fig.savefig(ROOT / "ground-phase.png", dpi=145, facecolor="white")

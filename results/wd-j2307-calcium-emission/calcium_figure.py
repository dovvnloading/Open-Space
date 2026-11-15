"""Plot the observed calcium triplet, leaving masked velocity bins empty."""
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze import BAND, CENTERS, C, HALF_WIDTH

ROOT = Path(__file__).resolve().parent


def normalized_spectrum(s):
    w, f, iv = s["wave"], s["flux"], s["ivar"]
    use = ((w > BAND[0]) & (w < BAND[1]) & (iv > 0) &
           (s["andmask"] == 0) & np.isfinite(f))
    w, f, iv = w[use], f[use], iv[use]
    v = C * (w[:, None] / CENTERS - 1)
    continuum = np.vander((w - 8580) / 160, 3, increasing=True)
    model = np.column_stack([continuum, np.abs(v) < HALF_WIDTH])
    coefs = np.linalg.lstsq(model * np.sqrt(iv[:, None]), f * np.sqrt(iv), rcond=None)[0]
    baseline = continuum @ coefs[:3]
    return v, f / baseline - 1, iv * baseline**2


def main():
    with gzip.open(ROOT / "spectra.csv.gz", "rt") as stream:
        raw = np.genfromtxt(stream, delimiter=",", names=True, dtype=None, encoding="utf8")
    fig, axes = plt.subplots(2, 3, figsize=(9.3, 5.6), sharex=True, sharey="row")
    edges = np.arange(-990, 991, 90)
    for row, (name, color) in enumerate([("BOSS-2010", "#19676d"), ("DESI-2021", "#853c69")]):
        velocity, flux, weights = normalized_spectrum(raw[raw["epoch"] == name])
        for line, ax in enumerate(axes[row]):
            values = []
            errors = []
            for lo, hi in zip(edges, edges[1:]):
                use = (velocity[:, line] >= lo) & (velocity[:, line] < hi)
                if np.any(use):
                    values.append(np.average(flux[use], weights=weights[use]))
                    errors.append(weights[use].sum()**-0.5)
                else:
                    values.append(np.nan)
                    errors.append(np.nan)
            ax.axvspan(-650, 650, color="#e9edef", zorder=0)
            ax.axhline(0, color="#62747a", lw=0.9)
            ax.errorbar((edges[:-1] + edges[1:]) / 2, values, yerr=errors,
                        fmt="o", ms=3.5, elinewidth=1.0, capsize=0, color=color)
            ax.set_xlim(-990, 990)
            ax.set_xticks([-800, -400, 0, 400, 800])
            ax.spines[["right", "top"]].set_visible(False)
            ax.tick_params(labelsize=8)
            if row == 0:
                ax.set_title(f"Ca II   {CENTERS[line]:.2f} Å", fontsize=11)
        axes[row][0].set_ylabel(name.replace("-", " ") + "\nContinuum-relative excess")
    fig.supxlabel("Velocity relative to the vacuum line center (km/s)", y=0.095, fontsize=11)
    fig.text(0.09, 0.016, "Shading marks the ±650 km/s fit windows. Each point uses unmasked pixels in a 90 km/s bin.\n"
             "Bars show propagated flux errors conditional on the fitted continuum. Empty bins carry no measurement.",
             fontsize=8.5, color="#394e55")
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.23, top=0.90, hspace=0.20, wspace=0.10)
    fig.savefig(ROOT / "calcium-profiles.png", dpi=145, bbox_inches="tight")


if __name__ == "__main__":
    main()

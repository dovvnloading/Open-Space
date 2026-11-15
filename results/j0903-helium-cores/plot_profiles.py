"""Plot the observed line cores and the reproducible fitted centroid separation."""
import csv

import matplotlib.pyplot as plt
import numpy as np

from measure import HERE, LINES, read_pixels


def main():
    pixels = read_pixels()
    with (HERE / "fits.csv").open() as handle:
        result = next(row for row in csv.DictReader(handle)
                      if row["source"] == "XSHOOTER" and row["split"] == "all")
    fig, ax = plt.subplots(figsize=(7, 3.7), layout="constrained")
    for n, (line, color, key) in enumerate(zip(LINES, ["#2a648c", "#a14920"], ["v4922", "v5016"])):
        v, f, w = pixels["XSHOOTER", line]
        keep = (np.abs(v) < 245) & (w.sum(0) > 0)
        wave, flux = v[keep], (f*w).sum(0)[keep]/w.sum(0)[keep]
        err = 1/np.sqrt(w.sum(0)[keep])
        ax.errorbar(wave, flux, yerr=err, color=color, drawstyle="steps-mid",
                    linewidth=1.15, elinewidth=0.7, alpha=0.9, label=f"He I {line:.3f} Å")
        ax.axvline(float(result[key]), color=color, ls="--", lw=1.2)
    ax.set(xlabel="Velocity relative to each air rest wavelength (km/s)",
           ylabel=r"Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)",
           title="SDSS J0903−0133   ·   VLT/X-shooter, 28 November 2024")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.17)
    fig.savefig(HERE/"line-profiles.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()

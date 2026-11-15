"""Render the ground and space photometric phase curves in the supplied period family."""

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from analyze import HERE, clean_ztf, read, subset, fit
from tess_check import extract, get_file, SECTORS


def bins(t, y, w, f, phase, n=12):
    phi = (t * f - phase / (2 * np.pi)) % 1
    ids = np.minimum((phi * n).astype(int), n - 1)
    points = []
    for j in range(n):
        mask = ids == j
        mu = np.average(y[mask], weights=w[mask])
        ww = w[mask]
        _, inv = np.unique(np.floor(t[mask]), return_inverse=True)
        # Inflate the error for residuals sharing a calendar day.
        cluster = np.bincount(inv, weights=ww * (y[mask] - mu))
        se = np.sqrt(sum(cluster**2)) / sum(ww)
        points.append(((j + 0.5) / n, mu, se))
    return np.array(points).T


def plot(cache, output):
    sol = json.loads((HERE / "measurement.json").read_text())
    all_data = clean_ztf(read("ztf.csv.gz"))
    red = subset(all_data, all_data[3] == "zr")
    tlist = []
    flist = []
    wlist = []
    for sector in SECTORS:
        t, cube, base, x, y, *_ = extract(get_file(cache, sector))
        flux = cube[:, y, x]
        detrended = flux - base @ np.linalg.lstsq(base, flux, rcond=None)[0]
        tlist.extend(t)
        flist.extend(detrended)
        wlist.extend(np.repeat(1 / np.var(detrended), len(t)))
    fig, axes = plt.subplots(
        2, 2, figsize=(8.2, 5.9), sharex=True, sharey="row", layout="constrained"
    )
    for col, label in enumerate(["signal", "plus_sidereal"]):
        f = sol["aliases"][label]["frequency"]
        co = fit(*red, f)[1]
        phase = np.arctan2(-co[2], -co[1])
        ax = axes[0, col]
        for band, color, shift in [("zg", "#127969", -0.008), ("zr", "#ab5031", 0.008)]:
            t, y, e, _ = subset(all_data, all_data[3] == band)
            ph, mean, se = bins(t, -1000 * y, 1 / e**2, f, phase)
            ax.errorbar(
                ph + shift,
                mean,
                yerr=se,
                fmt="o",
                ms=4,
                lw=1.1,
                color=color,
                label="ZTF " + band[-1],
            )
        ax.set_title(f"{24/f:.4f} h" + (" preferred" if col == 0 else " daily alias"))
        ax = axes[1, col]
        ph, mean, se = bins(np.array(tlist), np.array(flist), np.array(wlist), f, phase)
        ax.errorbar(
            ph,
            mean,
            yerr=se,
            fmt="o",
            ms=4,
            lw=1.2,
            color="#28578f",
            label="TESS sectors 42, 43, 70",
        )
        ax.set_xlabel("Phase from ZTF r brightness maximum")
    axes[0, 0].set_ylabel("ZTF brightening (mmag)")
    axes[1, 0].set_ylabel("TESS flux residual (electron/s)")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    axes[1, 0].legend(frameon=False, fontsize=8)
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="0.68", lw=0.65)
        ax.grid(axis="y", alpha=0.14)
        ax.set_xlim(0, 1)
    fig.suptitle("WD 0140+130", fontsize=17, fontweight="bold")
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(args.cache, args.output)

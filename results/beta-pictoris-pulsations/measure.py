"""Measure two resolved pulsation peaks of beta Pic in frozen SPOC photometry."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
BASE = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"
MODES = np.array([47.284, 47.439, 50.492, 48.919, 53.692, 54.237, 39.063,
                  46.543, 49.557, 43.527, 34.760, 50.827, 57.454, 38.126,
                  48.136, 44.684, 45.436, 53.007, 53.854, 54.964])


def verified_product(row, cache):
    path = cache / row["product"]
    if not path.exists():
        with urlopen(BASE + row["product"], timeout=120) as response:
            payload = response.read()
        if hashlib.sha256(payload).hexdigest() != row["sha256"]:
            raise ValueError(f"MAST product changed: {path.name}")
        path.write_bytes(payload)
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise ValueError(f"Cached product differs: {path.name}")
    return path


def prepare(data, key):
    time = data["TIME"]
    flux = data[key]
    keep = (data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(flux)
    keep &= np.isfinite(data["PDCSAP_FLUX"])
    time = time[keep]
    flux = flux[keep].astype(float)
    flux = 1e6 * (flux / np.median(flux) - 1)
    segments = np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > 0.05) + 1)
    # A local linear baseline per half day removes low-frequency variability.
    # A segment break prevents interpolation or slope fits across spacecraft gaps.
    for seg in segments:
        blocks = np.floor((time[seg] - time[seg[0]]) * 2)
        for indices in np.split(seg, np.flatnonzero(np.diff(blocks)) + 1):
            dt = time[indices] - np.mean(time[indices])
            flux[indices] -= np.mean(flux[indices])
            denom = np.dot(dt, dt)
            if denom:
                flux[indices] -= dt * np.dot(dt, flux[indices]) / denom
    return time, flux


def design(time, frequencies):
    angles = 2 * np.pi * time[:, None] * frequencies
    return np.column_stack([np.cos(angles), np.sin(angles)])


def coefficients(x, y):
    return np.linalg.lstsq(x, y, rcond=None)[0]


def amplitudes(time, flux, refine=False):
    time = time - np.mean(time)
    freqs = MODES.copy()
    if refine:
        # Remove the nuisance-mode fit once, then allow each local peak to move
        # by ±0.025 cycles/day. Bounds keep the two peak identifications disjoint.
        x = design(time, freqs)
        beta = coefficients(x, flux)
        idx = [0, 1, len(freqs), len(freqs) + 1]
        y_pair = flux - x @ beta + x[:, idx] @ beta[idx]

        def objective(pair):
            pair_x = design(time, pair)
            resid = y_pair - pair_x @ coefficients(pair_x, y_pair)
            return np.mean(resid ** 2)

        optimum = minimize(objective, freqs[:2], method="Nelder-Mead",
                           bounds=[(f - .025, f + .025) for f in freqs[:2]],
                           options={"xatol": 1e-8, "initial_simplex": [
                               freqs[:2], freqs[:2] + [0.001, 0],
                               freqs[:2] + [0, 0.001]]})
        if not optimum.success:
            raise RuntimeError(optimum.message)
        freqs[:2] = optimum.x
    n = len(freqs)
    b = coefficients(design(time, freqs), flux)
    return freqs[:2], np.hypot(b[:2], b[n:n + 2])


def measure(path, sector, key):
    data = fits.getdata(path, 1)
    t, y = prepare(data, key)
    pair, amp = amplitudes(t, y)
    free, amp_free = amplitudes(t, y, refine=True)
    midpoint = (t.min() + t.max()) / 2
    _, first = amplitudes(t[t <= midpoint], y[t <= midpoint])
    _, last = amplitudes(t[t > midpoint], y[t > midpoint])
    # Published close components can beat. Compare actual modeled band power on
    # observed cadences with the simultaneous 47.44/d fit.
    doublet = np.r_[47.2686, 47.2853, MODES[1:]]
    dx = design(t - np.mean(t), doublet)
    db = coefficients(dx, y)
    ni = len(doublet)
    band = dx[:, [0, 1, ni, ni+1]] @ db[[0, 1, ni, ni+1]]
    high = dx[:, [2, ni+2]] @ db[[2, ni+2]]
    band_ratio = np.sqrt(np.mean(band**2) / np.mean(high**2))
    return dict(doublet_band_ratio=band_ratio, sector=sector, flux=key, n=len(t), btjd_start=t.min(), btjd_end=t.max(),
                low_ppm=amp[0], high_ppm=amp[1], ratio=amp[0]/amp[1],
                free_low_freq=free[0], free_high_freq=free[1],
                free_low_ppm=amp_free[0], free_high_ppm=amp_free[1],
                first_half_ratio=first[0]/first[1], second_half_ratio=last[0]/last[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "measurements.csv")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (ROOT / "manifest.csv").open() as handle:
        manifest = list(csv.DictReader(handle))
    rows = []
    for entry in manifest:
        path = verified_product(entry, args.cache)
        rows.extend(measure(path, int(entry["sector"]), key)
                    for key in ("SAP_FLUX", "PDCSAP_FLUX"))
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) else v)
                             for k, v in row.items()})


if __name__ == "__main__":
    main()

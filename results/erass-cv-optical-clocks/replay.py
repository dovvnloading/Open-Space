"""Remeasure archived eRASS CV light curves and variable-flux centroids from MAST."""
import argparse
import csv
import hashlib
import json
import tempfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
from astropy.wcs import WCS
from scipy.ndimage import median_filter
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
MAST = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"


def product(row, cache):
    path = cache / row["product"]
    if not path.exists():
        with urlopen(MAST + row["product"], timeout=150) as response:
            path.write_bytes(response.read())
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != row["sha256"]:
        raise ValueError(f"Archived bytes differ for {path.name}: {digest}")
    return path


def cadence_series(path):
    with fits.open(path) as hdu:
        data = hdu[1].data
        time = data["TIME"].astype(float)
        flux = data["PDCSAP_FLUX"].astype(float)
        keep = (data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(flux)
        time, flux = time[keep], flux[keep]
    residual = flux - median_filter(flux, size=251)
    scale = np.median(np.abs(residual - np.median(residual))) * 1.4826
    return time, np.clip(residual, -5 * scale, 5 * scale)


def measure_clock(path):
    time, flux = cadence_series(path)
    grid = np.linspace(2, 200, round(np.ptp(time) * 198 * 6))
    model = LombScargle(time, flux)
    power = model.power(grid, method="fast", assume_regular_frequency=True)
    peaks = find_peaks(power)[0]
    idx = peaks[np.argmax(power[peaks])]
    peak = grid[idx]
    time -= time[0]

    def refine(t, y, center):
        def cost(freq):
            phase = 2 * np.pi * freq * t
            design = np.column_stack((np.sin(phase), np.sin(2 * phase),
                                      np.cos(phase), np.cos(2 * phase)))
            fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
            return -np.dot(fitted, fitted)
        return minimize_scalar(cost, bounds=(center - .022, center + .022),
                               method="bounded", options={"xatol": 1e-9}).x

    frequency = refine(time, flux, peak)
    days = time.astype(int)
    jack = []
    for block in np.array_split(np.unique(days), 8):
        keep = ~np.isin(days, block)
        jack.append(refine(time[keep], flux[keep], frequency))
    se = (len(jack) - 1) ** .5 * np.std(jack)
    return peak, dict(n=len(time), period_min=1440 / frequency,
                      jackknife_se_min=1440 * se / frequency**2, peak_power=power[idx])


def localize(path, frequency, coords):
    with fits.open(path) as hdu:
        data = hdu[1].data
        valid = (data["QUALITY"] == 0) & np.isfinite(data["TIME"])
        time = data["TIME"][valid].astype(float)
        cube = data["FLUX"][valid].astype(float)
        wcs = WCS(hdu[2].header, naxis=2)
    target = np.array(wcs.world_to_pixel_values(coords["ra"], coords["dec"]))
    height, width = cube.shape[1:]
    yy, xx = np.indices((height, width))
    flux = cube.reshape(len(time), -1)
    usable = np.isfinite(flux).all(axis=0)
    flux = flux[:, usable]
    phase = 2 * np.pi * frequency * time
    design = np.column_stack((np.sin(phase), np.cos(phase)))
    blocks = (2 * time).astype(int)
    for block in np.unique(blocks):
        use = blocks == block
        t = time[use]
        baseline = np.column_stack((np.ones(len(t)), t - np.mean(t)))
        for matrix in (design, flux):
            matrix[use] -= baseline @ np.linalg.lstsq(baseline, matrix[use], rcond=None)[0]
    amplitudes = np.linalg.lstsq(design, flux, rcond=None)[0]
    length = np.linalg.norm(amplitudes, axis=0)
    peak = np.argmax(length)
    signed = (amplitudes[:, peak] / length[peak]) @ amplitudes
    yinit, xinit = divmod(np.flatnonzero(usable)[peak], width)
    x, y = xx.ravel()[usable], yy.ravel()[usable]

    def residual(c):
        amp, x0, y0, log_sx, log_sy, angle, floor = c
        dx, dy = x - x0, y - y0
        u = (np.cos(angle) * dx + np.sin(angle) * dy) / np.exp(log_sx)
        v = (-np.sin(angle) * dx + np.cos(angle) * dy) / np.exp(log_sy)
        return floor + amp * np.exp(-.5 * (u*u + v*v)) - signed

    fit = least_squares(residual, [length[peak], xinit, yinit, 0, 0, 0, 0],
                        bounds=([0, -2, -2, -1.3, -1.3, -10, -np.inf],
                                [np.inf, width + 2, height + 2, 1.6, 1.6, 10, np.inf]))
    center = fit.x[1:3]
    return dict(offset_px=np.linalg.norm(center - target), source_xy=target.tolist(),
                modulation_xy=center.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", choices=("156", "110", "116", "106", "21"), default="156")
    parser.add_argument("--cache", type=Path)
    args = parser.parse_args()
    targets = json.loads((HERE / "catalog-checks.json").read_text())
    with (HERE / "products.csv").open() as stream:
        rows = [row for row in csv.DictReader(stream) if row["id"] == args.id]
    with tempfile.TemporaryDirectory() as temporary:
        cache = args.cache or Path(temporary)
        cache.mkdir(exist_ok=True, parents=True)
        for row in rows:
            if not row["product"].endswith("_lc.fits"):
                continue
            freq, metrics = measure_clock(product(row, cache))
            pixels = next((r for r in rows if r["sector"] == row["sector"]
                           and r["product"].endswith("_tp.fits")), None)
            if pixels:
                metrics.update(localize(product(pixels, cache), freq, targets[args.id]))
            print(json.dumps(dict(target=args.id, sector=row["sector"], **metrics)))


if __name__ == "__main__":
    main()

"""Replay the independent 2026 comet recovery search from the exposure manifest."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates, median_filter
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parent
Y, X = np.mgrid[-26:26.1:0.5, -26:26.1:0.5]
SEARCH = (abs(X) <= 11) & (abs(Y) <= 2)
CONTROL = (abs(X) > 10) | (abs(Y) > 7)


def normalized(image, mask=None):
    sample = image if mask is None else image[mask]
    center = np.nanmedian(sample)
    scale = 1.4826 * np.nanmedian(abs(sample - center))
    return (image - center) / scale


def project(row, cache):
    path = cache / (row["exposure"] + ".fits")
    if not path.exists():
        with urlopen(row["cutout_url"], timeout=90) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != row["sha256"]:
            raise ValueError(f"Archive checksum mismatch for {path.name}")
        path.write_bytes(content)
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise ValueError(f"Cached checksum mismatch for {path.name}")
    with fits.open(path) as h:
        header = h[0].header
        image = h[0].data.astype(float)
    image -= median_filter(image, size=15)
    y, x = np.mgrid[-6:7, -6:7]
    sigma = float(header["SEEING"]) / 2.38
    kernel = np.exp(-(x*x + y*y)/(2*sigma*sigma))
    smoothed = normalized(fftconvolve(image, kernel, mode="same"))
    ra, dec = float(row["comet_ra_deg"]), float(row["comet_dec_deg"])
    px, py = WCS(header).all_world2pix(ra + X/3600/np.cos(np.radians(dec)), dec + Y/3600, 0)
    return map_coordinates(smoothed, [py, px], order=1, mode="constant", cval=np.nan)


def stack(cube):
    return normalized(np.nanmean(np.clip(cube, -3, 3), axis=0), CONTROL)


def peak(stack_image):
    return np.unravel_index(np.nanargmax(np.where(SEARCH, stack_image, np.nan)), stack_image.shape)


def summarize(image, halves, count):
    split = []
    for i, first in enumerate(halves):
        idx = peak(first)
        split.append({"east_north_arcsec": [float(X[idx]), float(Y[idx])],
                      "selection_snr": round(float(first[idx]), 2),
                      "other_half_snr": round(float(halves[1-i][idx]), 2)})
    return {"frames": count, "nominal_snr": round(float(image[52, 52]), 2),
            "search_peak_snr": round(float(image[peak(image)]), 2),
            "control_peak_snr": round(float(np.nanmax(image[CONTROL])), 2),
            "chronological_halves": split}


def measure(cube):
    mid = len(cube) // 2
    return summarize(stack(cube), [stack(cube[:mid]), stack(cube[mid:])], len(cube))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, help="Rebuild the stacks from source images")
    args = parser.parse_args()
    with (ROOT / "recovery_frames.csv").open() as handle:
        rows = sorted(csv.DictReader(handle), key=lambda r: r["exposure"])
    deeper = np.array([float(row["maglim"]) >= 20.2 for row in rows])
    result = {"nights": len({r["exposure"][:8] for r in rows})}
    if args.cache:
        args.cache.mkdir(exist_ok=True, parents=True)
        images = np.array([project(row, args.cache) for row in rows])
        result.update(all=measure(images), deeper=measure(images[deeper]))
    else:
        # Normalized maps are stored in units of one hundredth of the
        # empirical spatial noise, the precision of the released summary.
        with np.load(ROOT / "recovery_stacks.npz", allow_pickle=False) as maps:
            for label, count in [("all", len(rows)), ("deeper", int(deeper.sum()))]:
                image, first, last = (maps[f"{label}_{part}"] / 100
                                      for part in ["full", "first", "last"])
                result[label] = summarize(image, [first, last], count)
    print(json.dumps(result, indent=2))

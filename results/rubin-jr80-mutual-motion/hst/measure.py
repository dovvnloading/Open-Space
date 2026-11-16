"""Reproduce June 13 HST relative astrometry from four frozen science cutouts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent


def download_cutout(obs):
    """Optional full-file verification and extraction. Requires astropy."""
    from tempfile import TemporaryFile
    from astropy.io import fits

    sha = hashlib.sha256()
    with urlopen(obs["url"], timeout=120) as response, TemporaryFile() as temp:
        for block in iter(lambda: response.read(4 * 1024 * 1024), b""):
            sha.update(block)
            temp.write(block)
        if sha.hexdigest() != obs["sha256_fits"]:
            raise ValueError("MAST product revision or download mismatch")
        temp.seek(0)
        with fits.open(temp, memmap=False) as hdus:
            x, y = obs["origin_zero_based"]
            return hdus["SCI"].data[y:y+73, x:x+73].astype(">f4")


def centroids(image, seed):
    x0, y0 = seed
    y, x = np.mgrid[y0-8:y0+9, x0-8:x0+9]
    stamp = image[y, x].astype(float)
    radius = np.hypot(x-x0, y-y0)
    background = np.median(stamp[radius >= 7])
    # Retain core pixels, allowing a separate fitted width on each image axis.
    core = radius < 3.5
    def residual(par):
        b, amp, xc, yc, sx, sy = par
        model = b + amp * np.exp(-0.5 * (((x-xc)/sx)**2 + ((y-yc)/sy)**2))
        return (model - stamp)[core]

    model = least_squares(residual, [background, stamp.max(), x0, y0, 1, 1],
                          bounds=([-2, 0, x0-2, y0-2, 0.4, 0.4],
                                  [3, 200, x0+2, y0+2, 3, 3]))
    gaussian = model.x[2:4]
    pos = np.array(seed, dtype=float)
    for _ in range(4):
        select = np.hypot(x-pos[0], y-pos[1]) <= 3.5
        weight = (stamp-background)*select
        pos = np.array([np.sum(weight*x), np.sum(weight*y)]) / np.sum(weight)
    return gaussian, pos, weight.sum()


def measure(obs, image):
    if hashlib.sha256(image.astype(">f4").tobytes()).hexdigest() != obs["sha256_pixels"]:
        raise ValueError("Science pixels mismatch")
    first, second = [centroids(image, seed) for seed in obs["seeds_zero_based_cutout"]]
    cd = np.asarray(obs["cd_deg_per_pixel"])
    result = {}
    for method, a, b in zip(("core", "aperture"), first[:2], second[:2]):
        east, north = 3600 * cd @ (b-a)
        result[method + "_east_arcsec"] = east
        result[method + "_north_arcsec"] = north
        result[method + "_separation_arcsec"] = np.hypot(east, north)
        result[method + "_pa_deg"] = np.degrees(np.arctan2(east, north)) % 360
    result["aperture_delta_mag"] = -2.5*np.log10(second[2]/first[2])
    return {k: f"{v:.6f}" for k, v in result.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", action="store_true", help="download and verify original DRC FITS")
    args = parser.parse_args()
    import sys
    includes = json.loads((ROOT/"manifest.json").read_text())
    with np.load(ROOT/"science_cutouts.npz") as cutouts:
        writer = None
        for obs in includes:
            image = download_cutout(obs) if args.archive else cutouts[obs["rootname"]]
            row = dict(rootname=obs["rootname"], midpoint_mjd_utc=obs["midpoint_mjd_utc"],
                       **measure(obs, image))
            if writer is None:
                writer = csv.DictWriter(sys.stdout, fieldnames=row, lineterminator="\n")
                writer.writeheader()
            writer.writerow(row)

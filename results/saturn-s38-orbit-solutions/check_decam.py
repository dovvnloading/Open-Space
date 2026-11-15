"""Fixed-position aperture measurements on two independent NOIRLab DECam images."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent


def aperture(image, x, y):
    ix, iy = int(x), int(y)
    stamp = image[iy-10:iy+11, ix-10:ix+11]
    yy, xx = np.mgrid[-10:11, -10:11]
    r2 = (xx + ix - x)**2 + (yy + iy - y)**2
    return np.sum(stamp[r2 <= 3.7**2] - np.median(stamp[(r2 >= 7**2) & (r2 < 10**2)]))


def measure(image, x, y):
    ix, iy = int(x)-105, int(y)-105
    patch = image[iy:iy+211, ix:ix+211].astype(float)
    grid = np.arange(12, 197, 8.5)
    references = [aperture(patch, a, b) for a in grid for b in grid]
    median = np.median(references)
    scatter = 1.4826 * np.median(abs(references - median))
    flux = aperture(patch, x-ix, y-iy)
    return dict(x=float(x), y=float(y), comparison_centers=len(references),
                aperture_counts=flux, median_counts=median, mad_sigma_counts=scatter,
                snr=(flux-median)/scatter)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, help="Directory containing the FITS files, named by archive_filename basename")
    args = p.parse_args()
    results = []
    for row in json.loads((ROOT / "decam_positions.json").read_text()):
        path = args.directory / Path(row["archive_filename"]).name
        with path.open("rb") as f:
            if hashlib.file_digest(f, "md5").hexdigest() != row["md5sum"]:
                raise ValueError(f"MD5 disagrees with archive for {path}")
        with fits.open(path) as hdus:
            for hdu in hdus[1:]:
                size = np.array([hdu.header["NAXIS1"], hdu.header["NAXIS2"]])
                coords = [point[:2] for point in row["positions"]]
                xy = WCS(hdu.header).all_world2pix(coords, 0, quiet=True)
                for branch, pos in zip(("moderate", "extreme"), xy):
                    if np.all(pos > 106) and np.all(pos < size-106):
                        results.append(dict(file=path.name, detector=hdu.name, branch=branch,
                                            **measure(hdu.data, *pos)))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

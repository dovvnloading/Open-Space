"""Measure the fixed Ca II triplet windows in each native BOSS red exposure."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

URL = ("https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/"
       "spectra/full/4209/spec-4209-55478-0518.fits")
SHA256 = "0829176abcc18a89b9ddb7e372622078c5edb9fde3e4007dd30c1b3d7428b303"
LINES = np.array([8500.36, 8544.44, 8664.52])


def measure(hdu):
    data = hdu.data
    wave = 10.0 ** data["loglam"].astype(float)
    flux = data["flux"]
    weight = data["ivar"]
    clean = ((wave > 8425) & (wave < 8735) & (weight > 0) &
             (data["mask"] == 0) & np.isfinite(flux))
    wave, flux, weight = wave[clean], flux[clean], weight[clean]
    coordinate = (wave - 8580) / 160
    line_windows = np.abs((wave[:, None] / LINES - 1) * 299792.458) < 650
    continuum = np.column_stack([np.ones_like(wave), coordinate, coordinate**2])
    design = np.column_stack([continuum, line_windows])
    normal = (design.T * weight) @ design
    covariance = np.linalg.inv(normal)
    amplitudes = np.linalg.solve(normal, (design.T * weight) @ flux)[-3:]
    line_cov = covariance[-3:, -3:]
    precision = np.linalg.inv(line_cov)
    common_z = (precision @ amplitudes).sum() / np.sqrt(precision.sum())
    individual_z = amplitudes / np.sqrt(line_cov.diagonal())
    return [hdu.name.split("-")[1], hdu.header["EXPTIME"], int(clean.sum()),
            *line_windows.sum(axis=0), common_z, *individual_z]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, default=Path("spec-4209-full.fits"))
    parser.add_argument("--output", type=Path, default=Path("boss-exposures.csv"))
    args = parser.parse_args()
    if not args.fits.exists():
        with urlopen(URL, timeout=90) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != SHA256:
            raise ValueError("SDSS full spectrum differs from the frozen release")
        args.fits.write_bytes(raw)
    if hashlib.sha256(args.fits.read_bytes()).hexdigest() != SHA256:
        raise ValueError("Unexpected BOSS full-spectrum bytes")
    with fits.open(args.fits) as hdus:
        result = [measure(hdu) for hdu in hdus if hdu.name.startswith("R2-")]
    if len(result) != 6:
        raise ValueError("Expected six native red exposures")
    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["exposure_id", "duration_s", "clean_pixels", "pixels_8500",
                         "pixels_8544", "pixels_8665", "triplet_z", "z_8500",
                         "z_8544", "z_8665"])
        writer.writerows(result)


if __name__ == "__main__":
    main()

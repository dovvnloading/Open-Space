"""Replay line widths or re-extract SPRAT SPEC_SS and a local-sky control."""
import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
C_KMS = 299792.458
LINES = {"Halpha": 6562.80, "Hbeta": 4861.33}
HALF_WIDTHS = (60, 80, 100)


def extract(cache, output):
    from astropy.io import fits

    spectra = {}
    wavelength = None
    with (HERE / "sprat-observations.csv").open() as stream:
        for row in csv.DictReader(stream):
            name = row["product"]
            contents = gzip.decompress((cache / name).read_bytes())
            if hashlib.sha256(contents).hexdigest() != row["fits_sha256"]:
                raise ValueError(f"Uncompressed FITS digest differs for {name}")
            with fits.open(io.BytesIO(contents)) as hdus:
                hdu = hdus["SPEC_SS"]
                header = hdu.header
                counts = hdu.data.ravel()
                pixels = np.arange(1, len(counts) + 1)
                grid = header["CRVAL1"] + (pixels - header["CRPIX1"]) * header["CDELT1"]
                if wavelength is not None and not np.array_equal(wavelength, grid):
                    raise ValueError("Exposure grids differ")
                wavelength = grid
                visit = row["exposure"]
                spectra[f"spec_ss_{visit}"] = counts
                # Fixed zero-based rows in the rectified long slit. A 9-row
                # extraction uses the median of 18 flanking rows as local sky.
                image = hdus["LSS_NONSS"].data
                spectra[f"slit_sky_{visit}"] = np.median(image[20:185], axis=0)
                sky = np.median(np.concatenate((image[75:84], image[103:112])), axis=0)
                spectra[f"local_sky_{visit}"] = image[89:98].sum(axis=0) - 9 * sky
    columns = {"wavelength_angstrom": wavelength, **spectra}
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        for values in zip(*columns.values()):
            writer.writerow([format(float(v), ".9g") for v in values])


def profile(x, baseline, slope, amplitude, center, sigma):
    return baseline + slope * x + amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def measure(wavelength, counts, rest, half_width):
    offset = wavelength - rest
    keep = np.abs(offset) <= half_width
    x, y = offset[keep], counts[keep]
    baseline = np.median(y[np.abs(x) > half_width * 0.6])
    fit, _ = curve_fit(profile, x, y, p0=[baseline, 0, max(y) - baseline, -4, 14],
                       bounds=([-np.inf, -np.inf, 0, -25, 2], [np.inf, np.inf, np.inf, 25, 40]))
    fwhm = np.sqrt(8 * np.log(2)) * fit[4]
    return rest + fit[3], fwhm, C_KMS * fwhm / rest


def replay():
    data = np.genfromtxt(HERE / "sprat-spectra.csv", names=True, delimiter=",")
    w = data["wavelength_angstrom"]
    print("line,extraction,exposure,fit_half_width_A,center_A,FWHM_A,FWHM_km_s")
    for line, rest in LINES.items():
        for name in data.dtype.names[1:]:
            method, exposure = name.rsplit("_", 1)
            if method == "slit_sky":
                continue
            for span in HALF_WIDTHS:
                center, width, velocity = measure(w, data[name], rest, span)
                print(f"{line},{method},{exposure},{span},{center:.2f},{width:.2f},{velocity:.1f}")

    for exposure in (1, 2, 3):
        sky = data[f"slit_sky_{exposure}"]
        for label, rest in (("OI5577", 5577.34), ("OI6300", 6300.30)):
            center, width, velocity = measure(w, sky, rest, 42)
            print(f"{label},slit_sky,{exposure},42,{center:.2f},{width:.2f},{velocity:.1f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, metavar="FITS_DIRECTORY")
    parser.add_argument("--output", type=Path, help="Regenerated spectrum CSV")
    args = parser.parse_args()
    if args.extract:
        if args.output is None:
            parser.error("--extract needs --output")
        extract(args.extract, args.output)
    else:
        replay()

"""Retrieve specified archive pixels and repeat a measurement. NumPy, SciPy, Astropy required."""

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from measure import measure

HERE = Path(__file__).resolve().parent


def download(url, path):
    with urlopen(url, timeout=300) as response, path.open("wb") as dst:
        shutil.copyfileobj(response, dst)
    return path


def file_checksum(path, algorithm):
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def reduce_decam(image, target_ra, target_dec, output):
    with fits.open(image) as hdus:
        for extension in hdus[1:]:
            wcs = WCS(extension.header)
            x, y = wcs.all_world2pix([[target_ra, target_dec]], 0)[0]
            ny, nx = extension.shape
            if 12 < x < nx - 12 and 12 < y < ny - 12:
                cut = Cutout2D(extension.data, (x, y), (590, 590), wcs=wcs, mode="trim")
                fits.PrimaryHDU(cut.data, header=cut.wcs.to_header(relax=True)).writeto(
                    output, overwrite=True
                )
                return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("image_id", nargs="+", help="IDs from inputs.csv, e.g. 2497047p")
    ap.add_argument("--cache", type=Path, default=Path("archive-pixels"))
    args = ap.parse_args()
    rows = {row["id"]: row for row in csv.DictReader((HERE / "inputs.csv").open())}
    args.cache.mkdir(exist_ok=True, parents=True)
    for key in args.image_id:
        r = rows[key]
        local = args.cache / (key + ".fits")
        if r["survey"] == "CFHT":
            if not local.exists():
                download(r["url"], local)
            reference = local
        else:
            reference = args.cache / (key + ".fits.fz")
            if not reference.exists():
                download(r["url"], reference)
        digest = file_checksum(reference, r["checksum_algorithm"])
        if digest != r["checksum"]:
            raise ValueError(f"Checksum mismatch for {key}: {digest}")
        if r["survey"] == "DECam" and not reduce_decam(
            reference, float(r["ra"]), float(r["dec"]), local
        ):
            print(key, "position outside illuminated detectors")
            continue
        metadata = {
            "MJD": r["mjd"],
            "Object_RA": r["ra"],
            "Object_Dec": r["dec"],
            "Exptime": r["exptime"],
            "Filter": r["band"],
        }
        local.with_suffix(".json").write_text(json.dumps(metadata) + "\n")
        print(json.dumps(measure(local), indent=2))


if __name__ == "__main__":
    main()

"""Redownload the four cutouts, verify their scientific arrays, and extract light curves."""
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def digest(hdus):
    h = hashlib.sha256()
    for name in ("TIME", "TIMECORR", "QUALITY", "FLUX"):
        a = hdus[1].data[name]
        h.update(np.ascontiguousarray(a.astype(a.dtype.newbyteorder("<"))).tobytes())
    for key in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2",
                "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CTYPE1", "CTYPE2"):
        h.update(str(hdus[2].header[key]).encode() + b"\n")
    return h.hexdigest()


def extract(hdus, sector):
    table = hdus[1].data
    good = ((table["QUALITY"] == 0) & np.isfinite(table["TIME"]) &
            np.isfinite(table["TIMECORR"]) & np.isfinite(table["FLUX"]).all(axis=(1, 2)))
    table = table[good]
    cube = table["FLUX"].astype(float)
    spacecraft = table["TIME"].astype(float) - table["TIMECORR"].astype(float)
    v = pd.read_csv(DATA / f"vectors-{sector}.csv")
    if spacecraft.min() < v.jd.min() - 2457000 or spacecraft.max() > v.jd.max() - 2457000:
        raise ValueError("Vector ephemeris does not span the observation")
    xyz = CubicHermiteSpline(v.jd - 2457000, v[["x", "y", "z"]], v[["vx", "vy", "vz"]])(spacecraft)
    ra, dec = np.deg2rad([74.68197144, -10.7558383])
    direction = [np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)]
    bary = spacecraft + xyz @ direction * 499.00478384 / 86400
    x, y = WCS(hdus[2].header).world_to_pixel_values(np.rad2deg(ra), np.rad2deg(dec))
    x, y = int(np.round(x)), int(np.round(y))
    background = np.percentile(cube, 15, axis=(1, 2))
    box = np.sum(cube[:, y-1:y+2, x-1:x+2], axis=(1, 2))
    data = pd.DataFrame({"time": bary, "central": cube[:, y, x] - background,
                         "box": box - 9 * background})
    data.to_csv(DATA / f"s{sector}.csv.gz", index=False, float_format="%.10g",
                compression={"method": "gzip", "mtime": 0})


if __name__ == "__main__":
    manifest = json.loads((ROOT / "sources.json").read_text())
    for row in manifest["tess"]:
        with urlopen(row["cutout"], timeout=240) as response:
            payload = response.read()
        if payload.startswith(b"PK"):
            with ZipFile(io.BytesIO(payload)) as archive:
                fits_name, = (n for n in archive.namelist() if n.endswith(".fits"))
                payload = archive.read(fits_name)
        with fits.open(io.BytesIO(payload)) as hdus:
            if digest(hdus) != row["arrays_sha256"]:
                raise ValueError(f"Scientific inputs changed for sector {row['sector']}")
            extract(hdus, row["sector"])

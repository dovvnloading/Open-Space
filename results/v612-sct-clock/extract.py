"""Download a 5 by 5 TESSCut (or reuse a supplied cutout) and reduce the pixels."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import median_abs_deviation

ROOT = Path(__file__).resolve().parent
RA, DEC = 277.94118, -14.31554
GAIA = [277.941193, -14.315539]


def reduce_sector(sector, local=None):
    params = dict(ra=RA, dec=DEC, y=5, x=5, units="px", sector=sector)
    url = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + urlencode(params)
    if local:
        content = Path(local).read_bytes()
    else:
        with urlopen(url, timeout=180) as response, ZipFile(io.BytesIO(response.read())) as archive:
            content = archive.read(next(n for n in archive.namelist() if n.endswith(".fits")))
    with fits.open(io.BytesIO(content)) as hdul:
        data = hdul[1].data
        xy = WCS(hdul[2].header).all_world2pix([GAIA], 0)[0]
        mask = (data["QUALITY"] == 0) & np.isfinite(data["TIME"]) & np.isfinite(data["FLUX"]).all(axis=(1, 2))
        time = data["TIME"][mask].astype(float)
        pixels = data["FLUX"][mask].reshape(-1, 25).astype(float)
        names = [str(v) for v in data["FFI_FILE"][mask]]
    frame = pd.DataFrame(pixels, index=pd.to_datetime(time, unit="D"))
    residual = pixels - frame.rolling("1D", center=True, min_periods=1).median().to_numpy()
    y, x = np.indices((5, 5))
    radius = np.hypot(x - xy[0], y - xy[1]).ravel()
    outer = residual[:, radius >= 2.1]
    scaled = outer / np.maximum(median_abs_deviation(outer, axis=0), .01)
    vectors = np.linalg.svd(scaled, full_matrices=False)[0][:, :2]
    corrected = residual - vectors @ (vectors.T @ residual)
    aperture = radius < 1.15
    light = corrected[:, aperture].sum(axis=1)
    raw_light = residual[:, aperture].sum(axis=1)
    keep = np.abs(light - np.median(light)) < 6 * median_abs_deviation(light)
    time, corrected, light, raw_light = time[keep], corrected[keep], light[keep], raw_light[keep]
    # Projected and unprojected series have very different outlier distributions.
    # Keep all projected-good cadences, recording the independent raw-aperture mask.
    raw_keep = np.abs(raw_light - np.median(raw_light)) < 6 * median_abs_deviation(raw_light)
    table = pd.DataFrame({"btjd": time, "corrected_e_per_s": light,
                          "unprojected_e_per_s": raw_light, "unprojected_good": raw_keep.astype(int)})
    table.to_csv(ROOT / f"data/tess-{sector}.csv.gz", index=False, float_format="%.8f",
                 compression={"method":"gzip","mtime":0})
    f = 1 / .15956199
    angle = 2 * np.pi * (time + 2457000 - 2459500) * f
    design = np.column_stack((np.ones(len(time)), np.cos(angle), np.sin(angle)))
    coeff = np.linalg.lstsq(design, corrected, rcond=None)[0]
    return dict(sector=sector, url=url, fits_sha256=hashlib.sha256(content).hexdigest(),
                first_ffi=names[0], last_ffi=names[-1], good_pixels_cadences=len(keep),
                clipped_cadences=int((~keep).sum()), target_xy=xy.tolist(),
                aperture_xy=[[int(i%5),int(i//5)] for i in np.flatnonzero(aperture)],
                pixel_semiamplitude_e_per_s=np.hypot(coeff[1],coeff[2]).reshape(5,5).round(5).tolist())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits-dir", type=Path, help="Directory containing tess80.fits and tess92.fits")
    args = parser.parse_args()
    rows = [reduce_sector(s, args.fits_dir / f"tess{s}.fits" if args.fits_dir else None) for s in (80, 92)]
    (ROOT / "data/tess-manifest.json").write_text(json.dumps(rows, indent=2) + "\n")

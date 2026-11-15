"""Repeat the early-sector aperture test on MAST full-frame-image cutouts."""
import argparse
import io
import json
import zipfile
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.wcs import WCS
from scipy.ndimage import median_filter


HERE = Path(__file__).resolve().parent
FREQUENCY = {"110": 11.323, "156": 16.385}
SECTORS = {"110": [11, 38], "156": [2, 29]}


def measure(path, coords, frequency):
    with fits.open(path) as hdu:
        rows = hdu[1].data
        keep = (rows["QUALITY"] == 0) & np.isfinite(rows["TIME"])
        time = rows["TIME"][keep].astype(float)
        cube = rows["FLUX"][keep].astype(float)
        x0, y0 = WCS(hdu[2].header, naxis=2).world_to_pixel_values(coords["ra"], coords["dec"])
    y, x = np.indices(cube.shape[1:])
    dist = np.hypot(x - x0, y - y0)
    aperture = dist < 1.6
    flux = cube[:, aperture].sum(axis=1) - aperture.sum() * np.median(cube[:, dist >= 2.2], axis=1)
    keep = np.isfinite(flux)
    time, flux = time[keep], flux[keep]
    step = np.median(np.diff(time))
    width = max(5, round(.6 / step) // 2 * 2 + 1)
    flux -= median_filter(flux, size=width)
    scale = 1.4826 * np.median(np.abs(flux - np.median(flux)))
    flux = np.clip(flux, -4 * scale, 4 * scale)
    freq = np.arange(2, .49 / step, 1 / np.ptp(time) / 6)
    model = LombScargle(time, flux)
    power = model.power(freq, method="fast", assume_regular_frequency=True)
    best = np.argmax(power)
    return dict(n=len(time), cadence_min=1440 * step, first_bjd_tdb=2457000+time[0],
                last_bjd_tdb=2457000+time[-1], first_date_tdb=Time(2457000+time[0], format="jd", scale="tdb").isot[:10],
                strongest_cpd=freq[best], strongest_power=power[best],
                strongest_period_min=1440/freq[best], power_at_prediction=model.power(frequency))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", choices=("110", "156"), default="156")
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    coords = json.loads((HERE / "catalog-checks.json").read_text())[args.id]
    for sector in SECTORS[args.id]:
        path = args.cache / f"cut{args.id}_{sector}.fits"
        query = urlencode(dict(ra=coords["ra"], dec=coords["dec"], sector=sector, x=5, y=5, units="px"))
        url = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + query
        if not path.exists():
            with urlopen(url, timeout=180) as response:
                with zipfile.ZipFile(io.BytesIO(response.read())) as zipped:
                    name = next(n for n in zipped.namelist() if n.endswith(".fits"))
                    path.write_bytes(zipped.read(name))
        values = measure(path, coords, FREQUENCY[args.id])
        print(json.dumps(dict(id=args.id, sector=sector, url=url, **values)))


if __name__ == "__main__":
    main()

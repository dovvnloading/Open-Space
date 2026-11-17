from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import warnings, pandas as pd, numpy as np
from pathlib import Path
from io import BytesIO
from math import prod, ceil
from get_stamps import section
from concurrent.futures import ThreadPoolExecutor

warnings.simplefilter("ignore", FITSFixedWarning)
site = pd.read_csv("onsets.csv").query("site==1004").iloc[0]


def patch(url):
    with section(url, 0, 57599) as src:
        b = BytesIO(src.read())
    p = fits.Header.fromfile(b)
    h = fits.Header.fromfile(b)
    off = b.tell()
    w = WCS(h, naxis=2)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    x, y = w.all_world2pix(site.ra, site.dec, 0)
    xi, yi = int(np.rint(x)), int(np.rint(y))
    r = 10

    def readpart(offset, dtype):
        lo = yi - r
        hi = yi + r + 1
        with section(url, offset + lo * nx * 4, offset + hi * nx * 4 - 1) as src:
            return (
                np.frombuffer(src.read(), dtype)
                .reshape(-1, nx)[:, xi - r : xi + r + 1]
                .copy()
            )

    a = readpart(off, ">f4")
    ygrid, xgrid = np.indices(a.shape)
    ok = np.hypot(xgrid + xi - r - x, ygrid + yi - r - y) < 7
    py, px = np.unravel_index(np.nanargmax(np.where(ok, a, -np.inf)), a.shape)
    bg = np.nanmedian(a)
    dt = off + 2880 * ceil(nx * ny * 4 / 2880)
    dq = None
    for i in range(5):
        with section(url, dt, dt + 17279) as src:
            b = BytesIO(src.read())
        hh = fits.Header.fromfile(b)
        start = dt + b.tell()
        n = (
            prod(hh.get(f"NAXIS{j}", 0) for j in range(1, hh["NAXIS"] + 1))
            * abs(hh["BITPIX"])
            // 8
        )
        if hh["EXTNAME"] == "DQ":
            dq = readpart(start, ">i4").astype("int64") + hh.get("BZERO", 0)
            break
        dt = start + 2880 * ceil(n / 2880)
    out = dict(
        url=url,
        utc=p["DATE-BEG"],
        filter=p["FILTER"],
        dither=p["PATT_NUM"],
        detector=p["DETECTOR"],
        predicted_x=x,
        predicted_y=y,
        peak_x=px + xi - r,
        peak_y=py + yi - r,
        peak_mjysr=a[py, px],
        peak_dq=dq[py, px] if dq is not None else np.nan,
        peak_fraction_3x3=(a[py, px] - bg)
        / np.nansum(a[py - 1 : py + 2, px - 1 : px + 2] - bg),
    )
    return out


if __name__ == "__main__":
    m = pd.read_csv("onset_dither_astrometry.csv")
    urls = [
        s.replace("_i2d.fits", suffix)
        for s in m[(m.label == 0) & m["filter"].isin(["F277W", "F356W", "F444W"])].url
        for suffix in ["_cal.fits", "_o001_crf.fits"]
    ]
    with ThreadPoolExecutor(9) as ex:
        res = list(ex.map(patch, urls))
    pd.DataFrame(res).to_csv("onset_detector.csv", index=False)
    print(
        pd.DataFrame(res)
        .iloc[::2][
            ["filter", "dither", "peak_x", "peak_y", "peak_dq", "peak_fraction_3x3"]
        ]
        .round(2)
        .to_string(index=False)
    )

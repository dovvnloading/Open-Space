"""Fetch SCI pixels for all nominated sites. Stage-3 FITS public uncompressed byte ranges."""

from pathlib import Path
from io import BytesIO
import hashlib, re, time, sys
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from urllib.request import Request, urlopen
from concurrent.futures import ThreadPoolExecutor
from catalog_screen import readcat, project
from scipy.spatial import cKDTree

BASE = "https://stpubdata.s3.amazonaws.com/"
root = Path("stamps")
center_file = Path("sites.csv")
root.mkdir(exist_ok=True)


def section(url, a, b):
    req = Request(url, headers={"Range": f"bytes={a}-{b}"})
    r = urlopen(req, timeout=150)
    if r.status != 206 or not r.headers.get("Content-Range", "").startswith(
        f"bytes {a}-{b}/"
    ):
        raise RuntimeError(("range", url, r.status, r.headers.get("Content-Range")))
    return r


def download(row):
    name, obs, band = re.search(r"jw(\d+)-o(\d+).*clear-f(\d+)w_", row.key).groups()
    dest = root / f"{int(name)}-{int(obs)}-{band}"
    if dest.with_suffix(".npz").exists():
        return
    url = BASE + row.key
    with section(url, 0, 57599) as conn:
        hdrs = BytesIO(conn.read())
    primary = fits.Header.fromfile(hdrs)
    h = fits.Header.fromfile(hdrs)
    off = hdrs.tell()
    w = WCS(h, naxis=2)
    nx = h["NAXIS1"]
    ny = h["NAXIS2"]
    pix = abs(h["CDELT1"]) * 3600
    assert h["EXTNAME"] == "SCI" and h["BITPIX"] == -32 and h["BUNIT"] == "MJy/sr"
    centers = pd.read_csv(center_file)
    xy = w.all_world2pix(centers[["ra", "dec"]].values, 0)
    radius = int(np.ceil(1.8 / pix))
    xc, yc = np.rint(xy).astype(int).T
    valid = (xc >= radius) & (yc >= radius) & (xc < nx - radius) & (yc < ny - radius)
    ids = centers.site.values[valid]
    xc = xc[valid]
    yc = yc[valid]
    xy = xy[valid]
    if not len(ids):
        return
    stamps = np.empty((len(ids), 2 * radius + 1, 2 * radius + 1), "f4")
    ranges = []
    for y in sorted(set(yc)):
        if ranges and y - radius < ranges[-1][1] + 70:
            ranges[-1][1] = max(ranges[-1][1], y + radius + 1)
        else:
            ranges.append([y - radius, y + radius + 1])
    for y0, y1 in ranges:
        # Slice response in strips to avoid loading full SW science planes into RAM.
        with section(url, off + y0 * nx * 4, off + y1 * nx * 4 - 1) as src:
            for first in range(y0, y1, 64):
                last = min(y1, first + 64)
                buf = src.read((last - first) * nx * 4)
                im = np.frombuffer(buf, ">f4").reshape(last - first, nx)
                for ix in np.flatnonzero((yc + radius >= first) & (yc - radius < last)):
                    bot = yc[ix] - radius
                    lo = max(first, bot)
                    hi = min(last, bot + 2 * radius + 1)
                    stamps[ix, lo - bot : hi - bot] = im[
                        lo - first : hi - first, xc[ix] - radius : xc[ix] + radius + 1
                    ]
    np.savez_compressed(
        dest.with_suffix(".npz"),
        site=ids,
        xy=xy,
        data=stamps,
        header=str(h.tostring(sep="\n", padding=False)),
        x0=xc - radius,
        y0=yc - radius,
        url=url,
        date=primary["DATE-BEG"],
    )


if __name__ == "__main__":
    q = (
        pd.read_csv("multiband.csv")
        .groupby("group")[["ra", "dec"]]
        .mean()
        .reset_index()
        .rename(columns={"group": "site"})
    )
    # Blank offset controls measure local difference-aperture noise on the same stamps.
    q.to_csv("sites.csv", index=False)
    products = pd.concat([pd.read_csv(f"{pid}_index.csv") for pid in [8060, 12577]])
    rows = list(products[products.key.str.endswith("i2d.fits")].itertuples())
    rows.sort(
        key=lambda r: (
            not ("o001" in r.key or "o002" in r.key),
            not ("f356" in r.key or "f150" in r.key),
            r.key,
        )
    )

    def retry(row):
        for i in range(4):
            try:
                return download(row)
            except Exception as e:
                if i == 3:
                    print(row.key, repr(e), flush=True)

    with ThreadPoolExecutor(5) as pool:
        list(pool.map(retry, rows))

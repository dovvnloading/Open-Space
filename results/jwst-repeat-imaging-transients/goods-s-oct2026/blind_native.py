"""Look for a fixed bright flagged native detector pixel in each signed two-band trigger."""

from pathlib import Path
from io import BytesIO
from math import ceil
import json, sys, warnings
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from get_stamps import section
from individual import associations, retry
from s3 import base
from concurrent.futures import ThreadPoolExecutor

warnings.simplefilter("ignore", FITSFixedWarning)
mode = "nuclear" if "--nuclear" in sys.argv else "blind"
q = (
    pd.read_csv("blind_pairs.csv")
    if mode == "blind"
    else pd.read_csv("nuclear_candidates.csv").rename(
        columns={"difference_f277w_ap_njy": "z277"}
    )
)


def run(row):
    pid, obs, band, asn, key = row
    url = base + key.replace("_i2d.fits", "_cal.fits")
    with section(url, 0, 57599) as src:
        b = BytesIO(src.read())
    p = fits.Header.fromfile(b)
    h = fits.Header.fromfile(b)
    off = b.tell()
    ww = WCS(h, naxis=2)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    x, y = ww.all_world2pix(q.ra.values, q.dec.values, 0)
    sel = (
        (x > 10)
        & (y > 10)
        & (x < nx - 10)
        & (y < ny - 10)
        & (((q.z277.values > 0) == (pid == 12577)) if mode == "blind" else True)
    )
    if not sel.any():
        return []
    x = x[sel]
    y = y[sel]
    sites = q.site.values[sel]
    lo = max(0, int(y.min()) - 10)
    hi = min(ny, int(y.max()) + 11)
    hoff = off + 2880 * ceil(nx * ny * 4 / 2880)
    dq = None
    for i in range(5):
        with section(url, hoff, hoff + 14399) as src:
            b = BytesIO(src.read())
            hh = fits.Header.fromfile(b)
        start = hoff + b.tell()
        n = (
            np.prod([hh[f"NAXIS{j}"] for j in range(1, hh["NAXIS"] + 1)])
            * abs(hh["BITPIX"])
            // 8
        )
        if hh["EXTNAME"] == "DQ":
            with section(url, start + lo * nx * 4, start + hi * nx * 4 - 1) as src:
                dq = np.frombuffer(src.read(), ">i4").reshape(hi - lo, nx).astype(
                    "int64"
                ) + hh.get("BZERO", 0)
            break
        hoff = start + 2880 * ceil(n / 2880)
    assert dq is not None
    with section(url, off + lo * nx * 4, off + hi * nx * 4 - 1) as src:
        im = np.frombuffer(src.read(), ">f4").reshape(hi - lo, nx)
    rows = []
    ys, xs = np.indices((17, 17))
    rad = 7
    for i, xx, yy in zip(sites, x, y):
        ix = int(np.rint(xx))
        iy = int(np.rint(yy))
        a = im[iy - 8 - lo : iy + 9 - lo, ix - 8 : ix + 9]
        ok = (np.hypot(xs + ix - 8 - xx, ys + iy - 8 - yy) < rad) & np.isfinite(a)
        if not ok.any():
            continue
        py, px = np.unravel_index(np.nanargmax(np.where(ok, a, -np.inf)), a.shape)
        bg = np.nanmedian(a)
        X = ix + px - 8
        Y = iy + py - 8
        f = a[py, px] - bg
        sum3 = np.nansum(a[max(0, py - 1) : py + 2, max(0, px - 1) : px + 2] - bg)
        rows.append(
            dict(
                site=i,
                pid=pid,
                filter=f"F{band}W",
                dither=p["PATT_NUM"],
                detector=p["DETECTOR"],
                x=X,
                y=Y,
                dq=int(dq[Y - lo, X]),
                pointiness=f / sum3 if sum3 > 0 else np.nan,
                url=url,
            )
        )
    return rows


jobs = []
for pid in [12577, 8060]:
    for obs in ([1] if mode == "blind" else [1, 2]):
        jobs.extend(associations((pid, obs)))
rows = []
with ThreadPoolExecutor(5) as ex:
    for r in ex.map(lambda j: retry(run, j), jobs):
        rows.extend(r)
df = pd.DataFrame(rows)
df.to_csv(mode + "_native_checks.csv", index=False)
r = (
    df[(df.dq != 0) & (df.pointiness > 0.70)]
    .groupby(["site", "detector", "x", "y"])
    .agg(
        n_filter=("filter", "nunique"),
        n_dither=("dither", "nunique"),
        n_exposures=("dq", "size"),
    )
    .reset_index()
)
r = r[(r.n_filter >= 2) & (r.n_dither >= 2)]
r.to_csv(mode + "_fixed_pixels.csv", index=False)
print("repeated flagged pixel groups", len(r), flush=True)

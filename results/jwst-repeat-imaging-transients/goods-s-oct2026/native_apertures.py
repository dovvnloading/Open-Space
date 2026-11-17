"""Check native SCI/DQ in small patches around individually measured nuclear positions."""

from io import BytesIO
from math import ceil, prod
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import warnings, time
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from get_stamps import section

warnings.simplefilter("ignore", FITSFixedWarning)
positions = pd.read_csv("sites.csv").set_index("site")
dithers = pd.read_csv("dither_fluxes.csv")
dithers = dithers[dithers.site < 1000].drop_duplicates(["site", "url"])


def read_group(entry):
    url, rows = entry
    url = url.replace("_i2d.fits", "_cal.fits")
    with section(url, 0, 57599) as r:
        b = BytesIO(r.read())
        p = fits.Header.fromfile(b)
        h = fits.Header.fromfile(b)
    off = b.tell()
    w = WCS(h, naxis=2)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    hoff = off + 2880 * ceil(nx * ny * 4 / 2880)
    for i in range(6):
        with section(url, hoff, hoff + 14399) as r:
            b = BytesIO(r.read())
            hh = fits.Header.fromfile(b)
        dqoff = hoff + b.tell()
        size = (
            prod(hh.get(f"NAXIS{j}", 0) for j in range(1, hh["NAXIS"] + 1))
            * abs(hh["BITPIX"])
            // 8
        )
        if hh["EXTNAME"] == "DQ":
            break
        hoff = dqoff + 2880 * ceil(size / 2880)
    assert hh["EXTNAME"] == "DQ"
    ret = []
    for row in rows.itertuples():
        s = positions.loc[row.site]
        x, y = w.all_world2pix(s.ra, s.dec, 0)
        xx, yy = int(np.rint(x)), int(np.rint(y))
        lo = yy - 5
        hi = yy + 6

        def read_patch(offset, t):
            with section(url, offset + lo * nx * 4, offset + hi * nx * 4 - 1) as r:
                return (
                    np.frombuffer(r.read(), t)
                    .reshape(11, nx)[:, xx - 5 : xx + 6]
                    .copy()
                )

        a = read_patch(off, ">f4")
        dq = read_patch(dqoff, ">i4").astype("int64") + hh.get("BZERO", 0)
        yg, xg = np.indices(a.shape)
        rad = np.hypot(xg + xx - 5 - x, yg + yy - 5 - y)
        core = rad <= 2.5
        ij = np.unravel_index(np.nanargmax(np.where(core, a, -np.inf)), a.shape)
        ret.append(
            dict(
                site=row.site,
                pid=row.pid,
                filter=row.filter,
                utc=row.utc,
                detector=p["DETECTOR"],
                dither=p["PATT_NUM"],
                predicted_x=x,
                predicted_y=y,
                peak_x=xx - 5 + ij[1],
                peak_y=yy - 5 + ij[0],
                peak_dq=dq[ij],
                n_core=np.count_nonzero(core),
                n_core_nonzero_dq=np.count_nonzero(dq[core]),
                max_sci_mjysr=a[ij],
                url=url,
            )
        )
    return ret


if __name__ == "__main__":
    res = []

    def retry(item):
        for i in range(3):
            try:
                return read_group(item)
            except Exception:
                if i == 2:
                    raise
                time.sleep(i + 1)

    with ThreadPoolExecutor(8) as ex:
        for group in ex.map(retry, dithers.groupby("url")):
            res.extend(group)
            pd.DataFrame(res).to_csv("nuclear_dq_apertures.csv", index=False)
    print(len(res), "native apertures", flush=True)

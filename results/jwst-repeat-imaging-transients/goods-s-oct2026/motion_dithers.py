"""Astrometry of the single-epoch emergence relative to nearby objects in the single-filter L3 frame."""

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from catalog_screen import readcat, project
from get_stamps import section
from s3 import keys, base
import numpy as np, pandas as pd, warnings, sys, re, json, time
from pathlib import Path
from io import BytesIO
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen

warnings.simplefilter("ignore", FITSFixedWarning)
siteid = int(sys.argv[1]) if len(sys.argv) > 1 else 1004
pid = int(sys.argv[2]) if len(sys.argv) > 2 else 12577
site = (
    pd.concat([pd.read_csv("onsets.csv"), pd.read_csv("sites.csv")])
    .query("site==@siteid")
    .iloc[0]
)
obs = (
    int(pd.read_csv("multiband.csv").query("group==@siteid").obs.iloc[0])
    if siteid < 1000
    else int(site.obs)
)
allrefs = {}
scale = 0.025
dy, dx = np.mgrid[-17:18, -17:18] * scale
R = np.hypot(dx, dy)
gap = (R > 0.30) & (R <= 0.4)
px = np.deg2rad(scale / 3600) ** 2 * 1e15
for band in ["115", "150", "200", "277", "356", "444"]:
    c = readcat(pid, obs, band)
    sep = (
        np.hypot(
            (c["sky_centroid.ra"] - site.ra) * np.cos(np.deg2rad(site.dec)),
            c["sky_centroid.dec"] - site.dec,
        )
        * 3600
    )
    sel = (
        (sep > 0.8)
        & (~c.is_extended)
        & (c.aper_total_flux / c.aper_total_flux_err > 28)
        & (c.nn_dist > 5)
    )
    ref = c.loc[sel].copy()
    ref["dist"] = sep[sel]
    ref = ref.sort_values("dist").head(8)
    allrefs[band] = [(0, site.ra, site.dec)] + [
        (int(r["label"]), r["sky_centroid.ra"], r["sky_centroid.dec"])
        for _, r in ref.iterrows()
    ]


def jobs():
    pref = f"jwst/public/jw{pid:05}/asn/jw{pid:05}-o{obs:03}"
    out = []
    for k, n in keys(pref):
        if "image3" not in k:
            continue
        d = json.load(urlopen(base + k, timeout=80))
        for p in d["products"]:
            band = re.search(r"clear-f(\d+)w$", p["name"]).group(1)
            for f in p["members"]:
                name = f["expname"].replace("_cal.fits", "_i2d.fits")
                out.append(
                    (band, f"jwst/public/jw{pid:05}/" + name.split("_")[0] + "/" + name)
                )
    return out


def run(item):
    band, key = item
    url = base + key
    with section(url, 0, 57599) as src:
        buf = BytesIO(src.read())
    ph = fits.Header.fromfile(buf)
    h = fits.Header.fromfile(buf)
    off = buf.tell()
    w = WCS(h, naxis=2)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    points = allrefs[band]
    ra, dec = site.ra, site.dec
    x, y = w.all_world2pix(ra, dec, 0)
    if not (15 < x < nx - 15 and 15 < y < ny - 15):
        return []
    grid = []
    for label, ra, dec in points:
        xx, yy = w.all_world2pix(
            ra + dx / 3600 / np.cos(np.deg2rad(dec)), dec + dy / 3600, 0
        )
        x0 = int(np.floor(xx.min()))
        x1 = int(np.ceil(xx.max())) + 1
        y0 = int(np.floor(yy.min()))
        y1 = int(np.ceil(yy.max())) + 1
        if min(x0, y0) < 1 or x1 >= nx or y1 >= ny:
            continue
        grid.append((label, xx, yy, x0, x1, y0, y1))
    if not grid or grid[0][0] != 0:
        return []
    lo = min(t[5] for t in grid)
    hi = max(t[6] for t in grid)
    with section(url, off + lo * nx * 4, off + hi * nx * 4 - 1) as src:
        raw = np.frombuffer(src.read(), ">f4").reshape(hi - lo, nx)
    rows = []
    for label, xx, yy, x0, x1, y0, y1 in grid:
        im = map_coordinates(raw, [yy - lo, xx], order=1, mode="constant", cval=np.nan)
        good = np.isfinite(im) & (R < 0.24)
        if (
            np.isfinite(im[gap]).mean() < 0.92
            or np.mean(np.isfinite(im[R < 0.24])) < 0.95
        ):
            continue
        bg = np.nanmedian(im[gap])
        e = max(1e-10, 1.4826 * np.nanmedian(abs(im[gap] - bg)))
        sig = (int(band) * 0.00034 + 0.015) / 2.355

        def res(par):
            x, y, am, b, sx, sy = par
            mod = (
                b
                + sx * dx[good]
                + sy * dy[good]
                + am
                * np.exp(-0.5 * ((dx[good] - x) ** 2 + (dy[good] - y) ** 2) / sig**2)
            )
            return (mod - im[good]) / e

        a = np.nanmax(im[good]) - bg
        p = least_squares(
            res,
            [0, 0, max(a, 1e-10), bg, 0, 0],
            bounds=(
                [-0.15, -0.15, 0, -np.inf, -np.inf, -np.inf],
                [0.15, 0.15, np.inf, np.inf, np.inf, np.inf],
            ),
            max_nfev=80,
        )
        mx, my, amp = p.x[:3]
        aper = R < 0.15
        apflux = np.sum(im[aper] - bg) * px
        rows.append(
            dict(
                label=label,
                filter=f"F{band}W",
                utc_start=ph["DATE-BEG"],
                utc_end=ph["DATE-END"],
                detector=ph["DETECTOR"],
                file=key.split("/")[-1],
                url=url,
                east_arcsec=mx,
                north_arcsec=my,
                peak_to_pixel_mad=amp / e,
                ap_njy=apflux,
            )
        )
    return rows


if __name__ == "__main__":
    j = jobs()
    out = []

    def checked(job):
        for i in range(3):
            try:
                return run(job)
            except Exception as e:
                if i == 2:
                    print(job, repr(e))
                    return []

    with ThreadPoolExecutor(10) as ex:
        for r in ex.map(checked, j):
            out.extend(r)
    df = pd.DataFrame(out)
    df.to_csv(
        (
            "onset_dither_astrometry.csv"
            if siteid == 1004 and pid == 12577
            else f"astrometry-{siteid}-{pid}.csv"
        ),
        index=False,
    )
    print(len(j), len(df), sum(df.label == 0), flush=True)

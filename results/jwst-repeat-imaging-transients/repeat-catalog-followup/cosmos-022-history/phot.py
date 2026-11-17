"""Coverage-verified fixed-sky aperture photometry for COSMOS candidate 022."""

from io import BytesIO
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request, urlopen
from hashlib import sha256
import time, json, sys, warnings
from threading import Lock
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates
from census import RA, DEC

warnings.simplefilter("ignore", FITSFixedWarning)
STEP = 0.025
yy, xx = np.mgrid[-61:62, -61:62] * STEP
r = np.hypot(xx, yy)
PIX = 1e15 * np.deg2rad(STEP / 3600) ** 2
sky = (r > 0.36) & (r < 0.64)
comparisons = [
    (
        np.hypot(xx - dx, yy - dy) < 0.16,
        (np.hypot(xx - dx, yy - dy) > 0.36) & (np.hypot(xx - dx, yy - dy) < 0.64),
    )
    for dx in np.arange(-0.9, 0.901, 0.3)
    for dy in np.arange(-0.9, 0.901, 0.3)
    if 0.55 < np.hypot(dx, dy) < 1.1
]


def byte_range(url, start, size):
    req = Request(url, headers={"Range": f"bytes={start}-{start+size-1}"})
    with urlopen(req, timeout=45) as res:
        if res.status != 206 or not res.headers.get("Content-Range", "").startswith(
            f"bytes {start}-{start+size-1}/"
        ):
            raise ValueError("range not honored")
        data = res.read()
    if len(data) != size:
        raise ValueError("truncated FITS range")
    return data


def url_for(oid, level):
    if level == 2:
        visit = oid.split("_")[0]
        return f"https://stpubdata.s3.amazonaws.com/jwst/public/{oid[:7]}/{visit}/{oid}_i2d.fits"
    return f'https://stpubdata.s3.amazonaws.com/jwst/public/{oid[:7]}/L3/t/{oid.split("-")[1].split("_")[0]}/{oid}_i2d.fits'


def aperture(image, area, ann):
    if np.isfinite(image[area]).mean() < 0.99 or np.isfinite(image[ann]).mean() < 0.95:
        return np.nan
    return np.nansum(image[area] - np.nanmedian(image[ann])) * PIX


def run(arg):
    obs, level = arg
    url = url_for(obs.obs_id, level)
    out = dict(
        obs_id=obs.obs_id, level=level, url=url, mjd_start=obs.t_min, mjd_end=obs.t_max
    )
    try:
        header_data = byte_range(url, 0, 57600)
        b = BytesIO(header_data)
        p = fits.Header.fromfile(b)
        h = fits.Header.fromfile(b)
        off = b.tell()
        assert h["EXTNAME"] == "SCI" and h["BITPIX"] == -32 and h["BUNIT"] == "MJy/sr"
        band = (
            p["FILTER"]
            if p["PUPIL"] in ["CLEAR", "CLEARP"] or not p["PUPIL"].startswith("F")
            else p["PUPIL"]
        )
        out.update(
            filter=band,
            optics=p["FILTER"] + "+" + p["PUPIL"],
            utc=p["DATE-BEG"],
            detector=p.get("DETECTOR"),
            dither=p.get("PATT_NUM"),
            header_sha256=sha256(header_data).hexdigest(),
        )
        w = WCS(h, naxis=2)
        x, y = w.all_world2pix(
            RA + xx / 3600 / np.cos(np.deg2rad(DEC)), DEC + yy / 3600, 0
        )
        nx, ny = h["NAXIS1"], h["NAXIS2"]
        if x.min() < 0 or x.max() >= nx - 1 or y.min() < 0 or y.max() >= ny - 1:
            return {**out, "status": "edge"}
        lo = int(np.floor(y.min()))
        hi = int(np.ceil(y.max())) + 1
        off += nx * 4 * lo
        size = nx * 4 * (hi - lo)
        raw = byte_range(url, off, size)
        image = map_coordinates(
            np.frombuffer(raw, ">f4").reshape(-1, nx), [y - lo, x], order=1, cval=np.nan
        )
        out.update(
            start_byte=off,
            bytes=size,
            pixels_sha256=sha256(raw).hexdigest(),
            finite_fraction=np.isfinite(image[r < 0.28]).mean(),
        )
        f = aperture(image, r < 0.16, sky)
        if not np.isfinite(f):
            return {**out, "status": "insufficient_finite_pixels"}
        ann = image[sky]
        weights = np.where(
            (r < 0.22) & np.isfinite(image), np.maximum(image - np.nanmedian(ann), 0), 0
        )
        if weights.sum() <= 0:
            return {**out, "status": "no_centroid"}
        ex = (xx * weights).sum() / weights.sum()
        no = (yy * weights).sum() / weights.sum()
        rr = np.hypot(xx - ex, yy - no)
        blank = np.array([aperture(image, a, bg) for a, bg in comparisons])
        sigma = (
            1.4826 * np.nanmedian(abs(blank - np.nanmedian(blank)))
            if np.isfinite(blank).sum() >= 10
            else np.nan
        )
        return {
            **out,
            "status": "measured",
            "fixed_njy": f,
            "centered_njy": aperture(image, rr < 0.16, (rr > 0.36) & (rr < 0.64)),
            "wide_njy": aperture(image, r < 0.28, sky),
            "local_scatter_njy": sigma,
            "east_arcsec": ex,
            "north_arcsec": no,
        }
    except Exception as e:
        return {**out, "status": type(e).__name__ + ": " + str(e)}


if __name__ == "__main__":
    level = int(sys.argv[1])
    obs = pd.read_csv(f"level{level}.csv")
    outfile = Path(f"phot_l{level}.csv")
    old = pd.read_csv(outfile).to_dict("records") if outfile.exists() else []
    good = {
        r["obs_id"]
        for r in old
        if not any(s in str(r["status"]) for s in ["Error", "Timeout"])
    }
    rows = [r for r in old if r["obs_id"] in good]
    inputs = [(o, level) for o in obs.itertuples() if o.obs_id not in good]
    from concurrent.futures import as_completed

    with ThreadPoolExecutor(8) as pool:
        futures = [pool.submit(run, o) for o in inputs]
        for i, f in enumerate(as_completed(futures)):
            rows.append(f.result())
            if i % 10 == 0:
                pd.DataFrame(rows).to_csv(outfile, index=False)
    pd.DataFrame(rows).sort_values("mjd_start").to_csv(outfile, index=False)
    print(
        level, len(rows), pd.DataFrame(rows).status.value_counts().to_dict(), flush=True
    )

"""Independent dither-level measurements from the uncombined Level-2 i2d products."""

from s3 import keys, base
from urllib.request import urlopen
from io import BytesIO
import numpy as np, pandas as pd, json, re, warnings, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from astropy.wcs import WCS, FITSFixedWarning
from astropy.io import fits
from scipy.ndimage import map_coordinates
from get_stamps import section

warnings.simplefilter("ignore", FITSFixedWarning)
keep = [93, 99, 185, 525, 829, 1004]
sites = (
    pd.concat(
        [pd.read_csv("sites.csv"), pd.read_csv("onsets.csv")[["site", "ra", "dec"]]]
    )
    .set_index("site")
    .loc[keep]
)
gy, gx = np.mgrid[-17:18, -17:18] * 0.035
rr = np.hypot(gx, gy)
ann = (rr >= 0.36) & (rr <= 0.51)
px = (np.deg2rad(0.035 / 3600)) ** 2 * 1e15


def associations(args):
    pid, obs = args
    prefix = f"jwst/public/jw{pid:05d}/asn/jw{pid:05d}-o{obs:03d}"
    ks = [k for k, n in keys(prefix) if "image3" in k]
    out = []
    for k in ks:
        doc = json.load(urlopen(base + k, timeout=50))
        for p in doc["products"]:
            m = re.search(r"clear-f(277|356)w$", p["name"])
            if not m:
                continue
            for f in p["members"]:
                if f["exptype"].lower() != "science":
                    continue
                name = f["expname"].replace("_cal.fits", "_i2d.fits")
                out.append(
                    (
                        pid,
                        obs,
                        m[1],
                        k,
                        "jwst/public/jw"
                        + f"{pid:05d}/"
                        + name.split("_")[0]
                        + "/"
                        + name,
                    )
                )
    return out


def cut(t):
    pid, obs, band, asn, key = t
    url = base + key
    with section(url, 0, 57599) as r:
        b = BytesIO(r.read())
    ph = fits.Header.fromfile(b)
    h = fits.Header.fromfile(b)
    off = b.tell()
    w = WCS(h, naxis=2)
    nx = h["NAXIS1"]
    ny = h["NAXIS2"]
    rows = []
    for site, row in sites.iterrows():
        xx, yy = w.all_world2pix(
            row.ra + gx / 3600 / np.cos(np.deg2rad(row.dec)), row.dec + gy / 3600, 0
        )
        x0 = int(np.floor(xx.min()))
        x1 = int(np.ceil(xx.max())) + 1
        y0 = int(np.floor(yy.min()))
        y1 = int(np.ceil(yy.max())) + 1
        if min(x0, y0) < 1 or x1 >= nx or y1 >= ny:
            continue
        with section(url, off + y0 * nx * 4, off + y1 * nx * 4 - 1) as src:
            arr = np.frombuffer(src.read(), ">f4").reshape(y1 - y0, nx)[:, x0:x1]
            im = map_coordinates(
                arr, [yy - y0, xx - x0], order=1, mode="constant", cval=np.nan
            )
        if not np.isfinite(im[(rr < 0.3) | ann]).all():
            continue
        bg = np.nanmedian(im[ann])
        out = dict(
            site=site,
            pid=pid,
            obs=obs,
            filter=f"F{band}W",
            utc=ph["DATE-BEG"],
            detector=ph["DETECTOR"],
            exposure=key.split("/")[-1],
            asn_key=asn,
            url=url,
            sha256_pixels=__import__("hashlib")
            .sha256(np.asarray(arr, dtype=">f4").tobytes())
            .hexdigest(),
        )
        for rad in [0.15, 0.22]:
            out[f"f{round(100*rad)}_njy"] = (im[rr <= rad] - bg).sum() * px
        rows.append(out)
    return rows


def retry(fn, args):
    for i in range(3):
        try:
            return fn(args)
        except Exception as e:
            if i == 2:
                print("failed", args, repr(e), flush=True)
                return []
            time.sleep(2 * i + 1)


if __name__ == "__main__":
    jobs = []
    with ThreadPoolExecutor(6) as pool:
        for a in pool.map(
            lambda a: retry(associations, a),
            [(p, o) for p in [8060, 12577] for o in [1, 2]],
        ):
            jobs.extend(a)
    print("exposures to check", len(jobs), flush=True)
    output = []
    with ThreadPoolExecutor(12) as pool:
        for a in pool.map(lambda a: retry(cut, a), jobs):
            output.extend(a)
            if output:
                pd.DataFrame(output).to_csv("dither_fluxes.csv", index=False)
    print("measurements", len(output), flush=True)

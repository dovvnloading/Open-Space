"""Independent native-detector check of the single F164N trigger in April 2025."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("repo", type=Path)
p.add_argument("cache", type=Path)
p.add_argument("--freeze", action="store_true")
args = p.parse_args()
src = args.repo / "results/jwst-nearby-companions"
sys.path.insert(0, str(src / "miri"))
sys.path.insert(0, str(src / "vetting"))
from fetch_inputs import download
from search import refine, smooth_response
from archival_images import aperture

out = Path(__file__).resolve().parent
manifest = out / "dither_inputs.json"
obs = "jw06584-o005_t001_nircam_f150w2-f164n"
with fits.open(args.cache / (obs + "_i2d.fits")) as f:
    names = [
        row["FILENAME"].replace("_o005_crf.fits", "_cal.fits")
        for row in f["HDRTAB"].data
        if row["DETECTOR"] == "NRCB4"
    ]
assert len(names) == 4
if args.freeze:
    if manifest.exists():
        raise FileExistsError(manifest)
    prefix = "jwst/public/jw06584/jw06584005001"
    url = "https://s3.amazonaws.com/stpubdata?" + urlencode(
        {"list-type": 2, "max-keys": 1000, "prefix": prefix}
    )
    body = urlopen(url, timeout=90).read()
    listing = ET.fromstring(body)
    assert listing.findtext("{*}IsTruncated") == "false"
    out.joinpath("dither_listing.xml").write_bytes(body)
    rows = []
    for item in listing.findall("{*}Contents"):
        key = item.findtext("{*}Key")
        if key.split("/")[-1] in names:
            rows.append(
                dict(
                    filename=key.split("/")[-1],
                    url="https://stpubdata.s3.amazonaws.com/" + key,
                    bytes=int(item.findtext("{*}Size")),
                    etag=item.findtext("{*}ETag"),
                )
            )
    assert len(rows) == 4
    data = dict(
        query=url,
        listing_sha256=hashlib.sha256(body).hexdigest(),
        retrieved_utc=datetime.now(timezone.utc).isoformat(),
        product_type="Stage 2 cal, NRCB4",
    )
else:
    data = json.loads(manifest.read_text())
    rows = data.pop("inputs")


def fetch(row):
    result = download(row, args.cache)
    if result["bytes"] != row["bytes"]:
        raise ValueError(row)
    with (args.cache / row["filename"]).open("rb") as f:
        parts = [hashlib.md5(b).digest() for b in iter(lambda: f.read(50 * 2**20), b"")]
    if hashlib.md5(b"".join(parts)).hexdigest() + "-" + str(len(parts)) != row[
        "etag"
    ].strip('"'):
        raise ValueError(row["filename"])
    return result


with ThreadPoolExecutor(4) as pool:
    rows = list(pool.map(fetch, rows))
if args.freeze:
    manifest.write_text(json.dumps(dict(data, inputs=rows), indent=2) + "\n")
epochs = list(csv.DictReader((out / "epochs.csv").open()))
info = next(t for t in epochs if t["obs"] == obs)
peaks = list(csv.DictReader((out / "peaks.csv").open()))
unmatched = [
    r
    for r in peaks
    if r["obs"] == obs and int(r["fixed_wins"]) == int(r["host_wins"]) == 0
]
assert len(unmatched) == 1
trigger = unmatched[0]
cosdec = np.cos(np.deg2rad(float(info["dec"])))
cat = Table.read(args.cache / (obs + "_cat.ecsv"))
ra_ref = cat["sky_centroid"].ra.deg
de_ref = cat["sky_centroid"].dec.deg
sky = np.c_[
    (ra_ref - float(info["ra"])) * 3600 * cosdec, (de_ref - float(info["dec"])) * 3600
]
dr = np.linalg.norm(sky, axis=1)
nn = cKDTree(sky).query(sky, k=2)[0][:, 1]
good = (
    (dr > 4)
    & (dr < 28)
    & (nn > 0.45)
    & (~np.array(cat["is_extended"]))
    & (cat["aper_total_flux"] > 15e-6)
)
refids = np.where(good)[0]
qra = float(trigger["ra_registered"]) - float(info["east_tie_shift"]) / 3600 / cosdec
qde = float(trigger["dec_registered"]) - float(info["north_tie_shift"]) / 3600


def local_peak(data, p, half):
    x, y = np.rint(p).astype(int)
    a = data[y - half : y + half + 1, x - half : x + half + 1]
    if a.shape != (2 * half + 1, 2 * half + 1) or not np.any(np.isfinite(a)):
        return np.array([np.nan, np.nan])
    sm = ndi.gaussian_filter(np.nan_to_num(a, nan=np.nanmedian(a)), 0.8)
    dy, dx = np.unravel_index(np.argmax(sm[1:-1, 1:-1]), sm[1:-1, 1:-1].shape)
    return np.array([x - half, y - half]) + refine(sm, dx + 1, dy + 1)


def write(name, rows):
    with (out / name).open("w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)


result = []
tieslog = []
for r in rows:
    with fits.open(args.cache / r["filename"]) as h:
        image = h["SCI"].data
        dq = h["DQ"].data
        bad = ((dq & 1) != 0) | ~np.isfinite(image)
        w = WCS(h["SCI"].header, naxis=2)
        hx, hy = map(
            float, w.world_to_pixel_values(float(info["ra"]), float(info["dec"]))
        )
        hostoff = local_peak(image, [hx, hy], 12) - [hx, hy]
        ref = np.array(w.world_to_pixel_values(ra_ref[good], de_ref[good])).T
        delta = np.array([local_peak(image, u + hostoff, 5) - u for u in ref])
        use = np.isfinite(delta).all(axis=1)
        for _ in range(3):
            off = np.median(delta[use], axis=0)
            dist = np.linalg.norm(delta - off, axis=1)
            use &= dist < max(0.7, 3 * np.median(dist[use]))
        off = np.median(delta[use], axis=0)
        rms = np.sqrt(np.mean(np.sum((delta[use] - off) ** 2, axis=1)))
        if use.sum() < 8 or rms > 1.2:
            raise ValueError((r["filename"], use.sum(), rms))
        for k in np.flatnonzero(use):
            tieslog.append(
                dict(
                    file=r["filename"],
                    catalog_label=int(cat["label"][refids[k]]),
                    dx=delta[k, 0],
                    dy=delta[k, 1],
                )
            )
        hostpt = np.array([hx, hy]) + off
        target = np.array(w.world_to_pixel_values(qra, qde)) + off
        # Mask DO_NOT_USE pixels before photometry. Inspect all native DQ bits separately.
        data = np.where(bad, np.nan, image)
        scale = (abs(np.linalg.det(w.pixel_scale_matrix)) ** 0.5) * 3600
        cal = float(info["host_flux_uJy"]) / aperture(data, *hostpt, 3)
        rawflux = cal * aperture(data, *target, 3)
        controls = []
        radius = np.linalg.norm(target - hostpt)
        for angle in np.linspace(0, 2 * np.pi, 32, endpoint=False):
            pos = hostpt + radius * np.array([np.cos(angle), np.sin(angle)])
            if np.linalg.norm(pos - target) > 20:
                controls.append(cal * aperture(data, *pos, 3))
        med = np.nanmedian(controls)
        flux = rawflux - med
        noise = 1.4826 * np.nanmedian(abs(np.array(controls) - med))
        err = np.hypot(noise, 0.05 * flux)
        x, y = np.rint(target).astype(int)
        yy, xx = np.mgrid[y - 3 : y + 4, x - 3 : x + 4]
        mask = np.hypot(xx - target[0], yy - target[1]) <= 3
        bits = int(np.bitwise_or.reduce(dq[yy[mask], xx[mask]]))
        cut = data[y - 18 : y + 19, x - 18 : x + 19]
        ff = smooth_response(np.nan_to_num(cut, nan=np.nanmedian(cut)))
        local = np.zeros(ff.shape, bool)
        local[1:-1, 1:-1] = True
        local &= ff == ndi.maximum_filter(ff, 5)
        iy, ix = np.where(local)
        xd = ix + x - 18
        yd = iy + y - 18
        distances = np.hypot(xd - target[0], yd - target[1])
        k = np.argmin(distances)
        result.append(
            dict(
                file=r["filename"],
                target_x=target[0],
                target_y=target[1],
                n_ties=int(use.sum()),
                rms_tie_pix=rms,
                pixscale_arcsec=scale,
                dx_pix=off[0],
                dy_pix=off[1],
                forced_flux_uJy=flux,
                azimuth_error_uJy=err,
                flux_over_error=flux / err,
                nearest_response_local_max_arcsec=distances[k] * scale,
                aperture_DQ_bitwise_or=bits,
                aperture_bad_pixels=int(bad[yy[mask], xx[mask]].sum()),
                aperture_pixels=int(mask.sum()),
            )
        )
write("dither_photometry.csv", result)
write("dither_registration.csv", tieslog)

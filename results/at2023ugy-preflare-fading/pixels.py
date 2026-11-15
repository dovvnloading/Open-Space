"""Verify the pre-event fading in NEOWISE L1b detector images."""
import argparse
import csv
import hashlib
import itertools
import math
import urllib.request
import warnings
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares

from measure import TARGET, ipac, write_csv

ARCHIVE = "https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/p1bm_frm"
OFFSETS = np.arange(-60, 61, dtype=float)
X, Y = np.meshgrid(OFFSETS, OFFSETS)
RA, DEC = TARGET
COS_DEC = math.cos(math.radians(DEC))
SKY_RA, SKY_DEC = RA + X / (3600 * COS_DEC), DEC + Y / 3600

def read_image(frame, cache):
    scan, number, band = frame["scan"], frame["frame"], frame["band"]
    filename = f"{scan}{number}-w{band}-int-1b.fits"
    path = cache / filename
    if not path.exists():
        url = f"{ARCHIVE}/{scan[-2:]}/{scan}/{number}/{filename}"
        with urllib.request.urlopen(url, timeout=120) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != frame["sha256"]:
            raise ValueError(f"Unexpected image digest for {filename}")
        path.write_bytes(data)
    if hashlib.sha256(path.read_bytes()).hexdigest() != frame["sha256"]:
        raise ValueError(f"Unexpected cached image digest for {filename}")
    with fits.open(path) as hdus:
        hdr = hdus[0].header
        image = hdus[0].data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            x, y = WCS(hdr).all_world2pix(SKY_RA, SKY_DEC, 0)
        inside = (x >= 1) & (x < image.shape[1]-2) & (y >= 1) & (y < image.shape[0]-2)
        if inside.mean() < .98:
            return None
        patch = map_coordinates(image, [y, x], order=1, mode="constant", cval=np.nan)
        return (patch - np.nanmedian(patch)) * 10 ** ((20 - hdr["MAGZP"]) / 2.5)

def psf_amplitude(image, sx, sy, free_center=False):
    mask = (X-sx)**2 + (Y-sy)**2 <= 8**2
    x, y, values = X[mask]-sx, Y[mask]-sy, image[mask]
    valid = np.isfinite(values)
    x, y, values = x[valid], y[valid], values[valid]
    if len(values) < 50:
        raise ValueError("Incomplete aperture")
    def model(dx, dy):
        profile = np.exp(-((x-dx)**2 + (y-dy)**2) / (2 * 2.85**2))
        return np.column_stack((profile, np.ones_like(x), x, y))
    def resid(center):
        basis = model(*center)
        return basis @ np.linalg.lstsq(basis, values, rcond=None)[0] - values
    center = least_squares(resid, [0, 0], bounds=(-4, 4)).x if free_center else (0., 0.)
    amplitude = np.linalg.lstsq(model(*center), values, rcond=None)[0][0]
    return float(amplitude), float(center[0]), float(center[1])

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True, type=Path, help="Cache already populated by measure.py")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    image_cache = args.cache / "images"
    image_cache.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    with Path(__file__).with_name("frames.csv").open() as handle:
        frames = list(csv.DictReader(handle))
    mode = defaultdict(list)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for row, data in zip(frames, pool.map(lambda r: read_image(r, image_cache), frames)):
            if data is not None:
                mode[int(row["year"]), row["season"], row["band"]].append(data)
    stars = list(ipac(args.cache / "allwise.tbl"))
    pos = [((float(s["ra"])-RA)*3600*COS_DEC, (float(s["dec"])-DEC)*3600) for s in stars]
    neighbor = pos[min((i for i, s in enumerate(stars) if float(s["dist"]) > 2.5),
                       key=lambda i: float(stars[i]["dist"]))]
    blank_xy = [(x,y) for x,y in itertools.product(range(-46,47,12), repeat=2)
                if all(math.hypot(x-u,y-v) > 12 for u,v in pos)]
    output = []
    for band, season in itertools.product(("1","2"), ("spring","fall")):
        stacks = [mode[y,season,band] for y in (2014,2015,2021,2022)]
        if min(map(len, stacks)) < 2:
            raise ValueError("Insufficient image coverage in an epoch")
        epochs = [np.nanmedian(v,axis=0) for v in stacks]
        diff = np.nanmedian(epochs[:2],axis=0) - np.nanmedian(epochs[2:],axis=0)
        blanks = [psf_amplitude(diff,x,y)[0] for x,y in blank_xy]
        target = psf_amplitude(diff,0,0)[0]
        _, dx, dy = psf_amplitude(diff,0,0,free_center=True)
        output.append(dict(band="W"+band, season=season, selected=16, used=sum(map(len,stacks)),
                           target_amplitude=target, delta_ra_cosdec_arcsec=dx, delta_dec_arcsec=dy,
                           centroid_offset_arcsec=math.hypot(dx,dy),
                           neighbor_amplitude=psf_amplitude(diff,*neighbor)[0],
                           blank_n=len(blanks), blank_median=float(np.median(blanks)),
                           blank_sd=float(np.std(blanks,ddof=1)), blank_max_abs=max(map(abs,blanks))))
    write_csv(args.out / "pixels.csv", output)

if __name__ == "__main__":
    main()

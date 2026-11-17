"""Blind annual F356W subtraction on two fixed 1536-pixel detector interiors."""

from pathlib import Path
from io import BytesIO
from math import ceil
import hashlib, json, sys, tempfile, subprocess, warnings
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates
import numpy as np, pandas as pd
from get_stamps import section
from s3 import base

warnings.simplefilter("ignore", FITSFixedWarning)


def read_map(args):
    key, target = args
    url = base + key
    with section(url, 0, 57599) as r:
        source = BytesIO(r.read())
    primary = fits.Header.fromfile(source)
    h = fits.Header.fromfile(source)
    off = source.tell()
    w = WCS(h, naxis=2)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    x, y = w.all_world2pix(*target, 0)
    low = max(0, int(np.floor(y.min())) - 1)
    high = min(ny, int(np.ceil(y.max())) + 2)
    planes = []
    hashes = []
    for ext in ["SCI", "ERR"]:
        if ext == "ERR":
            hoff = off + 2880 * ceil(nx * ny * 4 / 2880)
            for i in range(6):
                with section(url, hoff, hoff + 14399) as r:
                    hdr = BytesIO(r.read())
                hh = fits.Header.fromfile(hdr)
                off = hoff + hdr.tell()
                if hh["EXTNAME"] == "ERR":
                    break
                total = (
                    np.prod([hh[f"NAXIS{i}"] for i in range(1, hh["NAXIS"] + 1)])
                    * abs(hh["BITPIX"])
                    // 8
                )
                hoff = off + int(2880 * ceil(total / 2880))
            assert hh["EXTNAME"] == "ERR" and hh["BITPIX"] == -32
        with section(url, off + low * nx * 4, off + high * nx * 4 - 1) as r:
            raw = r.read()
        a = np.frombuffer(raw, ">f4").reshape(-1, nx)
        planes.append(map_coordinates(a, [y - low, x], order=1, cval=np.nan))
        hashes.append(
            dict(
                url=url,
                extname=ext,
                offset=off + low * nx * 4,
                bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest(),
            )
        )
        del a, raw
    return planes, hashes


def main(obs=1, band="356", side="west"):
    # Products are selected by coverage and one fixed detector-coordinate rectangle per side, independent of triggers.
    key = (
        lambda p: f"jwst/public/jw{p:05d}/L3/t/o{obs:03d}/jw{p:05d}-o{obs:03d}_t{obs:03d}_nircam_clear-f{band}w_i2d.fits"
    )
    with section(base + key(12577), 0, 57599) as r:
        b = BytesIO(r.read())
    p = fits.Header.fromfile(b)
    h = fits.Header.fromfile(b)
    w = WCS(h, naxis=2)
    ny = nx = 1536
    yy, xx = np.indices((ny, nx))
    x0 = 272 if side == "west" else h["NAXIS1"] - 272 - nx
    y0 = (h["NAXIS2"] - ny) // 2
    target = w.all_pix2world(xx + x0, yy + y0, 0)
    w.wcs.crpix -= np.array([x0, y0])
    hh = w.to_header()
    hh["PIXAR_SR"] = h["PIXAR_SR"]
    hh["BUNIT"] = "MJy/sr"
    out = Path(f"blind-o{obs}-f{band}-{side}")
    out.mkdir(exist_ok=True)
    paths = []
    snap = []
    with tempfile.TemporaryDirectory(dir=".") as tmp:
        for pid in [8060, 12577]:
            planes, hashes = read_map((key(pid), target))
            snap.extend(hashes)
            for ext, a in zip(["sci", "err"], planes):
                f = Path(tmp) / f"{pid}-{ext}.fits"
                fits.writeto(f, np.asarray(a, "f4"), hh)
                paths.append(str(f))
            del planes, a
        pd.DataFrame(snap).to_csv(out / "image_ranges.csv", index=False)
        script = next(
            p
            for p in [
                Path("../screen.py"),
                Path("../repo/results/jwst-repeat-imaging-transients/screen.py"),
            ]
            if p.is_file()
        )
        subprocess.run([sys.executable, str(script), str(out), *paths], check=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])

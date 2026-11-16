"""Measure local r-band depth and quality across the X-ray error circle from Legacy Survey maps."""

import argparse
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def analyze(root, brick, ra, dec, radius, flux, cache, band):
    cached = Path(cache)
    cached.mkdir(parents=True, exist_ok=True)

    def fetch(kind):
        suffix = kind if kind == "maskbits" else f"{kind}-{band}"
        name = f"legacysurvey-{brick}-{suffix}.fits.fz"
        url = f"https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/{root}/coadd/{brick[:3]}/{brick}/{name}"
        path = cached / name
        if not path.exists():
            with urlopen(url, timeout=90) as res:
                content = res.read()
            if not content.startswith(b"SIMPLE"):
                raise RuntimeError(url)
            path.write_bytes(content)
        with fits.open(path) as hdus:
            hdu = hdus[1]
            w = WCS(hdu.header).celestial
            x, y = w.world_to_pixel_values(ra, dec)
            if not (15 < x < hdu.shape[1] - 16 and 15 < y < hdu.shape[0] - 16):
                raise ValueError("Position outside usable brick " + str((x, y)))
            # Native brick 0.262 arcsec/pix. Use the sky transform for the precise selection.
            scale = min(
                abs(w.proj_plane_pixel_scales()[i].value) * 3600 for i in (0, 1)
            )
            d = int(math.ceil(radius / scale)) + 2
            i, j = round(float(x)), round(float(y))
            ys, xs = np.mgrid[j - d : j + d + 1, i - d : i + d + 1]
            ar, ad = w.pixel_to_world_values(xs, ys)
            sep = (
                np.hypot(
                    ((ar - ra + 180) % 360 - 180) * math.cos(math.radians(dec)),
                    ad - dec,
                )
                * 3600
            )
            vals = hdu.section[j - d : j + d + 1, i - d : i + d + 1][sep <= radius]
            if len(vals) == 0 or not np.all(np.isfinite(vals)):
                raise RuntimeError("Invalid image")
            return kind, {
                "url": url,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "samples": len(vals),
                "min": float(vals.min()),
                "max": float(vals.max()),
                "median": float(np.median(vals)),
                "bitwise_or": int(np.bitwise_or.reduce(vals.astype(np.int64)))
                if kind == "maskbits"
                else None,
            }

    with ThreadPoolExecutor(4) as p:
        maps = dict(p.map(fetch, ["depth", "maskbits", "nexp", "psfsize"]))
    ivar = maps["depth"]["min"]
    result = dict(ra=ra, dec=dec, radius_arcsec=radius, band=band, maps=maps)
    if ivar > 0:
        nmgy = 5 / math.sqrt(ivar)
        nu = 2.99792458e18 / ({"r": 6231.0, "g": 4770.0, "z": 9134.0}[band])
        result.update(
            five_sigma_nmgy=nmgy,
            ab_limit=22.5 - 2.5 * math.log10(nmgy),
            flux_x=flux,
            fx_over_nufnu=flux / (nu * nmgy * 3.631e-29),
        )
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("root", "brick"):
        p.add_argument(k)
    for k in ("ra", "dec", "radius", "flux"):
        p.add_argument(k, type=float)
    p.add_argument("cache")
    p.add_argument("output")
    p.add_argument("--band", choices=["r", "g", "z"], default="r")
    a = p.parse_args()
    result = analyze(a.root, a.brick, a.ra, a.dec, a.radius, a.flux, a.cache, a.band)
    Path(a.output).write_text(json.dumps(result, indent=2) + "\n")
    print(
        {k: v for k, v in result.items() if k != "maps"},
        {k: (v["min"], v["max"]) for k, v in result["maps"].items()},
    )

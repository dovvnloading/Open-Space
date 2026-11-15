"""Sample two published, fixed pulsar positions in the original RACS HiPS tiles."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

from replay import position

DIGESTS = {
    "low": "ad4906e7b229bd4765dd545c814df44c2aa6aa3b7011d2c63946ae1849267844",
    "mid": "371f593d214456440e1280940508ba1792297d804b553fda934986d3018e9991",
}


def sample(image, wcs, ra, dec):
    x, y = wcs.world_to_pixel_values(ra, dec)
    return map_coordinates(image, [np.atleast_1d(y), np.atleast_1d(x)],
                           order=1, mode="constant", cval=np.nan)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tile_dir", type=Path, help="Directory containing low.fits and mid.fits")
    args = parser.parse_args()
    ra, dec = position("MeerKAT 2024")
    alt_ra, alt_dec = position("MeerKAT thesis 2023")
    results, background_differences, signal_differences = {}, [], []
    for band, digest in DIGESTS.items():
        tile = args.tile_dir / f"{band}.fits"
        if hashlib.sha256(tile.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected {tile} digest")
        with fits.open(tile) as hdus:
            img, wcs = hdus[0].data, WCS(hdus[0].header)
            x, y = wcs.world_to_pixel_values(ra, dec)
            yy, xx = np.indices(img.shape)
            radius = np.hypot(xx-x, yy-y)
            mask = (radius > 25) & (radius < 90)
            # Both tiles have the same WCS. Use the same sky controls in each.
            cra, cdec = wcs.pixel_to_world_values(xx[mask][::4], yy[mask][::4])
            bg = sample(img, wcs, cra, cdec)
            bg_other = sample(img, wcs, cra + (alt_ra-ra), cdec + (alt_dec-dec))
            baseline = np.median(bg)
            scatter = 1.4826 * np.median(abs(bg-baseline))
            src = sample(img, wcs, [ra, alt_ra], [dec, alt_dec])
            a, b = (src-baseline)/scatter
            signal_differences.append(float(a-b))
            background_differences.append((bg-bg_other)/scatter)
            results[band] = {"normalized_at_2024": round(float(a), 3),
                             "normalized_at_2023": round(float(b), 3)}
    controls = np.sum(background_differences, axis=0)
    results["translated_positions"] = len(controls)
    results["contrast_at_least_as_large"] = int(np.sum(controls >= sum(signal_differences)))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

"""Inject a ten-noise-scale point source at eight field positions in eligible images."""

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS


def trial(row, cache):
    with fits.open(cache / "fits" / f"{row.id}.fits", memmap=False) as hdus:
        header = hdus[0].header
        image = hdus[0].data.squeeze().astype(float) * 1000
        wcs = WCS(header).celestial
    x0, y0 = wcs.world_to_pixel_values(row.ra_obs, row.dec_obs)
    yy, xx = np.indices(image.shape)
    radius = np.hypot(xx - x0, yy - y0)
    bg = np.nanmedian(image[(radius >= 8) & (radius <= 35)])
    pa = np.deg2rad(header["BPA"])
    result = []
    for i in range(8):
        x = x0 + 20 * np.cos(np.pi * i / 4)
        y = y0 + 20 * np.sin(np.pi * i / 4)
        choose = np.isfinite(image) & (np.hypot(xx - x, yy - y) <= 1.5)
        if not np.any(choose):
            continue
        east, north = (wcs.pixel_scale_matrix * 3600) @ np.stack(
            [xx[choose] - x, yy[choose] - y]
        )
        a = (east * np.sin(pa) + north * np.cos(pa)) / (header["BMAJ"] * 3600)
        b = (east * np.cos(pa) - north * np.sin(pa)) / (header["BMIN"] * 3600)
        signal = 10 * np.exp(-4 * np.log(2) * (a * a + b * b))
        result.append(np.max((image[choose] - bg) / row.rms_mjy + signal) >= 6)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    d = pd.read_csv(args.cache / "screen-final.csv")
    d = d[(d.sp >= 20) & (d.status == "measured")]
    totals = [v for r in d.itertuples() for v in trial(r, args.cache)]
    print(
        json.dumps(
            dict(
                injected_amplitude_in_local_noise_units=10,
                detection_threshold=6,
                TY_catalog_positions=len(d),
                injections=len(totals),
                recovered=sum(map(int, totals)),
                experimental_fraction=np.mean(totals),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

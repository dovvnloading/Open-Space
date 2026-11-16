"""Test whether the published 36.528 MHz map can distinguish ILTJ112022.08+502936.2.

Provide the NASA LAMBDA Mollweide FITS and an ASTRON LoTSS-DR3 0.6 deg cone.
The map statistic describes beam-scale sky contrast and is not a flux measurement
of the target. The neighbour diagnostic uses only measured 144 MHz fluxes.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

CENTER = SkyCoord(170.092166 * u.deg, 50.493238 * u.deg)


def contrast(centers, image, wcs):
    # Background at 35-50 arcmin; eight position angles per radius.
    controls = []
    for radius in (35, 42.5, 50):
        for pa in range(0, 360, 45):
            sample = centers.directional_offset_by(pa * u.deg, radius * u.arcmin)
            x, y = wcs.world_to_pixel(sample.galactic)
            controls.append(map_coordinates(image, [y, x], order=1, prefilter=False))
    x, y = wcs.world_to_pixel(centers.galactic)
    return map_coordinates(image, [y, x], order=1, prefilter=False) - np.median(
        controls, axis=0
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("image", type=Path)
    ap.add_argument("neighbours", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    with fits.open(args.image) as h:
        coords = CENTER.directional_offset_by(
            np.repeat(np.arange(0, 360, 15), 2) * u.deg, np.tile([1.4, 2.2], 24) * u.deg
        )
        coords = SkyCoord([CENTER] + list(coords))
        values = contrast(coords, h[1].data, WCS(h[1].header))
    data = pd.read_csv(args.neighbours)
    distance = CENTER.separation(
        SkyCoord(data.ra.to_numpy() * u.deg, data.dec.to_numpy() * u.deg)
    ).arcmin
    target = distance < 0.25
    # Geometric-mean circular FWHM of the published +45-degree beam at 36.528 MHz.
    width = np.sqrt(20.2 * 16.9)
    response = np.exp(-4 * np.log(2) * (distance / width) ** 2)
    # Table 1 gives thermal noise in K and Jy/beam at the same frequency.
    factor = 0.799 / 595
    with args.image.open("rb") as f:
        checksum = hashlib.file_digest(f, "sha256").hexdigest()
    result = {
        "frequency_mhz": 36.528,
        "source": "ILTJ112022.08+502936.2",
        "image_sha256": checksum,
        "map_brightness_unit": "K",
        "control_positions": len(values) - 1,
        "sky_contrast_K": values[0],
        "sky_contrast_control_mad_K": 1.4826
        * np.median(abs(values[1:] - np.median(values[1:]))),
        "approx_jy_per_K": factor,
        "published_thermal_noise_jy": 0.799,
        "published_confusion_estimate_jy": (36.528 / 74) ** (-0.7),
        "beam_fwhm_proxy_arcmin": width,
        "lotss_cone_rows": len(data),
        "target_144mhz_jy": data.total_flux[target].sum() / 1000,
        "weighted_neighbours_144mhz_jy": np.sum(
            data.total_flux[~target] * response[~target]
        )
        / 1000,
    }
    args.output.write_text(json.dumps(result, indent=2, default=float) + "\n")

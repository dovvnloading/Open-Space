"""Interpolate the HEASARC HI4PI Galactic atomic hydrogen column map at one position."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates


def sample(path, ra, dec):
    sky = SkyCoord(ra, dec, unit="deg").galactic
    with fits.open(path) as f:
        h = f[0]
        x, y = WCS(h.header).all_world2pix([[sky.l.deg, sky.b.deg]], 0)[0]
        c = float(map_coordinates(h.data, [[y], [x]], order=1)[0])
        return {
            "url": "https://heasarc.gsfc.nasa.gov/Images/w3nh/h1_nh_HI4PI.fits",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "ra": ra,
            "dec": dec,
            "l": sky.l.deg,
            "b": sky.b.deg,
            "N_HI_cm-2": c,
            "beam_FWHM_deg": h.header["BMAJ"],
            "method": "bilinear interpolation of projected HI4PI pixels",
        }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("map", type=Path)
    p.add_argument("ra", type=float)
    p.add_argument("dec", type=float)
    a = p.parse_args()
    print(json.dumps(sample(a.map, a.ra, a.dec), indent=2))

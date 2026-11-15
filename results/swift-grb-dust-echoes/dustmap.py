"""Sample the independent 10 pc EXPLORE differential-extinction cube."""
import argparse
import hashlib
from pathlib import Path
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from scipy.ndimage import map_coordinates

URL = ("https://zenodo.org/records/10405177/files/"
       "explore_cube_density_values_010pc_v2.fits?download=1")
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("cube", type=Path, help="FITS cube from " + URL)
parser.add_argument("output", type=Path)
args = parser.parse_args()
if hashlib.md5(args.cube.read_bytes()).hexdigest() != "b86cd29c41a9b8378994d2af923ad51d":
    raise ValueError("Mismatch with the publisher's MD5")
with fits.open(args.cube) as hd:
    h = hd[0].header
    direction = SkyCoord(253.5649583333*u.deg, -54.0445*u.deg).galactic.cartesian.xyz.value
    d = np.arange(0, 401, 5)
    center = (np.array(hd[0].data.shape[::-1]) - 1)/2   # centered, odd-sized cube
    voxels = center[:, None] + direction[:, None]*d/h["STEP"]
    values = map_coordinates(hd[0].data, voxels[::-1], order=1)
np.savetxt(args.output, np.stack([d,values],axis=1), header="distance_pc,A0_mag_per_pc",
           fmt=["%d","%.7f"], delimiter=",", comments="")

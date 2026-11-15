"""Optional Voyager-era model sensitivity using Stellarium's GUST86 shared library."""
import argparse
import ctypes
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation

from replay import PREDICTION, both, jd0, load_vectors

AU_KM = 149597870.7
VECTOR = ctypes.c_double * 3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path, help="Shared library exporting GetGust86Coor")
    args = parser.parse_args()
    states = load_vectors("prediction", PREDICTION)
    get_position = ctypes.CDLL(str(args.library.resolve())).GetGust86Coor
    get_position.argtypes = [ctypes.c_double, ctypes.c_int, VECTOR, VECTOR]
    get_position.restype = None
    # VSOP87 dynamical ecliptic to equatorial J2000, as specified in gust86.h.
    rotation = Rotation.from_euler("xz", [23.4392803055556, -0.0000275], degrees=True)

    grid = np.arange(0, 3 * 86400 + 1, 300.0)
    for body, index in [(701, 1), (702, 2)]:
        positions, velocities = [], []
        for t in grid:
            xyz, velocity = VECTOR(), VECTOR()
            get_position(jd0(PREDICTION) + t / 86400, index, xyz, velocity)
            positions.append(rotation.apply(np.array(xyz)) * AU_KM)
            velocities.append(rotation.apply(np.array(velocity)) * AU_KM / 86400)
        orbit = CubicHermiteSpline(grid, positions, velocities)
        states[body] = lambda t, orbit=orbit: states[799](t) + orbit(t)
    print(json.dumps(both(states), indent=2))


if __name__ == "__main__":
    main()

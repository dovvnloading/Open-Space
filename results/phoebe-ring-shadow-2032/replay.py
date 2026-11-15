"""Model the direct solar illumination of Phoebe using archived Horizons states and Voyager ring transmission."""
import argparse
import csv
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
JD0 = 2463292.5  # 2032 March 1, 00:00 TDB
C = 299792.458
POLE_RA, POLE_DEC = np.deg2rad([40.589 - .036 * .3217, 83.537 - .004 * .3217])
POLE = np.array([np.cos(POLE_DEC) * np.cos(POLE_RA), np.cos(POLE_DEC) * np.sin(POLE_RA), np.sin(POLE_DEC)])
REQ, RPOL, RSUN = 60268., 54364., 695700.
EYE = np.eye(3)
QUADRIC = EYE / REQ**2 + (1 / RPOL**2 - 1 / REQ**2) * np.outer(POLE, POLE)


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def angle(v, w):
    return np.arctan2(np.linalg.norm(np.cross(v, w)), np.dot(v, w))


def splines():
    curves = {}
    for path in sorted(INPUT.glob("*_from_*.csv")):
        with path.open() as f:
            rows = np.array([[float(row[k]) for k in (0, 2, 3, 4)] for row in csv.reader(f)])
        curves[path.stem] = CubicSpline(rows[:, 0] - JD0, rows[:, 1:])
    return curves


def sun_disk(n, limb=0.6):
    # Deterministic equal-area sampling on the apparent disk. Linear optical limb darkening.
    i = np.arange(n) + .5
    r = np.sqrt(i / n)
    theta = i * (np.pi * (3 - np.sqrt(5)))
    return r * np.cos(theta), r * np.sin(theta), 1 - limb * (1 - np.sqrt(1 - r**2))


def lighting(saturn, sun, disk, profile, ring_scale=1.):
    solar_direction = unit(sun)
    east = unit(np.cross(solar_direction, [0, 0, 1]))
    north = np.cross(solar_direction, east)
    x, y, weight = disk
    tangent_radius = RSUN / np.sqrt(sun @ sun - RSUN**2)
    rays = unit(solar_direction + tangent_radius * (x[:, None] * east + y[:, None] * north))

    # Each ray is tested at the nearest point to Saturn in the oblate body's ellipsoid metric.
    coefficient = np.einsum("ij,ij->i", rays @ QUADRIC, rays)
    distances = (rays @ (QUADRIC @ saturn)) / coefficient
    offsets = distances[:, None] * rays - saturn
    globe = (distances > 0) & (np.einsum("ij,ij->i", offsets @ QUADRIC, offsets) <= 1)
    cosines = rays @ POLE
    ring_distances = (saturn @ POLE) / cosines
    radius = np.linalg.norm(ring_distances[:, None] * rays - saturn, axis=1)
    radial_bins, normal_tau = profile
    tau = np.interp(radius, radial_bins, normal_tau, left=0., right=0.)
    transmission = np.where(ring_distances > 0, np.exp(-tau * ring_scale / np.abs(cosines)), 1.)
    return np.average((~globe) * transmission, weights=weight), np.average(~globe, weights=weight)


def utc_from_days(t):
    # Assume no additional leap seconds beyond the known TAI-UTC = 37 s.
    return (datetime(2032, 3, 1) + timedelta(days=t, seconds=-69.184)).isoformat(timespec="seconds")


def calculate(n=2048, limb=.6, ring_scale=1., step=5, alternate=False, cassini=False):
    curves = splines()
    disk = sun_disk(n, limb)
    tab = np.loadtxt(INPUT / "PS1P01.TAB", delimiter=",")
    # 99 flags optical depths above the measurement range. Either 99 or the published lower
    # limit already implies an essentially opaque ray at the geometry of this event.
    profile = (tab[:, 0], np.maximum(tab[:, 3], 0))
    if cassini:
        tab = np.loadtxt(INPUT / "cassini-profile.csv", delimiter=",", skiprows=1)
        profile = (tab[:, 0], np.maximum(tab[:, 1], 0))
    correction = None
    if alternate:
        d = np.loadtxt(INPUT / "ph20_minus_sat441.csv", delimiter=",", skiprows=1)
        correction = CubicSpline(d[:, 0] - JD0, d[:, 1:])
    out = []
    for i in range(0, 72 * 60 - 1, step):
        received = i / 1440 + 69.184 / 86400
        moon = curves["phoebe_from_earth"](received)
        emitted = received - np.linalg.norm(moon) / (C * 86400)
        if emitted < 0:
            continue
        saturn, sun = [curves[key](emitted) for key in ("saturn_from_phoebe", "sun_from_phoebe")]
        if correction:
            delta = correction(emitted)
            saturn, sun = saturn - delta, sun - delta
        fraction, globe = lighting(saturn, sun, disk, profile, ring_scale)
        separation = angle(moon, curves["saturn_from_earth"](received)) * 206264.806
        elongation = np.rad2deg(angle(moon, curves["sun_from_earth"](received)))
        out.append((utc_from_days(received), fraction, globe, separation, elongation))
    return out


def print_summary(rows):
    f, globe = np.array([r[1:3] for r in rows]).T
    for label, mask in [("Ring model first and last below 95%", f < .95), ("Globe-only first and last below 99%", globe < .99), ("Globe-only first and last below 0.1%", globe < .001)]:
        inds = np.flatnonzero(mask)
        print(label, rows[inds[0]][0], rows[inds[-1]][0])
    for j in find_peaks(f, prominence=.08)[0]:
        print("Local recovery", rows[j][0], round(f[j], 3), round(globe[j], 3))
    for j in find_peaks(-f, prominence=.08)[0]:
        print("Local minimum", rows[j][0], round(f[j], 3), round(globe[j], 3))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rays", type=int, default=2048)
    parser.add_argument("--limb", type=float, default=.6)
    parser.add_argument("--ring-scale", type=float, default=1.)
    parser.add_argument("--step", type=int, default=5, help="minutes between samples")
    parser.add_argument("--cassini", action="store_true", help="use the independent VIMS profile")
    parser.add_argument("--alternate", action="store_true", help="substitute the independent PH20 orbit")
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    rows = calculate(args.rays, args.limb, args.ring_scale, args.step, args.alternate, args.cassini)
    if args.csv:
        with args.csv.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["earth_utc", "direct_fraction", "globe_only_fraction", "saturn_separation_arcsec", "elongation_deg"])
            writer.writerows((r[0], *[f"{v:.6f}" for v in r[1:]]) for r in rows)
    print_summary(rows)

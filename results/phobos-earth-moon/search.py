"""Screen the two sampled rover routes. Recheck candidates with replay.py."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import spiceypy as spice
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent


def unit(v):
    return v / np.linalg.norm(v, axis=-1)[..., None]


def interpolator(body, spacing, start, end):
    epochs = np.arange(start - 2 * spacing, end + 3 * spacing, spacing)
    positions, _ = spice.spkpos(body, epochs, "IAU_MARS", "LT+S", "MARS")
    return CubicSpline(epochs, positions, axis=0)


def line_fit(direction, velocity, delta):
    projected_velocity = velocity - np.dot(velocity, direction) * direction
    distance = delta - (delta @ direction)[:, None] * direction
    time = -(distance @ projected_velocity) / np.dot(projected_velocity, projected_velocity)
    miss = np.linalg.norm(distance + time[:, None] * projected_velocity, axis=1)
    return time, miss, np.linalg.norm(projected_velocity)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    args = parser.parse_args()
    for name in ["naif0012.tls", "pck00011.tpc", "de440s.bsp", "mar099s.bsp"]:
        spice.furnsh(str(args.cache / name))
    start, end = map(spice.str2et, ["2026 NOV 15", "2037 JAN 1"])
    phobos = interpolator("PHOBOS", 600, start, end)
    earth, moon, sun = [interpolator(b, 3600, start, end) for b in ["EARTH", "MOON", "SUN"]]
    epochs = np.arange(start, end, 120)
    with (HERE / "search-sites.csv").open() as f:
        rows = list(csv.DictReader(f))
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["crater", "linearized_utc", "sol", "joint_impact_km", "earth_el_deg", "sun_el_deg"])
    for crater, radius in [("Gale", 3392.187), ("Jezero", 3392.53)]:
        points = [p for p in rows if p["crater"] == crater]
        sites = np.array([spice.latrec(radius, np.radians(float(p["lon"])), np.radians(float(p["lat"])))
                          for p in points])
        site = sites[-1]
        up = unit(site)
        for indices in np.array_split(np.arange(len(epochs)), 100):
            extended = np.arange(max(0, indices[0] - 1), min(len(epochs), indices[-1] + 2))
            ts = epochs[extended]
            pv = phobos(ts) - site
            directions = unit(earth(ts) - site)
            misses = np.linalg.norm(np.cross(pv, directions), axis=1)
            for idx in find_peaks(-misses)[0]:
                if not indices[0] <= extended[idx] <= indices[-1] or misses[idx] >= 300:
                    continue
                if pv[idx] @ directions[idx] <= 0 or directions[idx] @ up <= 0.07:
                    continue
                epoch = ts[idx]
                if unit(sun(epoch) - site) @ up >= -0.035:
                    continue
                pos = phobos(epoch)
                velocity = phobos(epoch + 1) - pos
                te, be, speed = line_fit(directions[idx], velocity, pos - sites)
                tm, bm, _ = line_fit(unit(moon(epoch) - site), velocity, pos - sites)
                offsets = te[:, None] + (tm - te)[:, None] * np.linspace(0, 1, 30)
                joint = np.minimum.reduce(np.maximum(
                    np.hypot(be[:, None], (offsets - te[:, None]) * speed),
                    np.hypot(bm[:, None], (offsets - tm[:, None]) * speed),
                ), axis=1)
                best = np.argmin(joint)
                if joint[best] >= 13:
                    continue
                time = epoch + te[best]
                el = np.degrees(np.arcsin(unit(earth(time) - site) @ up))
                solar = np.degrees(np.arcsin(unit(sun(time) - site) @ up))
                writer.writerow([crater, spice.et2utc(time, "ISOC", 2), points[best]["sol"],
                                 round(joint[best], 2), round(el, 2), round(solar, 2)])


if __name__ == "__main__":
    main()

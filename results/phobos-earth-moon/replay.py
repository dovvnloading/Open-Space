"""Reproduce finite-disc occultation contacts over Gale crater from archived kernels."""
import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import spiceypy as spice
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
BODIES = {"EARTH": 6378.137, "MOON": 1737.4}  # km; use Earth's larger, equatorial radius
CASES = [
    ("Gale sol 1399", 137.3562752, -4.68457822, "2035-02-05T06:25:00"),
    ("Gale sol 5073", 137.388886, -4.82602331, "2035-02-05T06:25:00"),
    ("Gale sol 5073", 137.388886, -4.82602331, "2029-10-27T01:32:38"),
    ("Gale sol 5073", 137.388886, -4.82602331, "2034-05-26T21:51:26"),
    ("Jezero sol 1907", 77.22012875, 18.41996953, "2026-07-02T10:14:00"),
]


def unit(vector):
    return vector / np.linalg.norm(vector)


def radius_at(grid, directions):
    directions = np.atleast_2d(directions)
    lat = np.degrees(np.arcsin(directions[:, 2]))
    lon = np.degrees(np.arctan2(directions[:, 1], directions[:, 0]))
    rows = -lat * 128 - 0.5
    columns = (lon - 90) * 128 - 0.5
    if (rows < 0).any() or (columns < 0).any() or (rows >= 5632).any() or (columns >= 11520).any():
        raise ValueError("Terrain samples fall outside MOLA tile megr00n090hb")
    samples = map_coordinates(grid, [rows, columns], output="float64", order=1)
    return 3396 + samples / 1000


def basis(up):
    east = unit(np.cross([0, 0, 1], up))
    return east, np.cross(up, east)


def horizon(grid, site, azimuth_deg):
    up = unit(site)
    east, north = basis(up)
    bearing = math.radians(azimuth_deg)
    along = math.sin(bearing) * east + math.cos(bearing) * north
    distance = np.arange(0.125, 100.01, 0.125)
    ground_radius = radius_at(grid, up)[0]
    arc = distance / ground_radius
    points = np.cos(arc)[:, None] * up + np.sin(arc)[:, None] * along
    rays = points * radius_at(grid, points)[:, None] - site
    elevations = np.degrees(np.arctan2(rays @ up, np.hypot(rays @ east, rays @ north)))
    index = np.argmax(elevations)
    return elevations[index], distance[index]


def vector(body, epoch, site, frame="J2000", correction="CN"):
    state, light_time = spice.spkcpo(body, epoch, frame, "OBSERVER", correction,
                                    site, "MARS", "IAU_MARS")
    return state[:3], light_time


def sky(body, epoch, site):
    v, _ = vector(body, epoch, site, "IAU_MARS", "CN+S")
    up = unit(site)
    east, north = basis(up)
    azimuth = np.degrees(np.arctan2(v @ east, v @ north)) % 360
    return azimuth, np.degrees(np.arcsin(unit(v) @ up))


def fully_covered(body, epoch, site, shape, segment, angles):
    p, light_time = vector("PHOBOS", epoch, site)
    b, _ = vector(body, epoch, site)
    if p @ b <= 0:
        return False
    # Ray directions in inertial space share the observer. Omit stellar
    # aberration consistently: it cannot change whether two rays coincide.
    rotation = spice.pxform("J2000", "IAU_PHOBOS", epoch - light_time)
    origin = -(rotation @ p)
    axis = unit(rotation @ b)
    u = unit(np.cross(axis, [0, 0, 1]))
    v = np.cross(axis, u)
    radius = BODIES[body]
    tangent = radius / np.sqrt(b @ b - radius**2)
    with spice.no_found_check():
        if not spice.dskx02(shape, segment, origin, axis)[2]:
            return False
        for cosine, sine in angles:
            ray = axis + tangent * (cosine * u + sine * v)
            if not spice.dskx02(shape, segment, origin, ray)[2]:
                return False
    return True


def contacts(body, epoch, site, shape, segment, angles):
    def covered(offset):
        return fully_covered(body, epoch + offset, site, shape, segment, angles)
    times = np.arange(-30.0, 30.25, 0.25)
    flags = [covered(t) for t in times]
    edges = []
    for t, first, last in zip(times, flags, flags[1:]):
        if first == last:
            continue
        lo, hi = t, t + 0.25
        while hi - lo > 0.002:
            mid = (lo + hi) / 2
            if covered(mid) == first:
                lo = mid
            else:
                hi = mid
        edges.append((lo + hi) / 2)
    if len(edges) not in (0, 2) or flags[0] or flags[-1]:
        raise ValueError(f"Contact window needs inspection: {edges}")
    return edges


def verify_file(path, details):
    with path.open("rb") as f:
        if "ranges" in details:
            if path.stat().st_size != details["bytes"]:
                raise ValueError(f"{path.name} has an unexpected file size")
            for span in details["ranges"]:
                f.seek(span["offset"])
                data = f.read(span["length"])
                if hashlib.sha256(data).hexdigest() != span["sha256"]:
                    raise ValueError(f"{path.name} record mismatch at {span['offset']}")
        elif hashlib.file_digest(f, "sha256").hexdigest() != details["sha256"]:
            raise ValueError(f"{path.name} checksum mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path, help="directory containing the inputs in sources.json")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--limb-rays", type=int, default=64)
    parser.add_argument("--shape", choices=("Willner", "Gaskell"), default="Willner")
    parser.add_argument("--spatial", action="store_true", help="also evaluate eight displaced Gale sites")
    parser.add_argument("--ephemeris", choices=("JPL", "IMCCE", "both"), default="both")
    args = parser.parse_args()
    if args.limb_rays < 8:
        parser.error("--limb-rays must be at least 8")
    sources = json.loads((HERE / "sources.json").read_text())
    args.cache.mkdir(exist_ok=True, parents=True)
    shape_name = "phobos512.bds" if args.shape == "Gaskell" else "phobos_2014_09_22.bds"
    for name, details in sources["kernels"].items():
        if name.endswith(".bds") and name != shape_name:
            continue
        if name == "NOE-4-2020.bsp" and args.ephemeris == "JPL":
            continue
        path = args.cache / name
        if args.download and not path.exists():
            with urlopen(details["url"], timeout=180) as incoming, path.with_suffix(".part").open("wb") as out:
                shutil.copyfileobj(incoming, out)
            path.with_suffix(".part").replace(path)
        verify_file(path, details)
    for name in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "mar099s.bsp"]:
        spice.furnsh(str(args.cache / name))
    shape = spice.dasopr(str(args.cache / shape_name))
    segment = spice.dlabfs(shape)
    grid = np.memmap(args.cache / "megr00n090hb.img", dtype=">i2", mode="r", shape=(5632, 11520))
    theta = np.arange(args.limb_rays) * (2 * np.pi / args.limb_rays)
    angles = np.column_stack([np.cos(theta), np.sin(theta)])
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["orbit", "site", "reference_utc", "east_km", "north_km", "site_radius_km",
                     "earth_in_s", "earth_out_s", "moon_in_s", "moon_out_s", "joint_s",
                     "earth_az_deg", "earth_el_deg", "sun_el_deg", "horizon_deg", "ridge_km"])
    for orbit in ["JPL", "IMCCE"]:
        if args.ephemeris not in ("both", orbit):
            continue
        if orbit == "IMCCE":
            spice.furnsh(str(args.cache / "NOE-4-2020.bsp"))
        for label, lon, lat, iso in CASES:
            reference = spice.latrec(1.0, math.radians(lon), math.radians(lat))
            east, north = basis(reference)
            terrain = label.startswith("Gale")
            offsets = [(0, 0)]
            if args.spatial and label == "Gale sol 1399":
                offsets += [(x, y) for r in (2, 5) for x, y in [(r, 0), (-r, 0), (0, r), (0, -r)]]
            for dx, dy in offsets:
                up = unit(reference + (dx * east + dy * north) / 3392)
                # Add a 2 m camera height. The Jezero waypoint lacks a direct
                # radius; 3392.53 km is estimated from the neighbouring areoid.
                radius = radius_at(grid, up)[0] if terrain else 3392.53
                site = (radius + 0.002) * up
                epoch = spice.str2et(iso)
                earth, moon = [contacts(b, epoch, site, shape, segment, angles) for b in BODIES]
                joint = max(0, min(earth[1], moon[1]) - max(earth[0], moon[0])) if earth and moon else 0
                az, el = sky("EARTH", epoch, site)
                _, solar_el = sky("SUN", epoch, site)
                ridge, distance = horizon(grid, site, az) if terrain else (None, None)
                values = [dx, dy, np.linalg.norm(site), *(earth or [None]*2), *(moon or [None]*2),
                          joint, az, el, solar_el, ridge, distance]
                writer.writerow([orbit, label, iso] + [round(v, 4) if v is not None else "" for v in values])


if __name__ == "__main__":
    main()

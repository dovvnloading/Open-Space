"""Reconstruct the spherical Mars tangent geometry with the pinned NAIF kernels."""

import argparse
import csv
from datetime import datetime, timedelta
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
import spiceypy as spice

ROOT = Path(__file__).resolve().parent
SYSTEM = "SOLAR SYSTEM BARYCENTER"
REFERENCE_RADIUS_KM = 3396.19
CLOCK_ORIGIN = datetime(1970, 1, 1)
PASSES = {
    "4432": ("2024-03-17T04:29:22", "2024-03-17T04:51:28"),
    "4433": ("2024-03-17T11:13:48.954", "2024-03-17T11:48:51.854"),
    "4434": ("2024-03-18T01:27:08", "2024-03-18T01:49:14"),
    "4436": ("2024-03-19T12:08:45.954", "2024-03-19T12:43:45.854"),
    "4443": ("2024-03-24T10:56:46.954", "2024-03-24T11:31:47.854"),
    "4445": ("2024-03-26T18:51:20.954", "2024-03-26T19:26:20.854"),
    "4451": ("2024-03-30T13:43:07.954", "2024-03-30T14:17:51.854"),
    "4452": ("2024-10-27T05:16:50.947", "2024-10-27T05:50:47.847"),
}


def barycentric(body, epoch):
    return spice.spkpos(body, epoch, "J2000", "NONE", SYSTEM)[0]


def ray_at_reception(receive):
    earth = barycentric("EARTH", receive)
    relative, light_time = spice.spkpos("MARS EXPRESS", receive, "J2000", "CN", "EARTH")
    spacecraft = earth + relative
    direction = spacecraft - earth
    squared_range = np.dot(direction, direction)
    fraction = 1.0
    for _ in range(4):
        epoch = receive - fraction * light_time
        mars = barycentric("MARS", epoch)
        fraction = np.dot(mars - earth, direction) / squared_range
    vector = earth + fraction * direction - mars
    radius = np.linalg.norm(vector)
    fixed = spice.pxform("J2000", "IAU_MARS", epoch) @ vector
    _, longitude, latitude = spice.reclat(fixed)
    sun = spice.spkpos("SUN", epoch, "J2000", "NONE", "MARS")[0]
    cosine = np.dot(vector, sun) / (radius * np.linalg.norm(sun))
    zenith = math.degrees(math.acos(np.clip(cosine, -1, 1)))
    return (
        round(radius, 5), round(radius - REFERENCE_RADIUS_KM, 5),
        round(math.degrees(latitude), 5), round(math.degrees(longitude) % 360, 5),
        round(zenith, 5), round(fraction * light_time, 7), round(fraction, 9),
    )


def run(kernels, output):
    receipt = json.loads((ROOT / "kernel-receipt.json").read_bytes())
    for item in receipt:
        path = kernels / item["name"]
        if path.stat().st_size != item["bytes"] or sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Kernel does not match the receipt: {path}")
        spice.furnsh(str(path))
    try:
        with output.open("w", newline="", encoding="ascii") as stream:
            writer = csv.writer(stream)
            writer.writerow([
                "volume", "utc_reception", "ray_tangent_radius_km", "spherical_altitude_km",
                "planetocentric_lat_deg", "east_lon_deg", "solar_zenith_deg",
                "earth_to_tangent_light_s", "ray_segment_fraction",
            ])
            for volume, (start, stop) in PASSES.items():
                lower = math.ceil((datetime.fromisoformat(start) - CLOCK_ORIGIN).total_seconds() / 2) * 2
                upper = int((datetime.fromisoformat(stop) - CLOCK_ORIGIN).total_seconds())
                for tick in range(lower, upper + 1, 2):
                    receive = spice.str2et((CLOCK_ORIGIN + timedelta(seconds=tick)).isoformat())
                    writer.writerow([volume, spice.et2utc(receive, "ISOC", 3), *ray_at_reception(receive)])
    finally:
        spice.kclear()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.kernels, args.output)

"""Reproduce the CTX footprint audit from its compact observational extract."""
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

RADIUS_KM = 3389.5
ROOT = Path(__file__).resolve().parent


def vector(lon, lat):
    lon, lat = map(math.radians, (lon, lat))
    return (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def angle(a, b):
    return math.atan2(math.hypot(*cross(a, b)), dot(a, b))


def distance_to_arc(point, a, b):
    normal = cross(a, b)
    length = math.hypot(*normal)
    ends = min(angle(point, a), angle(point, b))
    if length < 1e-14:
        return ends
    normal = [v/length for v in normal]
    component = dot(point, normal)
    foot = [x-component*y for x, y in zip(point, normal)]
    scale = math.hypot(*foot)
    if scale < 1e-14:
        return ends
    foot = [v/scale for v in foot]
    if abs(angle(a, foot) + angle(foot, b) - angle(a, b)) < 1e-10:
        return min(ends, angle(point, foot))
    return ends


def miss_distance(lon, lat, rings):
    # In each footprint's local longitude domain, ray crossings test inclusion.
    # Antimeridian-split coverage records are disjoint exterior rings.
    point = vector(lon, lat)
    nearest = math.inf
    for ring in rings:
        anchor = ring[0][0]
        x = anchor + (lon-anchor+180) % 360-180
        polygon = [(anchor + (u-anchor+180) % 360-180, v) for u, v in ring]
        inside = False
        for (u, v), (w, z) in zip(polygon, polygon[1:]):
            if (v > lat) != (z > lat) and x < u + (lat-v)*(w-u)/(z-v):
                inside = not inside
        if inside:
            return 0.0
        for a, b in zip(ring, ring[1:]):
            nearest = min(nearest, distance_to_arc(point, vector(*a), vector(*b)))
    return nearest * RADIUS_KM


def main():
    with (ROOT / "endpoints.csv").open() as f:
        rows = list(csv.DictReader(f))
    covering = defaultdict(set)
    outside = []
    for row in rows:
        distance = miss_distance(float(row["lon"]), float(row["lat"]),
                                 json.loads(row["footprint_rings"]))
        if distance > 0:
            outside.append((distance, row))
        else:
            covering[row["pds_product"]].add(row["site"])

    fields = ["site", "stage", "catalog_date", "pds_product", "miss_km", "also_covers_catalog_site"]
    writer = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for distance, row in sorted(outside, key=lambda it: -it[0]):
        writer.writerow({**{k: row[k] for k in fields[:4]}, "miss_km": f"{distance:.3f}",
                         "also_covers_catalog_site": "|".join(sorted(covering[row["pds_product"]]))})


if __name__ == "__main__":
    main()

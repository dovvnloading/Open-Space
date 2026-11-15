"""Intersect fixed pulsar coordinates with the moving, topocentric Venus disk."""
import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
START = datetime(2027, 2, 6, 16, 40, tzinfo=timezone.utc)
# Right ascension in hours, declination in degrees, errors in seconds of RA time / arcsec.
POSITIONS = {
    "MeerKAT 2024": ((18, 13, 39.9057931), (21, 13, 10.22721), 0.01081, 2.81633),
    "MeerKAT thesis 2023": ((18, 13, 39.87), (21, 13, 0.4), 0.03, 0.3),
    "Parkes 2002": ((18, 13, 39.87), (21, 13, 0.0), 0.03, 7.0),
}


def degrees(parts):
    return parts[0] + parts[1] / 60 + parts[2] / 3600


def position(name):
    ra, dec, _, _ = POSITIONS[name]
    return 15 * degrees(ra), -degrees(dec)


def timestamp(seconds):
    return (START + timedelta(seconds=round(seconds))).strftime("%H:%M:%S")


def intersect(track, ra, dec):
    # Gnomonic offsets in arcseconds centered on the pulsar, east/north.
    sin_d, cos_d = math.sin(math.radians(dec)), math.cos(math.radians(dec))
    points = []
    for epoch in track:
        dra = math.radians(float(epoch["ra_deg"]) - ra)
        d = math.radians(float(epoch["dec_deg"]))
        common = math.cos(d) * math.cos(dra)
        dot = sin_d * math.sin(d) + cos_d * common
        scale = (180 / math.pi) * 3600 / dot
        east = scale * math.cos(d) * math.sin(dra)
        north = scale * (cos_d * math.sin(d) - sin_d * common)
        points.append((float(epoch["seconds"]), east, north, float(epoch["diameter_arcsec"]) / 2))
    best = (float("inf"), None)
    contacts = []
    for p, q in zip(points, points[1:]):
        t, x, y, radius = p
        dt, dx, dy, dr = [b - a for a, b in zip(p, q)]
        f = min(1.0, max(0.0, -(x * dx + y * dy) / (dx * dx + dy * dy)))
        best = min(best, (math.hypot(x + f * dx, y + f * dy), t + f * dt))
        a = dx * dx + dy * dy - dr * dr
        b = 2 * (x * dx + y * dy - radius * dr)
        c = x * x + y * y - radius * radius
        disc = b * b - 4 * a * c
        if disc >= 0:
            for root in ((-b - math.sqrt(disc)) / (2*a), (-b + math.sqrt(disc)) / (2*a)):
                if 0 <= root < 1:
                    contacts.append(t + root * dt)
    if len(contacts) not in (0, 2):
        raise ValueError("Ephemeris interval does not bracket a full chord")
    return best[0], best[1], contacts


def main():
    sites = defaultdict(list)
    with (ROOT / "horizons.csv").open() as source:
        for row in csv.DictReader(source):
            sites[row["site"]].append(row)
    results = []
    for site, track in sites.items():
        for name, (_, _, ra_error, dec_error) in POSITIONS.items():
            ra, dec = position(name)
            sep, mid, times = intersect(track, ra, dec)
            # Boundary stress test across each coordinate's +/-1 quoted error,
            # including the center. No covariance or distribution assumed.
            sensitivity = []
            for east in (-1, 0, 1):
                for north in (-1, 0, 1):
                    *_, ends = intersect(track, ra + east * ra_error / 240,
                                          dec + north * dec_error / 3600)
                    sensitivity.append(ends[-1] - ends[0] if ends else 0)
            results.append(dict(site=site, solution=name, closest_utc=timestamp(mid),
                                separation_arcsec=round(sep, 3),
                                ingress_utc=timestamp(times[0]) if times else None,
                                egress_utc=timestamp(times[1]) if times else None,
                                duration_seconds=round(times[1]-times[0], 1) if times else 0,
                                stress_min_seconds=round(min(sensitivity), 1),
                                stress_max_seconds=round(max(sensitivity), 1)))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

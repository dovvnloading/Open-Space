"""Pairwise lower limits on bound, two-body mass from relative astrometry."""

import csv
import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AU = 149597870700.0
C = 299792458.0
G = 6.67430e-11
DAY = 86400.0
ARCSEC = np.deg2rad(1 / 3600)


def table(path):
    with path.open() as file:
        return list(csv.DictReader(file))


@dataclass
class Epoch:
    name: str
    mjd_utc: float
    ra: float
    dec: float
    distance_au: float
    east: float
    north: float

    def geometry(self):
        ra, dec = np.deg2rad([self.ra, self.dec])
        east = np.array([-np.sin(ra), np.cos(ra), 0])
        north = np.array(
            [-np.cos(ra) * np.sin(dec), -np.sin(ra) * np.sin(dec), np.cos(dec)]
        )
        scale = self.distance_au * AU * ARCSEC
        return (
            (self.east * east + self.north * north) * scale,
            np.cross(east, north),
            scale,
        )


def epochs():
    headers = table(ROOT / "hst" / "measurements.csv")
    geometry = (
        (HERE / "horizons-hst.txt").read_text().split("$$SOE")[1].split("$$EOE")[0]
    )
    row = next(csv.reader(geometry.strip().splitlines()))
    hst = Epoch(
        "HST June 13",
        np.mean([float(r["midpoint_mjd_utc"]) for r in headers]),
        float(row[3]),
        float(row[4]),
        float(row[5]),
        np.mean([float(r["core_east_arcsec"]) for r in headers]),
        np.mean([float(r["core_north_arcsec"]) for r in headers]),
    )
    sources = {r["measurement_id"]: r for r in table(ROOT / "data" / "target.csv")}
    result = [hst]
    for fit in table(ROOT / "output" / "fits.csv"):
        if (
            fit["id"] not in sources
            or fit["kind"] != "Science"
            or fit["variant"] != "nominal"
        ):
            continue
        source = sources[fit["id"]]
        result.append(
            Epoch(
                source["measurement_id"],
                float(source["mjd"]) - 37 / DAY,
                float(source["ephra"]),
                float(source["ephdec"]),
                float(source["toporange"]),
                float(fit["east_arcsec"]),
                float(fit["north_arcsec"]),
            )
        )
    return sorted(result, key=lambda x: x.mjd_utc)


def lower_limit(a, b, radius_a, radius_b):
    p, los1, scale1 = a.geometry()
    q, los2, scale2 = b.geometry()
    common_axis = np.cross(los1, los2)
    common_axis /= np.linalg.norm(common_axis)
    # A projected direction orthogonal to both sightlines is independent of
    # either unmeasured component along the line of sight.
    radii = max(0.0, np.linalg.norm(p) - radius_a * scale1) + max(
        0.0, np.linalg.norm(q) - radius_b * scale2
    )
    min_chord = min(abs(np.dot(p + sign * q, common_axis)) for sign in [-1, 1])
    chord = max(0.0, min_chord - radius_a * scale1 - radius_b * scale2)
    # Use the largest time interval allowed by the +/- 0.5% range envelope.
    # Distances also scale each projected vector, with independent errors.
    range_fraction = 0.005
    chord = max(
        0.0,
        chord
        - range_fraction
        * (
            np.linalg.norm(p)
            + np.linalg.norm(q)
            + radius_a * scale1
            + radius_b * scale2
        ),
    )
    radii *= 1 - range_fraction
    light_seconds = (b.distance_au - a.distance_au) * AU / C
    dt = (b.mjd_utc - a.mjd_utc) * DAY - light_seconds
    dt += range_fraction * (a.distance_au + b.distance_au) * AU / C + 2
    # Both bounds lower the parabolic flight time monotonically.
    mass = ((pow(radii + chord, 1.5) - pow(radii - chord, 1.5)) / (6 * dt)) ** 2 / G
    return mass, radii / 1000, chord / 1000, dt / DAY


if __name__ == "__main__":
    records = []
    for a, b in itertools.combinations(epochs(), 2):
        if b.mjd_utc - a.mjd_utc < 1:
            continue
        for radius_hst, radius_rubin in [
            (0.0, 0.0),
            (0.01, 0.03),
            (0.02, 0.06),
            (0.06, 0.06),
            (0.10, 0.10),
        ]:
            mass, radii, chord, days = lower_limit(
                a,
                b,
                radius_hst if a.name.startswith("HST") else radius_rubin,
                radius_rubin,
            )
            records.append(
                dict(
                    first=a.name,
                    second=b.name,
                    hst_radius_arcsec=radius_hst,
                    rubin_radius_arcsec=radius_rubin,
                    max_elapsed_days=f"{days:.6f}",
                    min_sum_radii_km=f"{radii:.2f}",
                    min_chord_km=f"{chord:.2f}",
                    mass_floor_kg=f"{mass:.8e}",
                )
            )
    with (HERE / "bounds.csv").open("w") as out:
        writer = csv.DictWriter(out, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)

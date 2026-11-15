"""Optional narrow-date reproduction of the Venus discovery screen."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree
from skyfield.api import load, load_file

DB_HASH = "bac32471fef3e60cfd734f1db85503a9bd0a565018cb98edff54475f87eb642f"
SPK_HASH = "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"


def sexagesimal(value):
    angle = sum(abs(float(v)) / 60**i for i, v in enumerate(value.split(":")))
    return -angle if value.startswith("-") else angle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalogue", type=Path)
    parser.add_argument("kernel", type=Path)
    args = parser.parse_args()
    for path, digest in [(args.catalogue, DB_HASH), (args.kernel, SPK_HASH)]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected input bytes for {path}")
    sources, vectors = [], []
    for record in args.catalogue.read_text().split("@---"):
        fields = {}
        for line in record.splitlines():
            parts = line.split()
            if len(parts) > 1 and parts[0] in {"PSRJ", "RAJ", "DECJ"}:
                fields.setdefault(parts[0], parts[1])
        if len(fields) != 3:
            continue
        ra = np.radians(15 * sexagesimal(fields["RAJ"]))
        dec = np.radians(sexagesimal(fields["DECJ"]))
        sources.append(fields["PSRJ"])
        vectors.append((np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)))
    stars = np.array(vectors)
    tree = cKDTree(stars)
    eph = load_file(str(args.kernel))
    ts = load.timescale()
    start = ts.utc(2026, 11, 15).tt
    end = ts.utc(2029, 1, 1).tt
    grid = np.arange(start, end + 0.1, 0.125)
    grid = np.unique(np.r_[grid[grid < end], end])
    earth, venus = eph["earth"], eph["venus"]

    def track(tt):
        vec = earth.at(ts.tt_jd(tt)).observe(venus).position.km.T
        dist = np.linalg.norm(vec, axis=-1)
        return vec / np.expand_dims(dist, -1), dist

    directions, _ = track(grid)
    segments = {}
    for t0, t1, a, b in zip(grid, grid[1:], directions, directions[1:]):
        # A generous extra 2 arcminutes covers disk size, parallax and curvature.
        radius = np.linalg.norm(b-a)/2 + np.radians(2/60)
        for j in tree.query_ball_point((a+b)/2, radius):
            segments.setdefault(j, []).append((t0, t1))
    results = []
    for j, intervals in segments.items():
        def impact(tt):
            direction, distance = track(tt)
            return np.linalg.norm(direction-stars[j]) * distance
        for t0, t1 in intervals:
            fit = minimize_scalar(lambda offset: impact(t0+offset), bounds=(0, t1-t0),
                                  method="bounded", options={"xatol": 1e-9})
            if fit.fun <= 6051.8 + 6378.2 and 1e-6 < fit.x < (t1-t0)-1e-6:
                results.append({"pulsar": sources[j], "closest_utc": ts.tt_jd(t0+fit.x).utc_iso(),
                                "geocentric_impact_km": round(fit.fun)})
    print(json.dumps({"coordinate_pairs": len(sources), "start_utc": ts.tt_jd(start).utc_iso(),
                      "end_utc": ts.tt_jd(end).utc_iso(), "earth_intercepts": results}, indent=2))


if __name__ == "__main__":
    main()

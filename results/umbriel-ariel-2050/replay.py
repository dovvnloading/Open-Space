"""Compute geocentric light-reception times for mutual occultations and eclipses."""
import argparse
import gzip
import hashlib
import io
import json
from contextlib import ExitStack
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from jplephem.spk import SPK
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parent
C = 299792.458
R_SUN = 695700.0
RADIUS = {701: 578.9, 702: 584.7, 704: 761.4}
PREDICTION = datetime(2050, 3, 11)
HISTORICAL = datetime(2007, 5, 4, 13)


def jd0(epoch):
    return 2451545.0 + (epoch - datetime(2000, 1, 1, 12)).total_seconds() / 86400


def load_vectors(label, epoch):
    raw = gzip.decompress((ROOT / "inputs" / f"{label}.csv.gz").read_bytes())
    expected = json.loads((ROOT / "horizons.json").read_text())[label + "_csv_sha256"]
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"Input checksum failed for {label}")
    data = np.loadtxt(io.BytesIO(raw), delimiter=",", skiprows=1)
    vectors = {}
    for number in np.unique(data[:, 0]):
        a = data[data[:, 0] == number]
        vectors[int(number)] = CubicHermiteSpline(
            (a[:, 1] - jd0(epoch)) * 86400, a[:, 2:5], a[:, 5:8], extrapolate=False)
    return vectors


def ray(states, t, background, foreground, source):
    """t is the epoch the light leaves the background moon toward Earth."""
    q = states[background](t)
    sign = 1 if source == 399 else -1
    h = t
    for _ in range(4):
        delta = states[source](h) - q
        distance = np.linalg.norm(delta)
        unit = delta / distance
        h = t + sign * distance / C

    crossing = t
    for _ in range(4):
        delta = states[foreground](crossing) - q
        depth = unit @ delta
        crossing = t + sign * depth / C
    scale = distance / (distance - depth)
    separation = np.linalg.norm(np.cross(delta, unit)) * scale
    image_radius = RADIUS[foreground] * scale
    penumbra = 0.0 if source == 399 else R_SUN * depth / (distance - depth)
    return separation, image_radius, penumbra, h, depth


def events(states, epoch, background, foreground, lo, hi, source):
    def geometry(t):
        return ray(states, t, background, foreground, source)

    def edge(t):
        sep, image_radius, penumbra, _, depth = geometry(t)
        if depth <= 0:
            return 1.e12
        return sep - (image_radius + RADIUS[background] + penumbra)

    def receive(t):
        tdb = ray(states, t, background, foreground, 399)[3]
        return (epoch + timedelta(seconds=float(tdb))).isoformat(timespec="milliseconds")

    grid = np.arange(lo, hi + 1, 60.0)
    values = np.array([edge(t) for t in grid])
    changes = np.flatnonzero((values[:-1] < 0) != (values[1:] < 0))
    if len(changes) % 2 or values[0] < 0 or values[-1] < 0:
        raise ValueError("Window truncates an event")
    roots = [brentq(edge, grid[k], grid[k + 1]) for k in changes]
    records = []
    for start, end in zip(roots[::2], roots[1::2]):
        loc = minimize_scalar(lambda t: geometry(t)[0], bounds=(start, end),
                              method="bounded", options={"xatol": 0.001}).x
        b, image_radius, penumbra, _, depth = geometry(loc)
        radius = RADIUS[background]
        foreground_fraction = None
        if source == 399:
            # All retained occultations are partial; guard the two complete cases.
            if b >= image_radius + radius:
                foreground_fraction = 0.0
            elif b <= abs(image_radius - radius):
                foreground_fraction = min(1.0, (image_radius / radius) ** 2)
            else:
                a = np.arccos((b*b + radius**2 - image_radius**2) / (2*b*radius))
                c = np.arccos((b*b + image_radius**2 - radius**2) / (2*b*image_radius))
                area = radius**2 * (a - np.sin(2*a)/2) + image_radius**2 * (c - np.sin(2*c)/2)
                foreground_fraction = round(area / (np.pi * radius**2), 4)
        records.append(dict(kind="occultation" if source == 399 else "eclipse",
                            start_received_tdb=receive(start), minimum_received_tdb=receive(loc),
                            end_received_tdb=receive(end), duration_min=round((end-start)/60, 2),
                            impact_km=round(b, 1), occulted_background_disk_fraction=foreground_fraction,
                            foreground_depth_km=round(depth),
                            projected_core_radius_km=round(image_radius - penumbra, 1)))
    return records


def both(states):
    return sorted([event for source in (399, 10) for event in
                   events(states, PREDICTION, 701, 702, 24*3600, 51*3600, source)],
                  key=lambda row: row["start_received_tdb"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    modern = load_vectors("prediction", PREDICTION)
    control = load_vectors("historical", HISTORICAL)
    modern_events = both(modern)
    with ExitStack() as stack:
        older = dict(modern)
        for b in (701, 702):
            file = SPK.open(ROOT / "inputs" / f"ura111xl-{b}.bsp")
            stack.callback(file.close)
            seg = file[7, b]
            older[b] = lambda t, seg=seg: modern[7](t) + seg.compute(jd0(PREDICTION), t / 86400)
        older_events = both(older)

    observed = events(control, HISTORICAL, 702, 704, 3*3600, 5*3600, 399)
    # TDB - UTC at 2007-05-04 is TT-UTC (65.184 s) plus ~0.0014 s.
    observed_mid = datetime.fromisoformat(observed[0]["minimum_received_tdb"]) - timedelta(seconds=65.1854)
    observed[0]["predicted_minimum_utc"] = observed_mid.isoformat(timespec="milliseconds")
    observed[0]["measured_minimum_utc"] = "2007-05-04T19:09:52"
    observed[0]["predicted_minus_measured_s"] = round(
        (observed_mid - datetime(2007, 5, 4, 19, 9, 52)).total_seconds(), 1)
    result = {"ura184": modern_events, "ura111": older_events, "historical_control": observed}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()

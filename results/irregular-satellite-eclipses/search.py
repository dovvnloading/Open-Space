"""Search planetary eclipses of nine distant irregular moons with JPL Horizons."""

import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
C = 299792.458
DAY = 86400.0
SUN_R = 695700.0
# J2000 ICRF. IAU pole approximations adequate for pole-angle sensitivity tests.
PLANETS = {
    599: (71492.0, 66854.0, 268.057, 64.495),
    699: (60268.0, 54364.0, 40.589, 83.537),
}
TARGETS = {
    506: "Himalia",
    507: "Elara",
    508: "Pasiphae",
    509: "Sinope",
    510: "Lysithea",
    511: "Carme",
    512: "Ananke",
    513: "Leda",
    609: "Phoebe",
}
J0 = 2458849.5  # JD TDB at 2020-01-01 00:00
EPOCH = datetime(2020, 1, 1)
DIRECTORY = Path(__file__).resolve().parent


def fetch(cache, name, params):
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / (name + ".txt")
    url = API + "?" + urlencode({"format": "json", **params})
    if not path.exists():
        for attempt in range(4):
            try:
                with urlopen(url, timeout=120) as response:
                    raw = json.load(response)["result"]
                if "$$SOE" not in raw:
                    raise ValueError(raw[-800:])
                path.write_text(raw)
                break
            except (OSError, ValueError, KeyError):
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    data = path.read_text()
    head, data = data.split("$$SOE")
    data = data.split("$$EOE")[0].strip()
    receipt = {
        "id": name,
        "url": url,
        "target": next(
            v.strip() for v in head.splitlines() if v.startswith("Target body name")
        ),
        "sha256_rows": hashlib.sha256(data.encode()).hexdigest(),
    }
    return data, receipt


def state_table(cache, target, start_year=2020, stop_year=2041):
    start = datetime(start_year, 1, 1) - timedelta(days=5)
    stop = datetime(stop_year, 1, 1) + timedelta(days=5)
    key = (
        str(target)
        if (start_year, stop_year) == (2020, 2041)
        else f"{target}-{start_year}-{stop_year}"
    )
    rows, receipt = fetch(
        cache,
        key,
        dict(
            COMMAND=f"'{target}'",
            CENTER="'500@0'",
            EPHEM_TYPE="VECTORS",
            START_TIME=start.strftime("'%Y-%b-%d'"),
            STOP_TIME=stop.strftime("'%Y-%b-%d'"),
            STEP_SIZE="'1 d'",
            OUT_UNITS="KM-S",
            VEC_TABLE="2",
            CSV_FORMAT="YES",
            REF_PLANE="FRAME",
            OBJ_DATA="NO",
        ),
    )
    data = np.array(
        [
            [float(r[i]) for i in (0, 2, 3, 4, 5, 6, 7)]
            for r in csv.reader(io.StringIO(rows))
        ]
    )
    return (
        CubicHermiteSpline(data[:, 0] - J0, data[:, 1:4], data[:, 4:7] * DAY),
        receipt,
    )


def norm(x):
    return np.linalg.norm(x, axis=-1)


def unit(x):
    return x / norm(x)[..., None]


def seen_from(point, target, t):
    """Target at emission for a photon arriving at point at barycentric t."""
    tm = t
    for _ in range(3):
        tm = t - norm(target(tm) - point) / (C * DAY)
    return target(tm) - point


def shadow(states, moon, planet, t):
    m = states[moon](t)
    p = seen_from(m, states[planet], t)
    s = seen_from(m, states[10], t)
    # Sky plane centered on planet as seen from the moon. Solar center offset.
    axis = unit(p)
    direction = unit(s)
    delta = direction - np.sum(direction * axis, axis=-1)[..., None] * axis
    return delta, axis, norm(p), norm(s), np.sum(direction * axis, axis=-1)


def spherical_gap(states, moon, planet, t):
    delta, axis, dp, ds, dot = shadow(states, moon, planet, t)
    return (
        np.arctan2(norm(delta), dot)
        - np.arcsin(PLANETS[planet][0] / dp)
        - np.arcsin(SUN_R / ds)
    )


def outline(states, moon, planet, t):
    """Orthographic shadow cross-section in km in plane at planet center."""
    a, b, ra, dec = PLANETS[planet]
    ra, dec = np.radians([ra, dec])
    pole = np.array([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)])
    delta, axis, dp, ds, dot = shadow(states, moon, planet, t)
    north = unit(pole - np.dot(axis, pole) * axis)
    east = np.cross(north, axis)
    minor = np.sqrt(b * b + (a * a - b * b) * np.dot(axis, pole) ** 2)
    # Tangent-plane projection (small angles, error <~3 km for this cohort).
    offset = delta * dp / dot
    return a, minor, np.dot(offset, east), np.dot(offset, north), dp * SUN_R / ds


def signed_distance(a, b, x, y):
    """Signed distance from point to the boundary of a filled ellipse."""
    x, y = abs(x), abs(y)

    def square(theta):
        return (a * np.cos(theta) - x) ** 2 + (b * np.sin(theta) - y) ** 2

    fit = minimize_scalar(square, bounds=(0, np.pi / 2), method="bounded")
    d = np.sqrt(min(square(0), square(np.pi / 2), fit.fun))
    return -d if (x / a) ** 2 + (y / b) ** 2 < 1 else d


def metrics(states, moon, planet, t):
    a, b, x, y, r = outline(states, moon, planet, t)
    d = signed_distance(a, b, x, y)
    return d - r, d + r


def fraction(states, moon, planet, t):
    """Occulted area of a uniform solar disc using vertical strips."""
    from scipy.integrate import quad

    a, b, x, y, r = outline(states, moon, planet, t)
    if signed_distance(a, b, x, y) <= -r:
        return 1.0

    def overlap(v):
        h = r * np.sqrt(max(0, 1 - v * v))
        k = b * np.sqrt(max(0, 1 - ((x + v * r) / a) ** 2))
        return max(0, min(k, y + h) - max(-k, y - h))

    lo, hi = max(-1, (-a - x) / r), min(1, (a - x) / r)
    if lo >= hi:
        return 0.0
    # Break into strips so adaptive quadrature cannot miss a narrow overlap.
    return quad(
        overlap, lo, hi, points=np.linspace(lo, hi, 33)[1:-1], limit=200, epsabs=1e-3
    )[0] / (np.pi * r)


def earth_time(states, moon, t):
    pt = states[moon](t)
    tt = t
    for _ in range(3):
        tt = t + norm(states[399](tt) - pt) / (C * DAY)
    return tt


def utc(t):
    # Effective UTC dates for leap seconds spanning this investigation.
    # Predicted UTC holds TAI-UTC at 37 seconds after 2017. |TDB-TT| < 2 ms.
    epoch = EPOCH + timedelta(days=float(t))
    tai_utc = 31
    for date, offset in [
        (datetime(1999, 1, 1), 32),
        (datetime(2006, 1, 1), 33),
        (datetime(2009, 1, 1), 34),
        (datetime(2012, 7, 1), 35),
        (datetime(2015, 7, 1), 36),
        (datetime(2017, 1, 1), 37),
    ]:
        if epoch >= date + timedelta(seconds=offset + 32.184):
            tai_utc = offset
    return (epoch - timedelta(seconds=tai_utc + 32.184)).strftime("%Y-%m-%d %H:%M:%S")


def angle(u, v):
    return np.degrees(np.arctan2(norm(np.cross(u, v)), np.dot(u, v)))


def run(cache, outdir, start_year=2020, stop_year=2041, step_hours=3.0):
    if not 1998 <= start_year < stop_year:
        raise ValueError(
            "Choose increasing years within the UTC conversion table (1998 onward)."
        )
    states, receipts = {}, []
    for target in [10, 399, 599, 699, *TARGETS]:
        states[target], r = state_table(cache, target, start_year, stop_year)
        receipts.append(r)
    step = step_hours / 24
    grid = np.arange(
        (datetime(start_year, 1, 1) - EPOCH).days,
        (datetime(stop_year, 1, 1) - EPOCH).days + step / 2,
        step,
    )
    events = []
    for moon, name in TARGETS.items():
        planet = 699 if moon == 609 else 599
        gaps = spherical_gap(states, moon, planet, grid)
        candidates = (
            np.flatnonzero((gaps[1:-1] < gaps[:-2]) & (gaps[1:-1] < gaps[2:])) + 1
        )
        for i in candidates:
            t = grid[i]
            if gaps[i] > 0.05:
                continue
            # Minimize the separation of the projected limb and the finite Sun in km.
            f = lambda tt: metrics(states, moon, planet, tt)[0]
            fit = minimize_scalar(
                f,
                bounds=(t - 0.20, t + 0.20),
                method="bounded",
                options={"xatol": 1e-9},
            )
            if fit.fun >= 0:
                continue
            mid = fit.x

            def contacts(fun):
                step = 0.035
                while fun(mid - step) < 0 or fun(mid + step) < 0:
                    step *= 2
                return brentq(fun, mid - step, mid), brentq(fun, mid, mid + step)

            start, end = contacts(f)
            full = lambda tt: metrics(states, moon, planet, tt)[1]
            full_start, full_end = contacts(full) if full(mid) < 0 else (None, None)
            # For partial eclipses maximize the actual occulted area.
            if full_start is None:
                mid = minimize_scalar(
                    lambda tt: -fraction(states, moon, planet, tt),
                    bounds=(start, end),
                    method="bounded",
                ).x
            arrival = earth_time(states, moon, mid)
            e = states[399](arrival)
            sat = states[moon](mid) - e
            p = seen_from(e, states[planet], arrival)
            sun = seen_from(e, states[10], arrival)
            seen = lambda v: utc(earth_time(states, moon, v)) if v is not None else ""
            events.append(
                dict(
                    target=name,
                    id=moon,
                    tdb_days_from_2020=round(mid, 8),
                    start_utc=seen(start),
                    full_start_utc=seen(full_start),
                    mid_utc=seen(mid),
                    full_end_utc=seen(full_end),
                    end_utc=seen(end),
                    max_blocked=round(fraction(states, moon, planet, mid), 6),
                    duration_min=round((end - start) * 1440, 2),
                    totality_min=(
                        round((full_end - full_start) * 1440, 2)
                        if full_start is not None
                        else 0.0
                    ),
                    elong_deg=round(angle(sat, sun), 2),
                    planet_sep_arcsec=round(angle(sat, p) * 3600, 1),
                )
            )
    outdir.mkdir(exist_ok=True, parents=True)
    with (outdir / "events.csv").open("w") as f:
        fields = (
            "target",
            "id",
            "tdb_days_from_2020",
            "start_utc",
            "full_start_utc",
            "mid_utc",
            "full_end_utc",
            "end_utc",
            "max_blocked",
            "duration_min",
            "totality_min",
            "elong_deg",
            "planet_sep_arcsec",
        )
        w = csv.DictWriter(f, fields, lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(events, key=lambda v: v["mid_utc"]))
    (outdir / "sources.json").write_text(json.dumps(receipts, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=DIRECTORY)
    p.add_argument("--start-year", type=int, default=2020)
    p.add_argument(
        "--stop-year",
        type=int,
        default=2041,
        help="Exclusive. Leap-second table starts in 1998.",
    )
    p.add_argument("--step-hours", type=float, choices=[1.5, 3.0], default=3.0)
    args = p.parse_args()
    run(args.cache, args.outdir, args.start_year, args.stop_year, args.step_hours)

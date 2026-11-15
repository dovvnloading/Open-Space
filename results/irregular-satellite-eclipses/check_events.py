"""Check every event by observer-table flags and test the sensitivity of contacts."""

import argparse
import csv
from datetime import datetime, timedelta
import io
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import search

# East longitude, latitude (degrees), height (km). A site test uses airless altitude.
SITES = {
    "Cerro Pachon": (-70.7366, -30.2407, 2.72),
    "La Palma": (-17.89, 28.76, 2.33),
    "Maunakea": (-155.4783, 19.8283, 4.20),
    "Siding Spring": (149.0661, -31.2733, 1.17),
    "Sutherland": (20.81, -32.38, 1.80),
    "Palomar": (-116.86, 33.36, 1.71),
    "La Silla": (-70.733, -29.258, 2.40),
    "Xinglong": (117.58, 40.40, 0.96),
    "Kottamia": (31.83, 29.93, 0.47),
}


def parse(s):
    return datetime.fromisoformat(s)


def site_visibility(states, result):
    windows = [
        ("Phoebe ring ingress", "2032-03-01 22:05", 609),
        ("Phoebe ring egress", "2032-03-03 10:30", 609),
        ("Sinope start of totality", "2037-10-17 15:01", 509),
        ("Sinope end of totality", "2037-10-18 00:28", 509),
    ]
    for label, dt, moon in windows:
        jd = search.J0 + (parse(dt) - search.EPOCH).total_seconds() / search.DAY
        alt = []
        for place, coords in SITES.items():
            # Altitude from celestial vector at the geocenter, rotating geodetic zenith.
            # GMST IAU approximated by 280.46 + 360.98564737 d; precession to date below.
            t = (jd - 2451545) / 36525
            zeta, z, theta = np.deg2rad(
                np.array(
                    [
                        2306.2181 * t + 0.30188 * t * t,
                        2306.2181 * t + 1.09468 * t * t,
                        2004.3109 * t - 0.42665 * t * t,
                    ]
                )
                / 3600
            )
            theta_sid = np.deg2rad(
                (
                    280.46061837
                    + 360.98564736629 * (jd - 2451545)
                    + 0.000387933 * t * t
                    + coords[0]
                )
                % 360
            )
            when = jd - search.J0 + 69.184 / search.DAY
            earth = states[399](when)

            def altitude(target):
                x, y, zz = search.unit(search.seen_from(earth, states[target], when))
                ra = np.arctan2(y, x)
                dec = np.arcsin(zz)
                ax = np.cos(dec) * np.sin(ra + zeta)
                bx = np.cos(theta) * np.cos(dec) * np.cos(ra + zeta) - np.sin(
                    theta
                ) * np.sin(dec)
                cx = np.sin(theta) * np.cos(dec) * np.cos(ra + zeta) + np.cos(
                    theta
                ) * np.sin(dec)
                hour = theta_sid - (np.arctan2(ax, bx) + z)
                lat = np.deg2rad(coords[1])
                return np.rad2deg(
                    np.arcsin(
                        np.sin(lat) * cx
                        + np.cos(lat) * np.sqrt(1 - cx * cx) * np.cos(hour)
                    )
                )

            alt.append(
                {
                    "site": place,
                    "target_deg": round(altitude(moon), 1),
                    "sun_deg": round(altitude(10), 1),
                }
            )
        result[label] = {"utc": dt, "altitudes": alt}


def run(cache, outdir):
    result = {"sites": {}, "events": {}}
    receipts = []
    states = {i: search.state_table(cache, i)[0] for i in [10, 399, 599, 699, 509, 609]}
    site_visibility(states, result["sites"])
    for row in csv.DictReader((outdir / "events.csv").open()):
        name = row["target"]
        moon = int(row["id"])
        planet = 699 if moon == 609 else 599
        mid = parse(row["mid_utc"])
        beg = mid - timedelta(hours=30)
        end = mid + timedelta(hours=30)
        data, receipt = search.fetch(
            cache,
            name + "-observer",
            dict(
                COMMAND=f"'{moon}'",
                CENTER="'500@399'",
                EPHEM_TYPE="OBSERVER",
                START_TIME=beg.strftime("'%Y-%b-%d %H:%M'"),
                STOP_TIME=end.strftime("'%Y-%b-%d %H:%M'"),
                STEP_SIZE="'1 m'",
                QUANTITIES="'1,9,12,23'",
                CSV_FORMAT="YES",
                OBJ_DATA="NO",
            ),
        )
        receipts.append(receipt)
        obs = list(csv.reader(io.StringIO(data)))
        codes = [r[8].strip() for r in obs]
        group = [r for r in obs if r[8].strip() == "u"]
        nearest = min(
            obs,
            key=lambda r: abs(datetime.strptime(r[0].strip(), "%Y-%b-%d %H:%M") - mid),
        )
        tm = float(row["tdb_days_from_2020"])

        def contacts(pad):
            func = lambda t: search.metrics(states, moon, planet, t)[1] - pad
            return [brentq(func, tm - 1, tm), brentq(func, tm, tm + 1)]

        original = np.array(contacts(0))
        # Isotropic signed limb-distance offset approximates projected position or limb error.
        margin = max(
            np.max(np.abs(np.array(contacts(pad)) - original)) for pad in [-1000, 1000]
        )
        result["events"][name] = {
            "observer_status_counts": {k: codes.count(k) for k in sorted(set(codes))},
            "first_full_minute_utc": group[0][0].strip(),
            "last_full_minute_utc": group[-1][0].strip(),
            "uneclipsed_model_mag": float(nearest[5]),
            "separation_arcsec": float(nearest[7]),
            "solar_elongation_deg": float(nearest[9]),
            "limb_clearance_during_totality_km": round(
                -search.metrics(states, moon, planet, tm)[1], 1
            ),
            "max_totality_contact_shift_1000km_min": round(margin * 1440, 2),
        }
    (outdir / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
    (outdir / "observer-sources.json").write_text(json.dumps(receipts, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=search.DIRECTORY)
    args = p.parse_args()
    run(args.cache, args.outdir)

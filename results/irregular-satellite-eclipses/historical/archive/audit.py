"""Rebuild Phoebe's May 2016 coverage search from public archive responses."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timedelta
from hashlib import sha256
from io import StringIO
import json
from math import acos, atan2, cos, degrees, radians, sin
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
EPOCH = datetime(2016, 5, 19)
MJD0 = datetime(1858, 11, 17)
LO, HI = datetime(2016, 5, 20), datetime(2016, 5, 26)
START, STOP = datetime(2016, 5, 22, 23, 39, 41), datetime(2016, 5, 23, 10, 59, 25)
LIGHT_SECONDS_AU = 499.004784
ISO = datetime.fromisoformat


def requests():
    horizons = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = dict(
        format="json",
        COMMAND="'609'",
        OBJ_DATA="NO",
        MAKE_EPHEM="YES",
        EPHEM_TYPE="OBSERVER",
        START_TIME="'2016-05-19'",
        STOP_TIME="'2016-05-27'",
        STEP_SIZE="'10m'",
        QUANTITIES="'1,20'",
        ANG_FORMAT="DEG",
        CSV_FORMAT="YES",
        EXTRA_PREC="YES",
    )

    def tap(query):
        return dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)

    services = {
        "cas.json": (horizons, dict(params, CENTER="'500@-82'")),
        "earth.json": (horizons, dict(params, CENTER="'500@399'")),
        "saturn-earth.json": (
            horizons,
            dict(params, COMMAND="'699'", CENTER="'500@399'", QUANTITIES="'1,13,20'"),
        ),
        "opus.json": (
            "https://opus.pds-rings.seti.org/opus/api/data.json",
            dict(
                time1="2016-05-20",
                time2="2016-05-26",
                limit=2000,
                cols="opusid,instrument,target,time1,time2,rightasc1,rightasc2,declination1,declination2,COISScamera",
            ),
        ),
        "cadca.out": (
            "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl",
            dict(
                search="bynameHorizons",
                object="609",
                epoch1="2016-05-20",
                epoch2="2016-05-26",
                eellipse="10",
                eunits="arcseconds",
                format="tsv",
                lang="en",
            ),
        ),
        "skymap.csv": (
            "https://api.skymapper.nci.org.au/public/tap/sync",
            tap(
                "SELECT image_id,date,exp_time,ra,decl,filter,object FROM dr4.images WHERE date BETWEEN 57528 AND 57534"
            ),
        ),
        "eso_raw.csv": (
            "https://archive.eso.org/tap_obs/sync",
            tap(
                "SELECT dp_id,instrument,object,mjd_obs,exposure,ra,dec,ra_pnt,dec_pnt FROM dbo.raw "
                "WHERE mjd_obs BETWEEN 57528 AND 57534 AND (INTERSECTS(s_region,CIRCLE('ICRS',252.55,-20.64,0.5))=1 "
                "OR (ra BETWEEN 249 AND 256 AND dec BETWEEN -24 AND -17))"
            ),
        ),
    }
    return {name: url + "?" + urlencode(p) for name, (url, p) in services.items()}


def hours(dt):
    return (dt - EPOCH).total_seconds() / 3600


def interpolator(raw, columns=(3, 4, 5)):
    data = json.loads(raw)["result"].split("$$SOE")[1].split("$$EOE")[0].strip()
    rows = list(csv.reader(data.splitlines()))
    t = [hours(datetime.strptime(r[0].strip(), "%Y-%b-%d %H:%M")) for r in rows]
    v = np.array([[float(r[i]) for i in columns] for r in rows]).T

    def interpolate(dt):
        return [float(np.interp(hours(dt), t, col)) for col in v]

    return interpolate


def to_spacecraft(t, earth, spacecraft):
    # Both range products are observer-to-emission-point light-time distances.
    emission = t - timedelta(seconds=earth(t)[2] * LIGHT_SECONDS_AU)
    local = t
    for _ in range(3):
        local = emission + timedelta(seconds=spacecraft(local)[2] * LIGHT_SECONDS_AU)
    return local


def distance(ra, dec, ra2, dec2):
    da, d1, d2 = map(radians, [ra - ra2, dec, dec2])
    return degrees(
        acos(max(-1.0, min(1.0, sin(d1) * sin(d2) + cos(d1) * cos(d2) * cos(da))))
    )


def box_distance(ra, dec, bounds):
    a, b, c, d = map(float, bounds)
    # Choose the nearest direction in the spherical RA/Dec bounding rectangle.
    candidates = [a, b]
    if (a <= ra <= b) or (a > b and (ra >= a or ra <= b)):
        candidates.append(ra)
    best = 180.0
    for lon in candidates:
        lat = degrees(
            atan2(sin(radians(dec)), cos(radians(dec)) * cos(radians(lon - ra)))
        )
        best = min(best, distance(ra, dec, lon, max(c, min(d, lat))))
    return best


def ring_intercept(t, earth, saturn):
    """Geocentric astrometric sky rays, accounting for light time in both target vectors."""

    def vec(point):
        ra, dec, r = point
        a, b = radians(ra), radians(dec)
        return np.array([cos(a) * cos(b), sin(a) * cos(b), sin(b)]) * r * 149597870.7

    y = (t - datetime(2000, 1, 1, 12)).total_seconds() / (36525 * 86400)
    pole = vec((40.589 - 0.036 * y, 83.537 - 0.004 * y, 1 / 149597870.7))
    phoebe, sat = vec(earth(t)), vec(saturn(t))
    ratio = np.dot(pole, sat) / np.dot(pole, phoebe)
    offset = phoebe * ratio - sat
    ra, dec, _ = earth(t)
    ra2, dec2, _ = saturn(t)
    return dict(
        utc=t.isoformat(),
        separation_arcsec=round(3600 * distance(ra, dec, ra2, dec2), 3),
        intercept_fraction=round(ratio, 8),
        ring_radius_km=round(float(np.linalg.norm(offset)), 0),
        phoebe_farther_by_km=round(
            float(np.linalg.norm(phoebe) - np.linalg.norm(sat)), 0
        ),
    )


def in_event(t, seconds, start=START, end=STOP):
    return t <= end and t + timedelta(seconds=float(seconds)) >= start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache", type=Path, required=True, help="Directory for public API responses"
    )
    parser.add_argument("--output", type=Path, default=HERE / "coverage.json")
    opts = parser.parse_args()
    opts.cache.mkdir(parents=True, exist_ok=True)
    blobs, manifest = {}, {}
    for name, url in requests().items():
        path = opts.cache / name
        if not path.exists():
            with urlopen(url, timeout=120) as response:
                raw = response.read()
            path.write_bytes(raw)
        raw = path.read_bytes()
        blobs[name] = raw.decode()
        manifest[name] = {"url": url, "sha256": sha256(raw).hexdigest()}
    cas = interpolator(blobs["cas.json"])
    earth = interpolator(blobs["earth.json"])
    saturn = interpolator(blobs["saturn-earth.json"], columns=(3, 4, 6))
    alpo_times = [
        "2016-05-22T13:50:24",
        "2016-05-22T14:09:06",
        "2016-05-22T14:16:46",
        "2016-05-22T15:25:58",
        "2016-05-22T15:30:11",
        "2016-05-23T16:02:30",
        "2016-05-23T16:12:15",
        "2016-05-23T16:29:53",
    ]
    first, last = [to_spacecraft(t, earth, cas) for t in (START, STOP)]
    opus = json.loads(blobs["opus.json"])
    rows = opus["page"]
    if opus["available"] != len(rows):
        raise ValueError("Truncated OPUS response")
    counts = Counter(r[1] for r in rows)
    overlaps = [
        r
        for r in rows
        if r[1].startswith("Cassini") and ISO(r[3]) <= last and ISO(r[4]) >= first
    ]
    iss = []
    for r in rows:
        if r[1] != "Cassini ISS":
            continue
        t = ISO(r[3]) + (ISO(r[4]) - ISO(r[3])) / 2
        ra, dec, _ = cas(t)
        sep = box_distance(ra, dec, r[5:9])
        iss.append(
            dict(
                id=r[0],
                start=r[3],
                end=r[4],
                target=r[2],
                camera=r[9],
                phoebe_ra_deg=round(ra, 6),
                phoebe_dec_deg=round(dec, 6),
                bounds_deg=r[5:9],
                gap_deg=round(sep, 4),
            )
        )
    sky = list(csv.DictReader(StringIO(blobs["skymap.csv"])))
    sky_gaps = []
    for r in sky:
        t = MJD0 + timedelta(days=float(r["date"]), seconds=float(r["exp_time"]) / 2)
        ra, dec, _ = earth(t)
        sep = distance(ra, dec, float(r["ra"]), float(r["decl"]))
        sky_gaps.append(
            (
                round(sep, 4),
                r["image_id"],
                t,
                in_event(
                    t - timedelta(seconds=float(r["exp_time"]) / 2), r["exp_time"]
                ),
            )
        )
    cadc = list(csv.DictReader(StringIO(blobs["cadca.out"]), delimiter="\t"))
    cadc = [(r, MJD0 + timedelta(days=float(r["MJD"]))) for r in cadc]
    cadc = [(r, t) for r, t in cadc if LO <= t < HI]
    esoraw = list(csv.DictReader(StringIO(blobs["eso_raw.csv"])))
    uv = []
    for r in overlaps:
        t = ISO(r[3]) + (ISO(r[4]) - ISO(r[3])) / 2
        ra, dec, _ = cas(t)
        uv.append(
            dict(
                id=r[0],
                target=r[2],
                start=r[3],
                end=r[4],
                bounds_gap_deg=round(box_distance(ra, dec, r[5:9]), 4),
            )
        )
    result = {
        "scope_utc": [LO, HI],
        "earth_disk_eclipse_utc": [START, STOP],
        "cassini_disk_eclipse_utc": [first, last],
        "opus_counts": dict(sorted(counts.items())),
        "opus_cassini_event_products": uv,
        "cassini_iss_minimum_gap": min(iss, key=lambda x: x["gap_deg"]),
        "iss_last_before": max(
            (r for r in iss if ISO(r["end"]) < first), key=lambda x: x["end"]
        ),
        "iss_first_after": min(
            (r for r in iss if ISO(r["start"]) > last), key=lambda x: x["start"]
        ),
        "cadc_products_by_date_and_instrument": dict(
            sorted(
                Counter(
                    str(t.date()) + " " + r["Telescope/Instrument"] for r, t in cadc
                ).items()
            )
        ),
        "cadc_distinct_mjd": len({r["MJD"] for r, t in cadc}),
        "cadc_event_products": sum(in_event(t, r["Exptime"] or 0) for r, t in cadc),
        "skymapper_count": len(sky),
        "skymapper_event_count": sum(r[3] for r in sky_gaps),
        "skymapper_minimum_center_gap_deg": min(r[0] for r in sky_gaps),
        "skymapper_event_minimum_center_gap_deg": min(r[0] for r in sky_gaps if r[3]),
        "eso_raw_preselection_count": len(esoraw),
        "eso_raw_event_products": sum(
            in_event(MJD0 + timedelta(days=float(r["mjd_obs"])), r["exposure"] or 0)
            for r in esoraw
        ),
        "alpo_geocentric_geometry": [
            ring_intercept(ISO(t), earth, saturn) for t in alpo_times
        ],
        "sources": manifest,
    }
    opts.output.write_text(
        json.dumps(result, indent=2, default=lambda d: d.isoformat(timespec="seconds"))
        + "\n"
    )


if __name__ == "__main__":
    main()

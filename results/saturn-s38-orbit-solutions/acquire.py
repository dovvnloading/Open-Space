"""Re-acquire the measurement tables and archive geometry from their primary sources."""

import csv
import gzip
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent / "inputs"
MPC = "https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/SatObs.txt.gz"
CIRCULAR = "https://minorplanetcenter.net/mpec/K25/K25EF4.html"
API = "https://ssd.jpl.nasa.gov/api/horizons.api"
EPOCH = 2460200.5
SEC_PER_DAY = 86400
TARGET = 38
# All objects from S30–S41 with measurements in this circular.
MOONS = (31, 35, 36, 37, 38, 39, 40, 41)


def read_url(url):
    for attempt in range(3):
        try:
            with urlopen(url, timeout=100) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def write_table(name, header, rows):
    with (ROOT / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def measurement(line):
    year, month = int(line[15:19]), int(line[20:22])
    day = float(line[23:32])
    midnight = datetime(year, month, int(day), tzinfo=timezone.utc)
    jd = midnight.timestamp() / SEC_PER_DAY + 2440587.5 + day - int(day)
    h, m, s = map(float, line[32:44].split())
    deg, arcmin, arcsec = map(float, line[45:56].split())
    ra = 15 * (h + m / 60 + s / 3600)
    dec = (-1 if line[44] == "-" else 1) * (deg + arcmin / 60 + arcsec / 3600)
    return jd, ra, dec


def horizons(params):
    encoded = urlencode(dict(format="json", **params))
    answer = json.loads(read_url(API + "?" + encoded))
    text = answer["result"]
    if "$$SOE" not in text:
        raise ValueError(text[-500:])
    rows = text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    return text, [row.split(",") for row in rows]


def jd_list(epochs):
    return "'" + ",".join(f"{value:.8f}" for value in epochs) + "'"


def vector_data(name, target, center, times):
    epochs = list(times)
    rows = []
    for offset in range(0, len(epochs), 36):
        _, chunk = horizons(dict(
            COMMAND=f"'{target}'", CENTER=f"'{center}'", EPHEM_TYPE="VECTORS",
            REF_PLANE="FRAME", REF_SYSTEM="ICRF", VEC_TABLE="2", CSV_FORMAT="YES",
            OUT_UNITS="KM-S", TLIST=jd_list(epochs[offset:offset+36]),
        ))
        rows.extend(chunk)
    write_table(name + ".csv", ["jd_tdb", "x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s"],
                [[row[0].strip(), *[x.strip() for x in row[2:8]]] for row in rows])


def main():
    ROOT.mkdir(exist_ok=True)
    circular = read_url(CIRCULAR).decode().splitlines()
    bulk = gzip.decompress(read_url(MPC)).decode().splitlines()
    all_observations = []
    for moon in MOONS:
        key = f"SK23S{moon:02d}0"
        rows = [line for line in circular if line[4:12] == key]
        if not rows or any(line[77:80] != "T14" for line in rows):
            raise ValueError(f"Unexpected discovery observations for {moon}")
        observed = [measurement(row) for row in rows]
        if moon == TARGET:
            archive = [measurement(row) for row in bulk if row[4:12] == key]
            if sorted(archive) != sorted(observed):
                raise ValueError("Public observation inventory changed. Reassess the full arc.")
            target_epochs = [row[0] for row in observed]
        text, model = horizons(dict(
            COMMAND=f"'{65235+moon}'", CENTER="'T14'", EPHEM_TYPE="OBSERVER",
            CSV_FORMAT="YES", QUANTITIES="'1'", ANG_FORMAT="DEG", EXTRA_PREC="YES",
            TLIST=jd_list(row[0] for row in observed),
        ))
        if "sat455_merged_DE440" not in text:
            raise ValueError("Ephemeris release changed. Review before updating the frozen data.")
        if len(model) != len(observed):
            raise ValueError("Horizons supplied an incomplete epoch list")
        for values, predicted in zip(observed, model):
            all_observations.append([moon, *values, predicted[3].strip(), predicted[4].strip()])
    write_table("astrometry.csv", ["s2023_s", "jd_utc", "obs_ra_deg", "obs_dec_deg",
                                   "sat455_ra_deg", "sat455_dec_deg"], all_observations)
    times = range(2460160, 2460251, 2)
    # TDB−UTC changes by under one millisecond over this arc.
    vector_data("ssb_from_t14", 0, "T14", [jd + 69.182 / SEC_PER_DAY for jd in target_epochs])
    vector_data("saturn_barycenter", 6, "500@0", times)
    vector_data("sun", 10, "500@6", times)
    vector_data("jupiter_barycenter", 5, "500@6", times)
    vector_data("sat455_state", 65273, "500@6", [EPOCH])

    preview = [t + 7 / 24 for t in [2460264.5, 2460293.5, 2460537.5, 2461375.5]]
    vector_data("prediction_observer", 0, "T14", [t + 69.182 / SEC_PER_DAY for t in preview])
    grid = [2460180.5 + 12 * k for k in range(102)]
    vector_data("prediction_saturn", 6, "500@0", grid)
    vector_data("prediction_sun", 10, "500@6", grid)
    vector_data("prediction_jupiter", 5, "500@6", grid)
    _, rows = horizons(dict(
        COMMAND="'65273'", CENTER="'T14'", EPHEM_TYPE="OBSERVER", CSV_FORMAT="YES",
        QUANTITIES="'1'", ANG_FORMAT="DEG", EXTRA_PREC="YES", TLIST=jd_list(preview),
    ))
    write_table("prediction_sat455.csv", ["jd_utc", "ra_deg", "dec_deg"],
                [(jd, row[3].strip(), row[4].strip()) for jd, row in zip(preview, rows)])


if __name__ == "__main__":
    main()

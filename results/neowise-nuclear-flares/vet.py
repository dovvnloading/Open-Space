"""Look up catalogued objects and reported transients at each flare position.

Usage: python vet.py FLARES.csv POSITIONS.csv OUT.json [--radius ARCSEC]

For every row of FLARES.csv (column `id`), the position from POSITIONS.csv is
searched in SIMBAD (object identifiers and types) and in the Transient Name
Server (reported and classified transients). Results are cached in OUT.json so
the command can be rerun as the flare list grows.
"""
import argparse
import csv
import io
import json
import pathlib
import time

import requests

SIMBAD = "https://simbad.harvard.edu/simbad/sim-tap/sync"
TNS = "https://www.wis-tns.org/search"
AGENT = {"User-Agent": "Mozilla/5.0 (compatible; open-space-research/0.2)"}


def simbad(ra, dec, radius):
    query = ("SELECT main_id, otype, ra, dec, rvz_redshift, "
             f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) * 3600 AS sep "
             f"FROM basic WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius / 3600})) = 1 "
             "ORDER BY sep")
    response = requests.post(SIMBAD, data={"request": "doQuery", "lang": "adql", "format": "csv", "query": query},
                             timeout=120)
    response.raise_for_status()
    return [{"name": r["main_id"], "type": r["otype"], "sep": round(float(r["sep"]), 2), "z": r["rvz_redshift"]}
            for r in csv.DictReader(io.StringIO(response.text))]


def tns(ra, dec, radius):
    """Cone search; the service allows two per minute and reports the reset time."""
    for attempt in range(8):
        response = requests.get(TNS, params={"ra": ra, "decl": dec, "radius": radius, "coords_unit": "arcsec",
                                             "format": "csv", "num_page": 50}, headers=AGENT, timeout=120)
        reset = response.headers.get("x-cone-rate-limit-reset", "30")
        wait = int(reset) + 2 if reset.isdigit() else 32
        if response.status_code == 200 and response.text.startswith('"ID"'):
            if response.headers.get("x-cone-rate-limit-remaining", "0") in ("0", "Exceeded"):
                time.sleep(wait)
            return [{"name": r["Name"], "type": r["Obj. Type"], "discovered": r["Discovery Date (UT)"],
                     "internal": r["Disc. Internal Name"], "group": r["Reporting Group/s"],
                     "ra": r["RA"], "dec": r["DEC"]}
                    for r in csv.DictReader(io.StringIO(response.text))]
        time.sleep(wait)
    raise RuntimeError(f"TNS search failed with status {response.status_code}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("flares")
    parser.add_argument("positions")
    parser.add_argument("out")
    parser.add_argument("--radius", type=float, default=10.0)
    args = parser.parse_args()
    positions = {r["id"]: (float(r["ra"]), float(r["dec"])) for r in csv.DictReader(open(args.positions))}
    out = pathlib.Path(args.out)
    cache = json.loads(out.read_text()) if out.exists() else {}
    for row in csv.DictReader(open(args.flares)):
        name = row["id"]
        if name in cache:
            continue
        ra, dec = positions[name]
        try:
            cache[name] = {"ra": ra, "dec": dec, "simbad": simbad(ra, dec, args.radius),
                           "tns": tns(ra, dec, args.radius)}
        except Exception as error:
            print(name, "FAILED", error, flush=True)
            continue
        out.write_text(json.dumps(cache, indent=1))
        entry = cache[name]
        print(name, [(o["name"], o["type"]) for o in entry["simbad"][:3]],
              [(o["name"], o["type"], o["discovered"][:10]) for o in entry["tns"]], flush=True)


if __name__ == "__main__":
    main()

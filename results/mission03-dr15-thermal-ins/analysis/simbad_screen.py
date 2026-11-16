"""Find published SIMBAD objects within 15 arcsec of the catalogue shortlist.

Small overlapping cone unions retain a query that users can inspect. Empty results
mean no returned index entry inside the cone, not an unidentified astrophysical source.
"""

from argparse import ArgumentParser
import csv
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor

SERVICE = "https://simbad.harvard.edu/simbad/sim-tap/sync"
RADIUS = 15.0


def separation(r1, d1, r2, d2):
    dra = math.radians(r1 - r2)
    ddec = math.radians(d1 - d2)
    a = (
        math.sin(ddec / 2) ** 2
        + math.cos(math.radians(d1))
        * math.cos(math.radians(d2))
        * math.sin(dra / 2) ** 2
    )
    return 7200 * math.degrees(math.asin(min(1, math.sqrt(a))))


def match_group(group):
    test = "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{r}))"
    clause = " OR ".join(
        test.format(ra=row["RA"], dec=row["DEC"], r=RADIUS / 3600.0) for row in group
    )
    query = "SELECT main_id, otype, ra, dec, nbref FROM basic WHERE " + clause
    params = urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=10000, QUERY=query)
    )
    with urlopen(SERVICE + "?" + params, timeout=60) as response:
        reader = csv.DictReader(io.StringIO(response.read().decode()))
        if reader.fieldnames != ["main_id", "otype", "ra", "dec", "nbref"]:
            raise ValueError(f"Unexpected SIMBAD schema: {reader.fieldnames}")
        matches = list(reader)
    output = []
    for row in group:
        nearby = []
        for m in matches:
            d = separation(
                float(row["RA"]), float(row["DEC"]), float(m["ra"]), float(m["dec"])
            )
            if d <= RADIUS:
                nearby.append(dict(separation_arcsec=round(d, 4), **m))
        output.append(
            dict(
                SRCID=row["SRCID"],
                IAUNAME=row["IAUNAME"],
                matches=sorted(nearby, key=lambda x: x["separation_arcsec"]),
            )
        )
    return output


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("shortlist", type=Path, nargs="+")
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    rows = []
    for path in args.shortlist:
        with path.open() as f:
            rows.extend(csv.DictReader(f))
    groups = [rows[s : s + 12] for s in range(0, len(rows), 12)]
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = [row for group in pool.map(match_group, groups) for row in group]
    args.output.write_text(
        json.dumps(
            dict(service=SERVICE, search_radius_arcsec=RADIUS, objects=results),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

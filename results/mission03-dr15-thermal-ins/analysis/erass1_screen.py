"""Record eRASS1 main catalogue neighbors within 30 arcsec of 5XMM entries.

The positional search radius is for retrieval. A returned eROSITA source can
have multiple XMM contributors. A missing entry is not a flux upper limit.
"""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree

import numpy as np
from astropy.io.votable import parse_single_table

from simbad_screen import separation

SERVICE = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
RADIUS = 30.0
FIELDS = (
    "name",
    "ra",
    "dec",
    "error_radius",
    "b0_detect_likelihood",
    "b1_counts",
    "b1_flux",
    "b1_flux_error",
    "b1_exposure",
    "source_extent",
    "p1_counts",
    "p2_counts",
    "p3_counts",
    "p4_counts",
)


def group_matches(rows):
    # The longitude boundary is from the eROSITA-DE DR1 release.
    covered = [row for row in rows if 179.94423568 < float(row["LII"]) < 359.94423568]
    matches = []
    if covered:
        test = "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.008333334))"
        clause = " OR ".join(test.format(**row) for row in covered)
        query = "SELECT " + ",".join(FIELDS) + " FROM erass1main WHERE " + clause
        params = urlencode(
            dict(REQUEST="doQuery", LANG="ADQL", MAXREC="1000", QUERY=query)
        )
        with urlopen(SERVICE + "?" + params, timeout=120) as f:
            response = f.read()
        for info in ElementTree.fromstring(response).findall(".//{*}INFO"):
            if info.get("name") == "QUERY_STATUS" and info.get("value") != "OK":
                raise RuntimeError(info.text or str(info.attrib))
        table = parse_single_table(BytesIO(response)).to_table()
        for rec in table:
            result = {}
            for col in FIELDS:
                val = rec[col]
                if np.ma.is_masked(val):
                    result[col] = None
                else:
                    val = val.item() if isinstance(val, np.generic) else val
                    result[col] = (
                        val.decode().strip() if isinstance(val, bytes) else val
                    )
            matches.append(result)
    out = []
    for row in rows:
        near = []
        is_covered = row in covered
        if is_covered:
            for m in matches:
                d = separation(float(row["RA"]), float(row["DEC"]), m["ra"], m["dec"])
                if d <= RADIUS:
                    near.append(dict(separation_arcsec=round(d, 4), **m))
        out.append(
            dict(
                SRCID=row["SRCID"],
                IAUNAME=row["IAUNAME"],
                in_western_hemisphere=is_covered,
                matches=sorted(near, key=lambda x: x["separation_arcsec"]),
            )
        )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("tables", nargs="+", type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.tables:
        with path.open() as f:
            rows.extend(csv.DictReader(f))
    groups = [rows[i : i + 8] for i in range(0, len(rows), 8)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = [r for group in pool.map(group_matches, groups) for r in group]
    args.output.write_text(
        json.dumps(
            dict(
                table="erass1main",
                service=SERVICE,
                search_radius_arcsec=RADIUS,
                objects=results,
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

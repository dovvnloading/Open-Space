"""Independently enumerate the parent and inspect final TAP status plus all local predicates."""
import csv
import hashlib
import json
import math
import pathlib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

from fetch import GAIA, ROOT


def main():
    query = "SELECT source_id\n" + (ROOT / "count.adql").read_text().split("\n", 1)[1].strip()
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable_plain", MAXREC=10000, QUERY=query)
    url = GAIA + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=180) as response:
        payload = response.read()
        date = response.getheader("Date")
    root = ET.fromstring(payload)
    statuses = [s.get("value") for s in root.findall(".//{*}INFO[@name='QUERY_STATUS']")]
    if not statuses or any(s != "OK" for s in statuses):
        raise RuntimeError(f"VOTable status {statuses}")
    ids = [tr.find("{*}TD").text.strip() for tr in root.findall(".//{*}TR")]
    if not ids:
        raise ValueError("No TABLEDATA rows")
    csv_rows = list(csv.DictReader((ROOT / "parent.csv").open()))
    count = int(next(csv.DictReader((ROOT / "count.csv").open()))["n"])
    if {r["source_id"] for r in csv_rows} != set(ids) or len(set(ids)) != len(ids) or len(ids) != count:
        raise ValueError("Cross-checks on census do not agree")
    for row in csv_rows:
        p = float(row["parallax"])
        mu = float(row["pm"])
        M = float(row["phot_g_mean_mag"]) + 5 * math.log10(p) - 10
        if not (mu > 35 and float(row["bp_rp"]) < 1.2 and p / float(row["parallax_error"]) > 3
                and 4.74047 * mu > 1000 * p and M > 3.5):
            raise ValueError("Source fails local predicate: " + row["source_id"])
    (ROOT / "parent_validation.xml").write_bytes(payload)
    report = dict(endpoint=GAIA, params=params, server_date=date,
                  finished_utc=datetime.now(timezone.utc).isoformat(), sha256=hashlib.sha256(payload).hexdigest(),
                  tap_statuses=statuses, distinct_ids=len(set(ids)), independent_count=count,
                  locally_meet_all_cuts=len(csv_rows))
    (ROOT / "parent_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report)

if __name__ == "__main__":
    main()

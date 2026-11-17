"""Fetch Bailer-Jones et al. EDR3 distance percentiles for the retained Gaia IDs."""

import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from fetch_gaia import ENDPOINT

COLS = (
    "source_id,r_med_geo,r_lo_geo,r_hi_geo,r_med_photogeo,r_lo_photogeo,r_hi_photogeo"
)


def query(block):
    q = (
        "SELECT "
        + COLS
        + " FROM external.gaiaedr3_distance WHERE source_id IN ("
        + ",".join(block)
        + ") ORDER BY source_id"
    )
    raw = urlopen(
        ENDPOINT
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)),
        timeout=120,
    ).read()
    text = raw.decode().splitlines()
    assert text[0].lower().replace('"', "") == COLS
    return raw, q, list(csv.DictReader(text))


if __name__ == "__main__":
    ids = sorted(
        {r["SOURCE_ID"] for r in csv.DictReader(open("desi_screen.csv"))}
        | (
            {r["source_id"] for r in csv.DictReader(open("desi/los/leads.csv"))}
            if Path("desi/los/leads.csv").exists()
            else set()
        ),
        key=int,
    )
    rows = []
    metadata = []
    with ThreadPoolExecutor(4) as pool:
        for raw, q, rs in pool.map(
            query, [ids[i : i + 60] for i in range(0, len(ids), 60)]
        ):
            metadata.append(
                dict(query=q, rows=len(rs), sha256=hashlib.sha256(raw).hexdigest())
            )
            rows += rs
    with open("bj.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS.split(","))
        w.writeheader()
        w.writerows(rows)
    Path("bj_receipt.json").write_text(json.dumps(metadata, indent=2) + "\n")

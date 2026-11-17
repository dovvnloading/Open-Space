"""Freeze ESA Hubble Science Archive image footprints intersecting the five hosts."""

import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
VISITS = HERE.parents[1] / "jwst-nearby-companions/inventory/archive/visits.csv"
HOSTS = {
    "lawd66": "4139531467491239680",
    "hd53143": "5479222240596469632",
    "hd38206": "2967287596507432704",
    "hd21997": "5084963312478488960",
    "hd202628": "6577318075351824512",
}
COLUMNS = "obs_collection,obs_id,obs_publisher_did,target_name,instrument_name,t_min,t_max,t_exptime,s_ra,s_dec,calib_level,dataproduct_type,access_url"


def fetch(item):
    target, source_id = item
    with VISITS.open() as stream:
        visit = next(
            row
            for row in csv.DictReader(stream)
            if row["gcns_edr3_source_id"] == source_id
        )
    query = (
        f"SELECT TOP 3000 {COLUMNS} FROM ivoa.obscore WHERE dataproduct_type='image' "
        f"AND INTERSECTS(s_region, CIRCLE('ICRS',{visit['targ_ra']},{visit['targ_dec']},0.0056))=1"
    )
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    url = "https://hst.esac.esa.int/tap-server/tap/sync?" + urlencode(params)
    with urlopen(url, timeout=90) as response:
        raw = response.read()
    table = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if table.fieldnames != COLUMNS.split(","):
        raise ValueError(f"{target}: {raw[:240]}")
    rows = [dict(target=target, **row) for row in table]
    if len(rows) >= 3000:
        raise ValueError(f"{target}: incomplete page")
    rawdir = HERE / "receipts"
    rawdir.mkdir(exist_ok=True)
    (rawdir / f"{target}.csv.gz").write_bytes(gzip.compress(raw, mtime=0))
    return rows, dict(
        target=target,
        url=url,
        retrieved_utc=datetime.now(timezone.utc).isoformat(),
        rows=len(rows),
        bytes=len(raw),
        uncompressed_sha256=hashlib.sha256(raw).hexdigest(),
    )


if __name__ == "__main__":
    with ThreadPoolExecutor(4) as pool:
        groups = list(pool.map(fetch, HOSTS.items()))
    with (HERE / "hst_footprints.csv").open("w") as stream:
        writer = csv.DictWriter(stream, ["target", *COLUMNS.split(",")])
        writer.writeheader()
        for rows, _ in groups:
            writer.writerows(rows)
    (HERE / "query_receipts.json").write_text(
        json.dumps([r for _, r in groups], indent=2) + "\n"
    )

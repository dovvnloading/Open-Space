"""Inspect the Gaia DR3 neighbourhood around the four highlighted DESI positions."""

from urllib.parse import urlencode
from urllib.request import urlopen
from fetch_gaia import ENDPOINT
from make_queue import rows, save
import csv, hashlib, json
from pathlib import Path

ids = {
    "1717410430709663616",
    "4444095982221534720",
    "2752545725694798336",
    "6340918677332110720",
}
allrows = []
meta = []
for g in rows("gaia_queue.csv"):
    if g["source_id"] not in ids:
        continue
    p = "POINT('ICRS',ra,dec)"
    c = f"POINT('ICRS',{g['ra']},{g['dec']})"
    q = f"SELECT source_id,ra,dec,phot_g_mean_mag,phot_bp_rp_excess_factor,DISTANCE({p},{c})*3600 AS sep_as FROM gaiadr3.gaia_source WHERE 1=CONTAINS({p},CIRCLE('ICRS',{g['ra']},{g['dec']},0.002778)) ORDER BY source_id"
    raw = urlopen(
        ENDPOINT
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)),
        timeout=60,
    ).read()
    rs = list(csv.DictReader(raw.decode().splitlines()))
    assert any(r["source_id"] == g["source_id"] for r in rs)
    allrows.extend(dict(center=g["source_id"], **r) for r in rs)
    meta.append(dict(query=q, sha256=hashlib.sha256(raw).hexdigest()))
save("neighbors.csv", allrows)
Path("neighbors_receipt.json").write_text(json.dumps(meta, indent=2) + "\n")

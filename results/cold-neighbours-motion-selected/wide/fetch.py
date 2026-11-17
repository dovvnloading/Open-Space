"""Independent bright search in 24 disjoint equal-radius spherical caps."""
import csv
import hashlib
import io
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "pilot"))
from fetch_catwise import COLUMNS, ENDPOINT, PREFILTER

def fetch(index):
    # Equal-area placement in latitude; the smallest pair separation exceeds 15 degrees.
    ra = index * (180 * (3 - math.sqrt(5))) % 360
    dec = math.degrees(math.asin(1 - 2 * (index + 0.5) / 24))
    circle = f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.7f},{dec:.7f},7.5))"
    limits = PREFILTER + " AND w2mpro_pm <= 15.5 AND w2rchi2_pm <= 3 AND w2ab_map = 0"
    query = f"SELECT {COLUMNS} FROM catwise_2020 WHERE {circle} AND {limits} ORDER BY source_id"
    uri = ENDPOINT + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=50000, QUERY=query))
    path = ROOT / "inputs" / f"cap{index:02d}.csv"
    if not path.exists():
        with urlopen(uri, timeout=300) as response:
            payload = response.read()
        if not payload.startswith(b"source_name,"):
            raise ValueError(payload[-500:])
        path.write_bytes(payload)
    data = path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if len(rows) >= 50000 or len({r["source_id"] for r in rows}) != len(rows):
        raise ValueError("Truncated or duplicate source IDs")
    print(index,len(rows),flush=True)
    return dict(cap=index, ra=ra, dec=dec, radius_deg=7.5, query=query, rows=len(rows), sha256=hashlib.sha256(data).hexdigest())

if __name__ == "__main__":
    (ROOT / "inputs").mkdir(exist_ok=True)
    with ThreadPoolExecutor(3) as pool:
        manifest=list(pool.map(fetch,range(24)))
    (ROOT/"inputs_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")

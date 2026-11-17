"""Recreate the bounded discovery sample from IRSA's public synchronous TAP endpoint."""
import argparse
import csv
import hashlib
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from time import sleep

ROOT = Path(__file__).resolve().parent
FIELDS = {"north": (269.733, 66.018), "south": (61.241, -48.423), "fornax": (52.932, -28.088)}
ENDPOINT = "https://irsa.ipac.caltech.edu/TAP/sync"
COLUMNS = (
    "source_name,source_id,ra,dec,ra_pm,dec_pm,meanobsmjd,pmra,pmdec,sigpmra,sigpmdec,"
    "w1snr_pm,w2snr_pm,w1mpro_pm,w2mpro_pm,w1sigmpro_pm,w2sigmpro_pm,"
    "w1flux_pm,w1sigflux_pm,w2flux_pm,w2sigflux_pm,w1rchi2_pm,w2rchi2_pm,"
    "w1cov,w2cov,pmcode,w2ab_map,cc_flags,nb,na,par_pm,par_pmsig"
)
PREFILTER = (
    "w2snr_pm >= 8 AND pmra*pmra + pmdec*pmdec >= 0.04 "
    "AND (w1snr_pm < 3 OR w1snr_pm IS NULL OR w1mpro_pm - w2mpro_pm >= 1.5)"
)

def queries():
    for field, (ra, dec) in FIELDS.items():
        footprint = f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},3.0))"
        for suffix, columns, extra in [
            ("rows", COLUMNS, f" AND {PREFILTER} ORDER BY source_id"),
            ("filtered_count", "COUNT(*) AS n", f" AND {PREFILTER}"),
            ("all_count", "COUNT(*) AS n", ""),
            ("snr_count", "COUNT(*) AS n", " AND w2snr_pm >= 8"),
        ]:
            yield f"{field}_{suffix}", f"SELECT {columns} FROM catwise_2020 WHERE {footprint}{extra}"

def download(job):
    label, query = job
    data_path = ROOT / "inputs" / f"{label}.csv"
    (ROOT / "queries" / f"{label}.adql").write_text(query+"\n")
    url = ENDPOINT + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=50000, QUERY=query))
    if not data_path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=180) as response:
                    data = response.read()
                if not data.startswith((b"n\n", b"n\r\n", b"source_name,")):
                    raise ValueError(f"non-CSV response for {label}: {data[:200]!r}")
                data_path.write_bytes(data)
                break
            except Exception:
                if attempt == 2: raise
                sleep(3 * (attempt+1))
    data = data_path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    print(label, len(rows) if label.endswith("_rows") else rows, flush=True)
    return dict(file=f"inputs/{data_path.name}", query=f"queries/{label}.adql", bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), data_rows=len(rows))

def main():
    for name in ["inputs","queries"]: (ROOT/name).mkdir(exist_ok=True)
    with ThreadPoolExecutor(4) as pool:
        manifest = list(pool.map(download,queries()))
    for field in FIELDS:
        p = ROOT/"inputs"
        with (p/f"{field}_rows.csv").open() as f: rows = list(csv.DictReader(f))
        with (p/f"{field}_filtered_count.csv").open() as f: count = int(next(csv.DictReader(f))["n"])
        if count != len(rows) or len({r["source_id"] for r in rows}) != count:
            raise ValueError(f"truncated or duplicated response in {field}")
    (ROOT/"catwise_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")

if __name__ == "__main__":
    main()

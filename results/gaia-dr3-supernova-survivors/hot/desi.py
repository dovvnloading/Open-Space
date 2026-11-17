"""Exact-ID crossmatch of the Culpan catalogue with DESI DR1 MWS on Data Lab."""
import csv
import hashlib
import io
import json
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from screen import catalogue, HERE

URL = "https://datalab.noirlab.edu/tap/sync"
COLS = "source_id,targetid,survey,program,healpix,rv_adop,rv_err,vrad,vrad_err,rvs_warn,zcat_primary"
BATCH = 300


def request(ids):
    q = f"SELECT {COLS} FROM desi_dr1.mws WHERE source_id IN (" + ",".join(ids) + ")"
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=10000, QUERY=q)
    raw = urllib.parse.urlencode(params).encode()
    with urllib.request.urlopen(URL, data=raw, timeout=180) as response:
        payload = response.read()
        if payload.startswith(b"<") or not payload.startswith(b"source_id,"):
            raise ValueError(payload[:300])
        rows = list(csv.DictReader(io.StringIO(payload.decode())))
        assert len(rows) < params["MAXREC"]
        assert all(r["source_id"] in ids for r in rows)
        receipt = dict(method="POST", url=URL, params=params, nbytes=len(payload),rows=len(rows),
                       date=response.getheader("Date"), sha256=hashlib.sha256(payload).hexdigest())
    return rows, receipt


def download():
    ids = sorted({r["GaiaEDR3"] for r in catalogue() if r["GaiaEDR3"]})
    with ThreadPoolExecutor(max_workers=3) as pool:
        chunks = list(pool.map(request, [ids[i:i+BATCH] for i in range(0,len(ids),BATCH)]))
    cat = [r for rows,_ in chunks for r in rows]
    with (HERE/"desi_rvs.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=COLS.split(","),lineterminator="\n");w.writeheader();w.writerows(cat)
    receipts = dict(finished_utc=datetime.now(timezone.utc).isoformat(), n_input=len(ids), chunks=[r for _,r in chunks])
    (HERE/"desi_receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")


def analyse():
    with (HERE/"desi_rvs.csv").open() as f:
        rows=list(csv.DictReader(f))
    good = [r for r in rows if r["rv_adop"] and r["rv_err"] and float(r["rv_err"])>0]
    clean = [r for r in good if int(r["rvs_warn"])==0]
    tails = {r["source_id"] for r in good if abs(float(r["rv_adop"]))>=700}
    summary=dict(input_distinct_ids=len({r["GaiaEDR3"] for r in catalogue() if r["GaiaEDR3"]}),
                 match_rows=len(rows),matched_ids=len({r["source_id"] for r in rows}),
                 usable_rv_ids=len({r["source_id"] for r in good}),
                 clean_rv_ids=len({r["source_id"] for r in clean}),
                 max_abs_clean_rv=max(abs(float(r["rv_adop"])) for r in clean),
                 raw_ge_700_ids=sorted(tails),clean_ge_700_ids=sorted({r["source_id"] for r in clean if abs(float(r["rv_adop"]))>=700}),
                 all_visits_of_raw_tail=[r for r in rows if r["source_id"] in tails])
    (HERE/"desi_summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__ == "__main__":
    if "--download" in sys.argv: download()
    analyse()

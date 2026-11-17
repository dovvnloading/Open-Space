"""Freeze ALeRCE TAP CSV replies following results/rubin-alerce-2026/rapid/retrieve.py."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
TAP = "https://tap.alerce.online/tap/sync"
CAP = 100000
OFFSET = 37 / 86400
FIELDS = """o.oid,o.meanra,o.meandec,o.n_det,o.firstmjd,d.measurement_id,d.mjd,d.band,p.visit,p.detector,
d.ra,d.dec,p.psfflux,p.psffluxerr,p.reliability,p.ssobjectid,p.diaobjectid,
p.isnegative,p.isdipole,p.glint_trail,p.psfflux_flag,p.pixelflags_edge,p.pixelflags_crcenter,
p.pixelflags_saturatedcenter,p.pixelflags_streakcenter,p.pixelflags_injectedcenter,
p.pixelflags_injected_templatecenter"""
def cohort(day):
    return (f"o.sid=1 AND o.firstmjd>={day + OFFSET:.12f} "
            f"AND o.firstmjd<{day+1 + OFFSET:.12f} AND o.n_det>=2")
def queries():
    return {f"day_{day}": f"""SELECT {FIELDS}
FROM alerce_tap.object AS o
JOIN alerce_tap.detection AS d ON d.oid=o.oid AND d.sid=o.sid
JOIN alerce_tap.lsst_detection AS p ON p.measurement_id=d.measurement_id AND p.oid=d.oid AND p.sid=d.sid
WHERE {cohort(day)} AND d.mjd>={61345+OFFSET:.12f} AND d.mjd<{61361+OFFSET:.12f}"""
            for day in range(61345,61361)}
def fetch(item, target):
    name, sql=item
    params=dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",MAXREC=CAP,QUERY=sql)
    url=TAP+"?"+urlencode(params)
    receipt=dict(key=name,adql=sql,url=url,utc_start=datetime.now(timezone.utc).isoformat())
    try:
        with urlopen(url,timeout=240) as response:
            body=response.read()
            receipt.update(status=response.status, http_date=response.headers.get("Date"))
        rows=csv.DictReader(io.StringIO(body.decode()))
        if "oid" not in (rows.fieldnames or []):
            raise ValueError("TAP CSV schema mismatch: "+body[:300].decode(errors="replace"))
        n=sum(1 for _ in rows)
        receipt.update(rows=n,bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
        if n>=CAP: raise ValueError("MAXREC ceiling reached; partition is incomplete")
        (target/(name+".csv.gz")).write_bytes(gzip.compress(body,compresslevel=9,mtime=0))
    except Exception as e:
        receipt["error"]=repr(e)
        if isinstance(e,HTTPError):
            receipt.update(status=e.code,error_body=e.read(500).decode(errors="replace"))
    receipt["utc_end"]=datetime.now(timezone.utc).isoformat()
    (target/(name+".json")).write_text(json.dumps(receipt,indent=2)+"\n")
    print(name,receipt.get("rows"),receipt.get("error"),flush=True)
    return receipt
if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    args=parser.parse_args()
    if args.directory.resolve().is_relative_to(ROOT.parent.resolve()):
        raise ValueError('New acquisition must be outside the frozen result')
    args.directory.mkdir(parents=True,exist_ok=False)
    with ThreadPoolExecutor(3) as executor:
        receipts=list(executor.map(lambda item:fetch(item,args.directory),queries().items()))
    if any("error" in item for item in receipts):
        raise SystemExit("At least one Rubin partition is incomplete")

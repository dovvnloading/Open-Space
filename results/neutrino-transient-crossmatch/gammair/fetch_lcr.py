"""Fetch the NASA LCR source list and all eligible source responses; resume safely."""
import argparse
import concurrent.futures
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

BASE = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/queryDB.php"
LIST_URL = BASE + "?typeOfRequest=SourceList&catalog=4FGL"
HERE = Path(__file__).resolve().parent

def eligible(s):
    return s["CLASS1"].strip().lower() in {"bll","fsrq","bcu"} and float(s["Variability_Index"]) > 21.67

def light_curve_url(name):
    return BASE + "?typeOfRequest=lightCurveData&source_name=" + urllib.parse.quote(name, safe="") + "&cadence=weekly&flux_type=photon&index_type=fixed&ts_min=4"

def fetch(path, url, kind):
    receipt = path.with_suffix(".receipt")
    if path.exists() and receipt.exists():
        return
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=40) as response:
                raw = response.read()
                stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
                obj = json.loads(raw)
                if kind == "curve" and not {"ts","flux","flux_error","fit_convergence"}.issubset(obj):
                    raise ValueError("unexpected light curve schema")
                if kind == "list" and not (isinstance(obj,list) and len(obj)>1000):
                    raise ValueError("unexpected source list")
                path.with_suffix(".part").write_bytes(gzip.compress(raw, compresslevel=5, mtime=0))
                path.with_suffix(".part").replace(path)
                record = dict(url=url, retrieved_utc=stamp, http_status=response.status,
                              http_date=response.headers.get("Date"), bytes=len(raw),
                              sha256=hashlib.sha256(raw).hexdigest())
                receipt.write_text(json.dumps(record)+"\n")
                return
        except (OSError, ValueError) as exc:
            if attempt == 2:
                raise
            time.sleep(2**attempt)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path)
    p.add_argument("--workers",type=int,default=6)
    args=p.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    manifest=args.cache / "4FGL.json.gz"
    fetch(manifest, LIST_URL, "list")
    sources=[s for s in json.loads(gzip.decompress(manifest.read_bytes())) if eligible(s)]
    def one(s):
        name=s["Source_Name"]
        path=args.cache / (name.replace(" ","_")+".json.gz")
        try:
            fetch(path, light_curve_url(name),"curve")
        except Exception as exc:
            return f"{name}: {exc}"
    failed=0
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for idx,err in enumerate(pool.map(one,sources),1):
            if err:
                failed+=1
                print(err,flush=True)
            if idx%100==0:
                print(f"Processed {idx}/{len(sources)}",flush=True)
    if failed:
        raise SystemExit(f"{failed} retrieval failures; rerun to resume")

if __name__=="__main__":
    main()

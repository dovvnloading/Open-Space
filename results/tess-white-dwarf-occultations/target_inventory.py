"""Inventory every released MAST bulk target-pixel roster; never execute a roster."""
import argparse
import csv
import hashlib
import json
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

INDEX = "https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html"
FILE = re.compile(r"\btess\d+-s\d+-(\d+)-\d+-[a-z]+_(?:fast-)?tp\.fits\b")
SCRIPT = re.compile(r'href="([^"]*/tesscurl_sector_(\d+)_(fast-tp|tp)\.sh)"')

def get(url, dest):
    if not dest.exists():
        with urllib.request.urlopen(url, timeout=100) as response:
            data = response.read()
        tmp = dest.with_suffix(".tmp")
        tmp.write_bytes(data)
        tmp.replace(dest)
    data = dest.read_bytes()
    return data, dict(url=url, bytes=len(data), sha256=hashlib.sha256(data).hexdigest())

def inventory(out, cache, targets):
    out.mkdir(exist_ok=True, parents=True)
    cache.mkdir(exist_ok=True, parents=True)
    index, receipt = get(INDEX, cache / "bulk-index.html")
    index = index.decode()
    rows = SCRIPT.findall(index)
    ffi = sorted({int(s) for s in re.findall(r'tesscurl_sector_(\d+)_ffic\.sh',index)})
    def parse(entry):
        link, sector, mode = entry
        url = urljoin(INDEX, link)
        cadence = 20 if mode == "fast-tp" else 120
        try:
            data, meta = get(url, cache / link.split("/")[-1])
            names = re.findall(r'-o\s+(tess\S+\.fits)\s+', data.decode())
            matches = [FILE.fullmatch(n) for n in names]
            if not names or any(m is None for m in matches):
                raise ValueError("unexpected filename schema")
            selected = [(int(m[1]), int(sector), cadence, n) for n,m in zip(names,matches) if int(m[1]) in targets]
            return {**meta,"sector":int(sector),"cadence_s":cadence,"products":len(names)}, selected
        except Exception as exc:
            return dict(url=url,sector=int(sector),cadence_s=cadence,error=str(exc)), []
    with ThreadPoolExecutor(8) as pool:
        results = list(pool.map(parse, rows))
    catalog = sorted((r[0] for r in results), key=lambda r:(r["sector"],r["cadence_s"]))
    with (out / "target_pixel_products.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["tic","sector","cadence_s","filename"])
        w.writerows(sorted(row for _,products in results for row in products))
    (out / "bulk_receipts.json").write_text(json.dumps(dict(index=receipt,read_utc=datetime.now(timezone.utc).isoformat(),ffi_sectors=ffi,scripts=catalog),indent=2)+"\n")
    print("scripts",len(catalog),"errors",sum("error" in x for x in catalog),"products",sum(x.get("products",0) for x in catalog))

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets",type=Path,help="CSV with a tic column")
    p.add_argument("--out",type=Path,default=Path("record"))
    p.add_argument("--cache",type=Path,default=Path("cache"))
    a=p.parse_args()
    targets=set()
    if a.targets:
        with a.targets.open() as f:
            targets={int(r["tic"]) for r in csv.DictReader(f) if r["tic"]}
    inventory(a.out,a.cache,targets)

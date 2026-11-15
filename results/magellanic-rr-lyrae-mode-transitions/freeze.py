"""Freeze the eight responses used to reproduce the target measurement offline."""
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

root = Path(__file__).parent
m = json.loads((root/"measurements.json").read_text())
c = json.loads((root/"independent-checks.json").read_text())
urls = [*m["ogle_sha256"], *m["tap_queries"], c["eros"]["url"]]


def retrieve(url):
    for attempt in range(3):
        try:
            with urlopen(url,timeout=45) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise


with ThreadPoolExecutor(8) as pool:
    inputs = list(pool.map(retrieve,urls))
manifest = []
with ZipFile(root/"target-inputs.zip","w") as dest:
    def add(name, data):
        entry = ZipInfo(name, (2026,11,15,0,0,0))
        dest.writestr(entry,data,compress_type=ZIP_DEFLATED,compresslevel=9)
    for url, raw in zip(urls,inputs):
        key = hashlib.sha256(url.encode()).hexdigest()
        digest = hashlib.sha256(raw).hexdigest()
        if m["ogle_sha256"].get(url,digest) != digest:
            raise ValueError(f"OGLE digest changed for {url}")
        add(key,raw)
        manifest.append(dict(member=key,url=url,bytes=len(raw),sha256=digest))
    add("manifest.json",json.dumps(manifest,indent=2)+"\n")

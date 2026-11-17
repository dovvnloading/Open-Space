"""Measure telescope-to-barycentric factors from DESI exposure primary FITS headers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import subprocess
from hashlib import sha256
from pathlib import Path
import zlib

import fitsio
from astropy.io import fits

BASE = "https://desidata.s3.amazonaws.com/dr1/spectro/redux/iron/exposures"


def get_header(key):
    night, expid, camera = key
    url = f"{BASE}/{night}/{expid:08d}/cframe-{camera}-{expid:08d}.fits.gz"
    payload = subprocess.check_output(["curl","-fsSL","--retry","4","--retry-all-errors",
        "--max-time","25","--range","0-32767","--max-filesize","32768",url])
    if len(payload) != 32768:
        raise ValueError(f"short header prefix: {url}")
    raw = zlib.decompressobj(wbits=31).decompress(payload)
    cards = [raw[i:i+80] for i in range(0,len(raw)-79,80)]
    end = next(i for i,c in enumerate(cards) if c[:8].strip() == b"END")
    hdr = fits.Header.fromstring(b"".join(cards[:end+1]).decode("ascii"), sep="")
    if (hdr["NIGHT"], hdr["EXPID"], hdr["CAMERA"].strip()) != key:
        raise ValueError(f"mismatched frame header {url}")
    return {"night":night, "expid":expid, "camera":camera,
            "heliocor":hdr["HELIOCOR"], "sha256_prefix32768":sha256(payload).hexdigest(),
            "url":url}


def retrieve(tables, output):
    exposures=set()
    for tab in tables:
        fib=fitsio.read(str(tab),ext="FIBERMAP")
        exposures.update((int(x["NIGHT"]),int(x["EXPID"]),a+str(x["PETAL_LOC"]))
                         for x in fib for a in "brz")
    work = output.with_suffix(".parts")
    work.mkdir(exist_ok=True)
    def cached(key):
        path = work / f"{key[0]}-{key[1]}-{key[2]}.json"
        if path.exists():
            return json.loads(path.read_text())
        row = get_header(key)
        path.write_text(json.dumps(row))
        return row
    with ThreadPoolExecutor(max_workers=5) as ex:
        rows=list(ex.map(cached, sorted(exposures)))
    with output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("rvtab",type=Path,nargs="+")
    a=parser.parse_args()
    retrieve(a.rvtab,a.output)

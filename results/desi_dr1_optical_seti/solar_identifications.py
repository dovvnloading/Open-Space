"""Context-only solar line identities and displaced-centre controls for DESI candidates."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from vetting import air_to_vacuum

ATLAS = "https://raw.githubusercontent.com/AlexPietrow/HelioSpectrotron5000/f32e5417383128593c76ee4f9552f1e0e5621c02/"
PRODUCTS = (
    ("moore_clean02012026.csv","b84107af4d1d6c06fba1e0a16a34f419c39002420eb0d9e14634b912bca65edb","wavelength","ew"),
    ("babcock_clean_02012026.csv","d2071760dc9dfb609e4923d9dbd019ea16ac26bfe52918fd30d49ad07cb076a4","wav","strength"),
)
OFFSETS = (-35,-27,-19,-11,11,19,27,35)


def read_lines(cache):
    lines=[]
    cache.mkdir(parents=True,exist_ok=True)
    for file,h,col,strength in PRODUCTS:
        path=cache/file
        raw=path.read_bytes() if path.exists() else urlopen(ATLAS+file,timeout=45).read()
        if hashlib.sha256(raw).hexdigest()!=h:
            raise ValueError("input atlas digest differs")
        if not path.exists():
            path.write_bytes(raw)
        for r in csv.DictReader(io.StringIO(raw.decode())):
            air=float(r[col])
            species=r["id"].strip()
            if not species:
                continue
            lines.append((air_to_vacuum(air),species,r[strength],r["flag"],file))
    lines.sort()
    return lines,np.array([r[0] for r in lines])


def report(catalog, cache, output):
    lines, waves=read_lines(cache)
    def near(x, half):
        lo,hi=np.searchsorted(waves,[x-half,x+half])
        return lines[lo:hi]
    with gzip.open(catalog,"rt") as source, gzip.open(output,"wt") as dest:
        for s in source:
            r=json.loads(s)
            center=r["rest_A"]
            window=min(v["fwhm_rest_A"] for v in r["visits"])/2
            match=near(center,window)
            support=[(lo,hi) for lo,hi in ((waves[0],waves[-1]),) if lo<=center<=hi]
            controls=[int(bool(near(center+offset,window))) for offset in OFFSETS
                      if any(lo<=center+offset<=hi for lo,hi in support)]
            out={k:r[k] for k in ("targetid","healpix","arm","rest_A")}
            out.update(half_window_A=window,
                # Each source book uses its own strength scale. Preserve the original field.
                matches=[dict(vacuum_A=round(w,4),id=name,source_strength=strength,source_flag=flag,
                              source="Moore1966" if src.startswith("moore") else "Babcock1947") for w,name,strength,flag,src in match],
                shifted_positions_with_match=sum(controls),shifted_positions=len(controls))
            print(json.dumps(out,separators=(",",":")),file=dest)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("catalog",type=Path)
    p.add_argument("cache",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    report(a.catalog,a.cache,a.output)

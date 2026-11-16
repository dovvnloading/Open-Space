"""Verify original Fink cone responses, recreate the sample and validate image digests."""
import csv
import gzip
import hashlib
import json
import math
from urllib.request import urlopen
from cutouts import DATA, rows, read_hdus

STARS=DATA/"stellar"
def select(t):
    found=json.loads((STARS/f"fink_star{t['measurement_id']}.json").read_bytes())
    candidates=[]
    for r in found:
        if str(r["r:visit"])!=t["visit"] or str(r["r:detector"])!=t["detector"]:
            continue
        s=r["r:scienceFlux"] or 0
        template=r["r:templateFlux"] or 0
        sn=s/(r["r:scienceFluxErr"] or 1)
        if sn<35 or s>300000 or not .5*s<=template<=1.5*s or str(r["r:diaSourceId"])==t["measurement_id"]:
            continue
        dra=(r["r:ra"]-float(t["ra"]))*math.cos(math.radians(float(t["dec"])))
        separation=3600*math.hypot(dra,r["r:dec"]-float(t["dec"]))
        candidates.append((sn,str(r["r:diaSourceId"]),separation))
    return sorted(candidates,key=lambda s:-s[0])[:12]

def main():
    with (STARS/"candidates.csv").open() as f:
        candidates=list(csv.DictReader(f))
    selected=[mid for t in rows("target.csv") for sn,mid,dist in select(t)]
    if selected != [r["mid"] for r in candidates]:
        raise ValueError("Frozen point-source selection not reproduced")
    with (STARS/"stamps.csv").open() as f:
        manifest=list(csv.DictReader(f))
    if set(selected)!={r["id"] for r in manifest}:
        raise ValueError("Incomplete comparison-star images")
    for r in manifest:
        b=gzip.decompress((STARS/r["file"]).read_bytes())
        if hashlib.sha256(b).hexdigest()!=r["sha256_fits"]:
            raise ValueError(f"Unexpected image bytes in {r['file']}")
        read_hdus(b)
    print(f"Verified {len(selected)} comparison-star images and the selection from the original queries.")
if __name__=="__main__":
    main()

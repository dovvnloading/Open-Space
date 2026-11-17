"""Query CADC's image index with explicit ephemeris samples for C/2025 V1."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode

out=Path(__file__).parent/"image-search"
response=json.loads(gzip.decompress((out/"horizons.json.gz").read_bytes()))
before,body=response["result"].split("$$SOE",1)
header=next(csv.reader([next(line for line in reversed(before.splitlines()) if "R.A." in line)]))
ri=next(i for i,c in enumerate(header) if "R.A." in c)
di=next(i for i,c in enumerate(header) if "DEC_" in c)
# MJD, J2000 RA in degrees, J2000 Dec in degrees. All telescope sites use
# the Palomar ephemeris here. A one-arcminute padding accommodates parallax.
ephem="\n".join(f"{float(r[0])-2400000.5:.6f} {float(r[ri]):.8f} {float(r[di]):.8f}" for r in csv.reader(body.split("$$EOE")[0].strip().splitlines()))
data={"search":"userephem","ephem":ephem,"format":"tsv","lang":"en","eelipse":"1","eunits":"arcminutes","xyres":"no"}
uri="https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"
post=urlencode(data).encode()
record={"url":uri,"method":"POST","parameters":{k:v for k,v in data.items() if k!="ephem"},"n_ephemeris_rows":len(ephem.splitlines()),"query_bytes":len(post),"query_sha256":hashlib.sha256(post).hexdigest()}
try:
    req=Request(uri,data=post,headers={"User-Agent":"Open-Space-footprints/1.0 (+https://github.com/dovvnloading/Open-Space)"})
    with urlopen(req,timeout=100) as resp:
        raw=resp.read()
        record.update(status=resp.status,date=resp.headers.get("Date"),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        (out/"cadc-retry-coordinate.tsv").write_bytes(raw)
        print("CADC",raw[:2200],flush=True)
except Exception as e:
    record["error"]=str(e)
(out/"cadc-retry-coordinate-receipt.json").write_text(json.dumps(record,indent=2)+"\n")

# The alternative local orbital-element mode keeps ephemeris generation at
# the host archive. The supplied orientation is heliocentric ecliptic J2000.
import math
jpl=json.loads(gzip.decompress((out.parent/"sbdb.json.gz").read_bytes()))
r=next(dict(zip(jpl["fields"],v)) for v in jpl["data"] if str(v[0])=="1004097")
a=float(r["q"])/(1-float(r["e"]))
m=0.9856076686*(float(r["epoch"])-float(r["tp"]))/abs(a)**1.5
params={"search":"orbel","epoch":float(r["epoch"])-2400000.5,"a":a,"e":r["e"],"i":r["i"],
        "longi":r["om"],"peri":r["w"],"anom":m,"epoch1":"2025-07-01","epoch2":"2025-11-02",
        "format":"tsv","lang":"en","xyres":"no","eelipse":1,"eunits":"arcminutes"}
url=uri+"?"+urlencode(params)
receipt={"url":url,"method":"GET"}
try:
    with urlopen(Request(url,headers={"User-Agent":"Open-Space-footprints/1.0 (+https://github.com/dovvnloading/Open-Space)"}),timeout=75) as resp:
        raw=resp.read()
        receipt.update(status=resp.status,date=resp.headers.get("Date"),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        (out/"cadc-retry-elements.tsv").write_bytes(raw)
        print("CADC orbital elements",raw[:1200],flush=True)
except Exception as e:
    receipt["error"]=str(e)
(out/"cadc-retry-elements-receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")

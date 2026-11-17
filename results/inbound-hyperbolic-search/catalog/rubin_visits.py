"""Cross-match nominal and three-sigma C/2025 V1 positions with the public Rubin DP2 visit database."""
from bisect import bisect_right
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import tempfile
from urllib.parse import urlencode
from urllib.request import Request, urlopen

OUT=Path("image-search")
DB="https://s3df.slac.stanford.edu/data/rubin/sim-data/prelsst/dp2_visits.db"
PUBLISHED="3ff666cbfe4052db82b864cafdd85332fb55b47e7d7f80e5276866811e3afba6"
doc="https://survey-strategy.lsst.io/progress/prelsst/prelsst_dp2.html"
receipts=[]
def get(url):
    req=Request(url,headers={"User-Agent":"Open-Space-DP2/1.0 (+https://github.com/dovvnloading/Open-Space)"})
    with urlopen(req,timeout=90) as response:
        raw=response.read()
        receipts.append({"url":url,"bytes":len(raw),"status":response.status,"http_date":response.headers.get("Date"),
                         "sha256":hashlib.sha256(raw).hexdigest()})
        return raw

database=get(DB)
if receipts[-1]["sha256"]!=PUBLISHED:
    raise ValueError("DP2 visit database differs from the official SHA-256")
params={"format":"json","COMMAND":"'DES=C/2025 V1;'","EPHEM_TYPE":"'OBSERVER'","CENTER":"'X05'",
        "START_TIME":"'2025-04-24'","STOP_TIME":"'2025-11-09'","STEP_SIZE":"'6 h'","QUANTITIES":"'1,37'",
        "ANG_FORMAT":"'DEG'","CSV_FORMAT":"'YES'","CAL_FORMAT":"'JD'","EXTRA_PREC":"'YES'"}
raw=get("https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(params))
(OUT/"rubin-horizons.json.gz").write_bytes(gzip.compress(raw,mtime=0))
hdr,body=json.loads(raw)["result"].split("$$SOE",1)
header=next(csv.reader([next(l for l in reversed(hdr.splitlines()) if "R.A." in l)]))
i0=next(i for i,h in enumerate(header) if "R.A." in h)
i1=next(i for i,h in enumerate(header) if "DEC" in h)
i2=next(i for i,h in enumerate(header) if "SMAA" in h)
points=[(float(r[0]),float(r[i0]),float(r[i1]),float(r[i2])) for r in csv.reader(body.split("$$EOE")[0].strip().splitlines())]
times=[p[0] for p in points]
def predict(t):
    i=bisect_right(times,t)-1
    t0,ra0,de0,ma0=points[i]
    t1,ra1,de1,ma1=points[i+1]
    frac=(t-t0)/(t1-t0)
    return (ra0+frac*((ra1-ra0+180)%360-180))%360,de0+frac*(de1-de0),max(ma0,ma1)

data=sqlite3.connect(":memory:")
data.deserialize(database)
data.row_factory=sqlite3.Row
tables=[r[0] for r in data.execute("select name from sqlite_master where type='table'")]
rows=data.execute("select * from observations").fetchall()
columns=list(rows[0].keys())
ids=[f for f in columns if f.lower() in {"observationid","observation_id","visit_id","visitid","visit"}]
print("DP2",tables,len(rows),columns,flush=True)
kept,closest=0,None
matches=[]
for r in rows:
    # Opsim schema timestamps use TAI. The 0.02 degree padding also exceeds
    # the along-track error from a UTC/TAI convention difference in this interval.
    exp=float(r["visitExposureTime"] or 0)
    jd=float(r["observationStartMJD"])+2400000.5+(exp/2-37)/86400
    if not times[0]<=jd<times[-1]:
        continue
    kept+=1
    ra,dec,sma=predict(jd)
    a,b,h=map(math.radians,(ra,dec,r["fieldDec"]))
    d=math.radians(r["fieldRA"])-a
    sep=math.degrees(math.acos(max(-1,min(1,math.sin(b)*math.sin(h)+math.cos(b)*math.cos(h)*math.cos(d)))))
    row={"visit":r[ids[0]] if ids else None,"mid_jd_utc":jd,"ra_pred":ra,"dec_pred":dec,
         "ra_center":r["fieldRA"],"dec_center":r["fieldDec"],"smaa_arcsec_conservative":sma,"separation_deg":sep}
    if closest is None or sep<closest["separation_deg"]:
        closest=row
    # 2.1 degrees covers the roughly 1.75 degree LSSTCam field radius including corners.
    if sep<=2.1+sma/3600+.02:
        matches.append(row)
result={"documentation":doc,"published_sha256":PUBLISHED,"total_visits":len(rows),
        "interval_visits":kept,"interval_jd_utc":[times[0],times[-1]],"observatory":"X05",
        "closest":closest,"coverage_radius_deg":2.1,"interpolation_padding_deg":.02,
        "possible_visits":matches,"receipts":receipts}
(OUT/"rubin-visits.json").write_text(json.dumps(result,indent=2)+"\n")
print("RUBIN",len(rows),kept,closest,"possible",len(matches),flush=True)

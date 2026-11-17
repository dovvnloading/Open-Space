"""Search DECam raw exposure centers over the pre-discovery interval, using one archive query per month."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from bisect import bisect_right
from pathlib import Path
import sys
from urllib.request import urlopen, Request

out = Path(sys.argv[1])
horizons = json.loads((out / "horizons.json").read_text())["result"]
before, body = horizons.split("$$SOE", 1)
header = next(csv.reader([next(line for line in reversed(before.splitlines()) if "R.A." in line)]))
ri = next(i for i,col in enumerate(header) if "R.A." in col)
di = next(i for i,col in enumerate(header) if "DEC" in col.upper())
ephem = [(float(r[0]), float(r[ri]), float(r[di])) for r in csv.reader(body.split("$$EOE")[0].strip().splitlines())]
dates = [p[0] for p in ephem]
uri = "https://astroarchive.noirlab.edu/api/adv_search/find/?rectype=file&limit=100000"
receipts, near = [], []
count, missing, closest = 0, 0, None
intervals = ["2025-06-30","2025-08-01","2025-09-01","2025-10-01","2025-11-02"]
for start, stop in zip(intervals, intervals[1:]):
    fields = ["md5sum","archive_filename","dateobs_center","ra_center","dec_center","exposure","ifilter","caldat"]
    q = {"outfields":fields,"search":[["caldat",start,stop],["instrument","decam"],["proc_type","raw"],["obs_type","object"]]}
    req = Request(uri, data=json.dumps(q).encode(), headers={"Content-Type":"application/json","User-Agent":"Open-Space/1.0 (+https://github.com/dovvnloading/Open-Space)"})
    rec = {"url":uri,"method":"POST","query":q}
    receipts.append(rec)
    try:
        with urlopen(req, timeout=120) as resp:
            b = resp.read()
            rec.update(status=resp.status,http_date=resp.headers.get("Date"),bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
        rows = json.loads(b)
        if not isinstance(rows,list) or not rows or "HEADER" not in rows[0] or len(rows) > 100000:
            raise ValueError("Malformed or truncated NOIRLab response")
        rec["rows"] = len(rows)-1
        for row in rows[1:]:
            r={k.split(":")[-1]:v for k,v in row.items()}
            if not r.get("dateobs_center") or r.get("ra_center") is None or r.get("dec_center") is None:
                missing += 1
                continue
            t=datetime.fromisoformat(str(r["dateobs_center"]).strip().replace("Z","+00:00"))
            if t.tzinfo is None:
                t=t.replace(tzinfo=timezone.utc)
            jd=t.timestamp()/86400+2440587.5
            # caldat is a local observing date, often the previous UTC day.
            if not (dates[0] <= jd < dates[-1]) or not (start<=str(r["caldat"])[:10]<stop):
                continue
            i=bisect_right(dates,jd)-1
            t0,ra0,dec0=ephem[i]
            t1,ra1,dec1=ephem[i+1]
            f=(jd-t0)/(t1-t0)
            ra=(ra0+f*((ra1-ra0+180)%360-180))%360
            dec=dec0+f*(dec1-dec0)
            a,b,c=map(math.radians,[float(r["ra_center"])-ra,float(r["dec_center"]),dec])
            distance=math.degrees(math.acos(max(-1,min(1,math.sin(b)*math.sin(c)+math.cos(b)*math.cos(c)*math.cos(a)))))
            count+=1
            closest=distance if closest is None else min(closest,distance)
            # Padded well beyond the approximately 1.1 degree focal-plane radius.
            if distance<1.5:
                near.append({**r,"separation_deg":distance,"ra_pred":ra,"dec_pred":dec})
    except Exception as exc:
        rec["error"]=str(exc)
    print("DECam",start,rec.get("rows"), rec.get("error"),flush=True)
results={"source":"NOIRLab raw DECam object-exposure centers. Palomar topocentric comet coordinates used with a 1.5-degree padded radial prefilter.",
         "rows_compared":count,"missing_coordinates_or_date":missing,"closest_center_deg":closest,"nearby_exposures":near,
         "receipts":receipts}
(out/"decam.json").write_text(json.dumps(results,indent=2)+"\n")
print("DECam compared",count,"near",len(near),"closest",closest)

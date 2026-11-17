"""Independent replay of the historical ASAS-SN excursion check for TIC 316937197."""
import json
import hashlib
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from urllib.request import urlopen

URL='http://asassn-lb01.ifa.hawaii.edu:9006/get_lightcurve/541166428682'
with urlopen(URL,timeout=60) as response:
    raw=response.read()
curve=json.loads(raw)['light_curve']
times=[]
cameras=defaultdict(list)
for values in curve['data']:
    row=dict(zip(curve['columns'],values))
    times.append(row['jd'])
    flux,error=row['flux'],row['flux_err']
    if row['quality']=='G' and flux and error and flux>0 and 0<error/flux<=.03:
        cameras[row['image_id'][:2]].append(flux)
depths=[]
for values in cameras.values():
    reference=median(values)
    depths.extend(1-flux/reference for flux in values)
summary=dict(url=URL,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),
             rows=len(times),first_jd=min(times),last_jd=max(times),
             sector_points=sum(61204<=jd-2400000.5<=61234 for jd in times),
             selected=len(depths),n_deeper_5pct=sum(x>.05 for x in depths),
             n_deeper_10pct=sum(x>.10 for x in depths),
             mad=median(abs(x-median(depths)) for x in depths))
Path(sys.argv[1]).write_text(json.dumps(summary,indent=2)+'\n')

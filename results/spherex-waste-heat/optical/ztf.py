"""Search final-year Hephaistos neighborhoods in public ZTF r-band light curves."""
import concurrent.futures
import csv
from collections import defaultdict
import datetime
import gzip
import hashlib
import io
import json
from pathlib import Path
from statistics import median
import sys
import time
from urllib.request import urlopen
from urllib.parse import urlencode
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'wise'))
from fetch import targets, position, separation

def inspect(t):
    ra,dec=position(t,2022)
    params=dict(POS=f'CIRCLE {ra:.7f} {dec:.7f} 0.0006944444',BANDNAME='r',FORMAT='csv',BAD_CATFLAGS_MASK=32768,COLLECTION='ztf_dr24')
    url='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?'+urlencode(params)
    cache=ROOT/'cache';cache.mkdir(exist_ok=True)
    path=cache/(t['label']+'.csv.gz')
    meta={'url':url,'label':t['label'],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    for attempt in range(2):
        try:
            if path.exists():
                body=gzip.decompress(path.read_bytes())
            else:
                with urlopen(url,timeout=150) as response: body=response.read()
            reader=csv.DictReader(io.StringIO(body.decode('utf-8-sig')))
            if 'mjd' not in (reader.fieldnames or []):
                if b'no records' in body.lower() or b'no data' in body.lower(): rows=[]
                else: raise ValueError(body[:360])
            else: rows=list(reader)
            path.write_bytes(gzip.compress(body,mtime=0))
            meta.update(sha256=hashlib.sha256(body).hexdigest(),bytes=len(body),rows=len(rows))
            break
        except Exception as e:
            if attempt==1:
                meta['error']=str(e)[:450]
                return meta,[]
            time.sleep(2)
    groups=defaultdict(list)
    for r in rows:
        if int(r['catflags'])!=0 or float(r['magerr'])>0.1:continue
        mjd=float(r['mjd']); when=2000+(mjd-51544.5)/365.25
        if separation(float(r['ra']),float(r['dec']),*position(t,when))>1.5:continue
        groups[r['oid']].append((mjd,float(r['mag']),float(r['magerr'])))
    result=[]
    for oid,pts in sorted(groups.items()):
        if len(pts)<5:continue
        nights=defaultdict(list)
        for mjd,mag,err in pts:nights[int(mjd)].append(mag)
        nightly=np.array([median(x) for x in nights.values()]); baseline=np.median(nightly)
        count=sum(x>=baseline+.2 for x in nightly)
        result.append(dict(label=t['label'],oid=oid,points=len(pts),nights=len(nights),mjd_min=min(nights),mjd_max=max(nights),median_r=round(baseline,4),p95_minus_median=round(np.percentile(nightly,95)-baseline,4),faintest_minus_median=round(max(nightly)-baseline,4),nights_fainter_by_0p2=int(count)))
    return meta,result

if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        values=list(pool.map(inspect,targets()))
    (ROOT/'query_receipts.json').write_text(json.dumps([r for r,_ in values],indent=2)+'\n')
    fields='label oid points nights mjd_min mjd_max median_r p95_minus_median faintest_minus_median nights_fainter_by_0p2'.split()
    with (ROOT/'r_lightcurves.csv').open('w') as out:
        w=csv.DictWriter(out,fields,lineterminator='\n');w.writeheader()
        for _,rows in values:w.writerows(rows)

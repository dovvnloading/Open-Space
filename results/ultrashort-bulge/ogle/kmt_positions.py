"""Check public KMT positional DIAPL data for independent measurements on OGLE peak nights."""
import argparse
import concurrent.futures as cf
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

import numpy as np
from counterparts import HERE,coords


BASE='https://kmtnet.kasi.re.kr/search/?'
TARGETS=['2025-BLG-1151','2025-BLG-1201','2025-BLG-0643']


def get(path,url):
    if not path.exists():
        with urllib.request.urlopen(url,timeout=120) as resp: data=resp.read()
        path.write_bytes(data)
    return path.read_bytes()


def one(task):
    cache,row=task
    event=row['event']; ra,dec=coords(row['ra_hms'],row['dec_dms']); t0=float(row['t0_hjd'])
    url=BASE+urllib.parse.urlencode(dict(search_type='radec',ra=f'{ra:.10f}',dec=f'{dec:.10f}',ra_type='decimal_degree',arcsec='2'))
    out=dict(event=event,query=url,query_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    try:
        html=get(cache/(event+'.html'),url)
        out.update(query_sha256=hashlib.sha256(html).hexdigest(),matches=[])
        markers=json.loads(re.search(r'const kmtnetSearchMarkers = (.*?);',html.decode(),re.S)[1])
        for mark in markers:
            item={k:mark.get(k) for k in ('lcid','ra','dec','dist')}
            lcid=mark['lcid']
            if not re.fullmatch(r'[A-Za-z0-9.]+',lcid):raise ValueError(lcid)
            link=BASE+urllib.parse.urlencode(dict(action='download_lc',lcid=lcid))
            item['url']=link; out['matches'].append(item)
            data=get(cache/(lcid+'.csv'),link)
            item['sha256']=hashlib.sha256(data).hexdigest()
            rows=list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))
            t=np.array([float(p['date'])+2400000 for p in rows])
            fl=np.array([float(p['flux']) for p in rows]); er=np.array([float(p['err']) for p in rows])
            site=np.array([p['site'] for p in rows]);yrs=np.array([p['year'] for p in rows])
            item.update(n=len(t),first_hjd=float(t.min()),last_hjd=float(t.max()),years=sorted(set(yrs.tolist())),sites=[])
            for s in sorted(set(site.tolist())):
                a=(site==s)&(yrs==event[:4])&np.isfinite(fl)&np.isfinite(er)&(er>0)
                local=a&(abs(t-t0)<.5); off=a&(abs(t-t0)>5)
                # Same-site seasonal baseline, flux zero point as downloaded.
                stat=None;base=None;sc=None
                if np.sum(off)>4:
                    base=float(np.median(fl[off]));sc=max(1.,float(np.median(abs((fl[off]-base)/er[off]))/.67449))
                    if local.any():stat=float(np.max((fl[local]-base)/(sc*er[local])))
                item['sites'].append(dict(site=s,n_season=int(a.sum()),n_halfday=int(local.sum()),
                    n_hour=int(np.sum(a&(abs(t-t0)<1/24))),peak_max_z=stat,baseline=base,scale=sc,
                    closest_dt_days=float((t[a]-t0)[np.argmin(abs(t[a]-t0))]) if a.any() else None,
                    local_epochs=[dict(dt_days=round(float(tt-t0),5),flux=float(ff),err=float(ee))
                                  for tt,ff,ee in zip(t[local],fl[local],er[local])]))
        out['status']='ok' if markers else 'no light-curve ID returned'
        out['search_center_returned']=bool(re.search(rb'const kmtnetSearchCenter = \{',html))
    except Exception as e:out['status']=str(e)
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);cache=ap.parse_args().cache
    cache.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((HERE/'catalog.csv').open()))
    with cf.ThreadPoolExecutor(3) as pool:
        out=list(pool.map(one,[(cache,r) for r in rows if r['event'] in TARGETS]))
    (HERE/'kmt-positions.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()

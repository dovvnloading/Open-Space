"""Compare the event position with nearby stars on the same field, filter and site."""
import argparse
import csv
import io
import json
import re
import hashlib
from pathlib import Path
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
import numpy as np

from freeze import HERE,get
from neighbours import angle

PORTAL='https://kmtnet.kasi.re.kr/search/?'
TARGETS={'KMT-2024-BLG-3326':60570.59,'KMT-2024-BLG-0456':60404.63,
         'KMT-2025-BLG-0220':60884.503,'KMT-2025-BLG-0300':60951.28,'KMT-2025-BLG-2510':60849.48}


def summarize(raw,center):
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    out=[]
    for site in ['BLG','SAO','SSO']:
        sub=[r for r in rows if r['site']==site]
        if not sub:continue
        a=np.array([[float(r[k]) for k in ('date','flux','err','see','year')] for r in sub]);a=a[np.isfinite(a).all(axis=1)&(a[:,2]>0)&(a[:,3]>0)]
        if not len(a):continue
        t,f,e,see,y=a.T
        # Separate calibration for every site/source/year. Suppress worst-seeing 5 percent.
        zz=np.full(len(t),np.nan)
        for yr in np.unique(y):
            sel=y==yr
            quiet=sel&(abs(t-center)>15)
            if quiet.sum()<8:continue
            hi=np.quantile(see[quiet],.95)
            valid=sel&(see<=hi)
            ctrl=quiet&valid
            floor=np.quantile(e[ctrl],.05)
            e0=np.maximum(e,floor);z=(f-np.median(f[ctrl]))/e0
            zz[valid]=z[valid]/max(1.,1.4826*np.median(abs(z[ctrl]-np.median(z[ctrl]))))
        good=np.isfinite(zz)
        peak=good&(abs(t-center)<.15)
        days=np.floor(t[good]).astype(int);n=np.bincount(days);hit=np.bincount(days,weights=zz[good]>5,minlength=len(n))
        repeats=np.flatnonzero((hit>=3)&(2*hit>=n)&(abs(np.arange(len(n))+.5-center)>15))
        out.append(dict(site=site,years='|'.join(map(str,map(int,np.unique(y)))),usable=int(good.sum()),first=round(t.min(),4),
             last=round(t.max(),4),n_near=int(peak.sum()),peak_median_z=round(float(np.median(zz[peak])),3) if peak.any() else None,
             n5_near=int((peak&(zz>5)).sum()),recurring_nights=len(repeats),recurrence_dates='|'.join(map(str,repeats))))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args();out=args.cache/'field-vet';out.mkdir(exist_ok=True)
    ix={r['event']:r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines())}
    records=[];rows=[]
    def search(url):
        key=hashlib.sha256(url.encode()).hexdigest()
        local=out/(key+'.html');rec=out/(key+'.receipt')
        if local.exists() and rec.exists():return local.read_bytes(),json.loads(rec.read_text())
        raw,r=get(url);local.write_bytes(raw);rec.write_text(json.dumps(r))
        return raw,r

    for event,center in TARGETS.items():
        r=ix[event];ra,dec=angle(r['ra'],True),angle(r['dec'])
        query=PORTAL+urlencode(dict(search_type='radec',ra=f'{ra:.9f}',dec=f'{dec:.9f}',ra_type='decimal_degree',arcsec=14))
        raw,rec=search(query);records.append(dict(event=event,kind='cone',**rec))
        m=re.search(rb'const kmtnetSearchMarkers = (\[.*?\]);',raw)
        if not m:continue
        markers=sorted(json.loads(m.group(1)),key=lambda x:x['dist'])
        if not markers:continue
        # Anchor on catalog field code; all controls must share the target detector/patch.
        field='BLG'+r['field'][3:];target=next((v for v in markers if v['lcid'].split('.')[0]==field),markers[0])
        field=target['lcid'].split('.')[0]
        same=[v for v in markers if v['lcid'].startswith(field+'.')]
        # Different projected positions beyond the usual stellar PSF and within 14 arcsec;
        # choose nearby and distant controls. Disclose duplicates in the crowded inner 1.5 arcsec.
        neighbours=[v for v in same if v['dist']<=1.5]
        controls=[v for v in same if v['dist']>=4]
        if len(controls)<2:
            # The portal truncates crowded cones to its first few matches.
            # Place separate cones 7 arcsec east/west for independent local controls.
            for offset in (-7.,7.):
                ura=ra+offset/(3600*np.cos(np.deg2rad(dec)))
                uq=PORTAL+urlencode(dict(search_type='radec',ra=f'{ura:.9f}',dec=f'{dec:.9f}',ra_type='decimal_degree',arcsec=2.))
                control_raw,control_rec=search(uq);records.append(dict(event=event,kind='displaced_cone',**control_rec))
                cm=re.search(rb'const kmtnetSearchMarkers = (\[.*?\]);',control_raw)
                if not cm:continue
                nearby=sorted([c for c in json.loads(cm.group(1)) if c['lcid'].startswith(field+'.')],key=lambda c:c['dist'])
                if nearby:
                    c=nearby[0]
                    c['dist']=3600*np.hypot((c['ra']-ra)*np.cos(np.deg2rad(dec)),c['dec']-dec)
                    if c['dist']>=4:controls.append(c)
        picks=neighbours[:2]+(controls[:2]+controls[-1:] if controls else [])
        unique={v['lcid']:v for v in picks}
        def download(v):
            lcid=v['lcid'];cached=out/(lcid+'.json')
            if cached.exists():return json.loads(cached.read_text())
            url=PORTAL+urlencode(dict(action='download_lc',lcid=lcid))
            try:
                data,receipt=get(url)
            except OSError as exc:
                return dict(receipt=dict(url=url,status=getattr(exc,'code',None),error=str(exc)),data=[])
            value=dict(receipt=receipt,data=summarize(data,center))
            cached.write_text(json.dumps(value)+'\n')
            return value
        with ThreadPoolExecutor(2) as pool:
            vals=list(pool.map(download,unique.values()))
        for marker,val in zip(unique.values(),vals):
            records.append(dict(event=event,kind='light_curve',lcid=marker['lcid'],**val['receipt']))
            rows.extend(dict(event=event,lcid=marker['lcid'],role='target' if marker['lcid']==target['lcid'] else 'neighbor',
                sep_arcsec=round(marker['dist'],4),**s) for s in val['data'])
    with (HERE/'field-controls.tsv').open('w') as h:
        w=csv.DictWriter(h,list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    (HERE/'field-control-receipts.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()

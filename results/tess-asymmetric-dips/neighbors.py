"""Fit each event in the six nearest QLP stars outside a 2 arcmin radius."""
import argparse
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from acquire import fetch, product_url
from pixels import flux_fit, read_stamp


def sky_vectors(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)], axis=-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('triggers','sample','out'): parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    sample=json.loads(args.sample.read_text())
    data=fetch(sample['url'])
    if hashlib.sha256(data).hexdigest()!=sample['sha256']:
        raise ValueError('roster changed')
    ids, ras, decs=np.loadtxt(io.BytesIO(data),delimiter=',',unpack=True)
    vectors=sky_vectors(ras,decs)
    jobs={}
    for ev in csv.DictReader(args.triggers.open()):
        jobs.setdefault(ev['tic'],[]).append(ev)
    tasks=[]
    for tic,evs in jobs.items():
        target=evs[0]
        v=sky_vectors(float(target['ra_obj']),float(target['dec_obj']))
        arcmin=np.degrees(np.arccos(np.clip(vectors@v,-1,1)))*60
        dist=np.where(arcmin>=2,arcmin,np.inf)
        for i in np.argsort(dist)[:6]:
            tasks.append((tic,evs,int(ids[i]),dist[i]))
    def run(item):
        tic,evs,neighbor,sep=item
        url=product_url(int(evs[0]['sector']),neighbor)
        r=dict(target=tic,neighbor=neighbor,sep_arcmin=sep,url=url)
        try:
            raw=fetch(url)
            h,_,tab=read_stamp(raw)
            if int(h['TICID'])!=neighbor:raise ValueError('wrong star')
            r.update(sha256=hashlib.sha256(raw).hexdigest(),camera=int(h['CAMERA']),ccd=int(h['CCD']))
            same=(int(evs[0]['camera']),int(evs[0]['ccd']))==(r['camera'],r['ccd'])
            r.update(same_detector=same)
            good=(tab['QUALITY']==0)&np.isfinite(tab['TIME'])
            t=tab['TIME'][good]
            flux=np.column_stack([tab['SAP_FLUX'][good],tab['DET_FLUX'][good],tab['SAP_BKG'][good]])
            r['events']=[]
            for ev in evs:
                peak=float(ev['peak_btjd'])
                fit=flux_fit(t,flux,peak,float(ev['duration_h'])/24,float(ev['ingress_fraction']))
                e=dict(peak_btjd=peak)
                if fit:
                    baseline,depth,count=fit
                    e.update(sap_depth=float(depth[0]/baseline[0]),det_depth=float(depth[1]/baseline[1]),
                             bkg_deficit=float(depth[2]),n=count)
                r['events'].append(e)
            r['status']='ok'
        except Exception as err:
            r.update(status='error',error=str(err))
        return r
    with ThreadPoolExecutor(5) as pool:
        rows=list(pool.map(run,tasks))
    args.out.write_text(json.dumps(rows,indent=1)+'\n')
    print(len(rows),sum(r['status']=='ok' for r in rows))


if __name__=='__main__':
    main()

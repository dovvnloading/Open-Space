"""Paired W1 tests at observed W2 positions, reading unTimely catalogs in memory."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import math
from pathlib import Path
import re
import time
from urllib.request import urlopen

from astropy.io import fits
import numpy as np
from scipy.spatial import cKDTree

BASE="https://portal.nersc.gov/project/cosmo/data/unwise/neo7/untimely-catalog/"

def fetch(url, cache=None):
    if cache and cache.exists():
        return cache.read_bytes()
    for attempt in range(3):
        try:
            with urlopen(url,timeout=90) as r:return r.read()
        except OSError:
            if attempt==2: raise
            time.sleep(1+attempt)

def table(payload):
    return fits.getdata(io.BytesIO(gzip.decompress(payload)),1)

def write(path,rows):
    with path.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run_tile(job):
    tile,tracks,idx,cache=job
    sel=idx[idx['COADD_ID']==tile]
    w1=sel[sel['BAND']==1];w2=sel[sel['BAND']==2]
    look={int(s['EPOCH']):s for s in w2}
    groups=defaultdict(list);receipts=[];matched_w2={}
    for track in tracks:
        choices=defaultdict(list)
        for name in track['ids'].split():
            epoch=int(re.search(r'e(\d+)$',name).group(1))
            wr=look[epoch];same=w1[w1['FORWARD']==wr['FORWARD']]
            if not len(same):continue
            w=same[np.argmin(abs(same['MJDMEAN']-wr['MJDMEAN']))]
            if abs(w['MJDMEAN']-wr['MJDMEAN'])<=3:
                choices[int(w['FORWARD'])].append((float(wr['MJDMEAN']),name,wr,w))
        for parity,options in choices.items():
            ordered=sorted(options,key=lambda row:row[0])
            # All same-direction synchronized visits of tracks with >=10 linked W2 epochs.
            count=len(ordered) if int(track['n'])>=10 else min(3,len(ordered))
            for j in np.linspace(0,len(ordered)-1,count,dtype=int):
                t,name,wr,w=ordered[j]
                groups[int(w['EPOCH'])].append((track,name,wr,w))
    result=[]
    for ep,targets in sorted(groups.items()):
        wn=targets[0][3]['CATALOG_FILENAME'].removeprefix('./')
        data=fetch(BASE+wn);receipts.append(dict(band=1,tile=tile,file=wn,sha256=hashlib.sha256(data).hexdigest()))
        a=table(data)
        ref_dec=float(targets[0][3]['DEC']);fx=math.cos(math.radians(ref_dec))*3600
        xy=np.column_stack([a['ra']*fx,a['dec']*3600]);tr=cKDTree(xy)
        for track,name,wr,w in targets:
            eid=int(wr['EPOCH'])
            if eid not in matched_w2:
                nm=wr['CATALOG_FILENAME'].removeprefix('./');d=fetch(BASE+nm,cache/Path(nm).name)
                receipts.append(dict(band=2,tile=tile,file=nm,sha256=hashlib.sha256(d).hexdigest()))
                rows=table(d);all_ids=set(q for tt in tracks for q in tt['ids'].split())
                s=rows[np.isin(rows['unwise_detid'].astype(str),list(all_ids))]
                matched_w2[eid]={str(v['unwise_detid']):tuple(float(v[k]) for k in ['ra','dec','flux','dflux','fracflux','rchi2','dx','dy']) for v in s}
            ra,de,f2,df2,ff2,c2,dx2,dy2=matched_w2[eid][name]
            dt=(float(w['MJDMEAN'])-float(wr['MJDMEAN']))/365.25
            ra+=float(track['pmra'])*dt/fx;de+=float(track['pmdec'])*dt/3600
            target=[ra*fx,de*3600]
            dist,ind=tr.query(target)
            neighbors=np.array(tr.query_ball_point(target,180))
            local_sep=np.linalg.norm(xy[neighbors]-target,axis=1)
            r=a[neighbors];good=(local_sep>=12)&(r['qf']>=.85)&(r['flags_unwise']==0)&(r['dflux']>0)&(r['flux']>4*r['dflux'])&(r['flux']<15*r['dflux'])&(r['fracflux']>=.6)
            ds=r['dflux'][good];noise=np.percentile(ds,90) if len(ds)>=10 else np.nan
            m=a[ind]
            result.append(dict(track=track['source_name'],tile=tile,w2_id=name,w1_epoch=ep,w1_mjd=float(w['MJDMEAN']),parity=int(w['FORWARD']),
                w2_mjd=float(wr['MJDMEAN']),ra=ra,dec=de,w2_flux=f2,w2_sigma=df2,w2_fracflux=ff2,w2_rchi2=c2,w2_dx_pix=dx2,w2_dy_pix=dy2,
                nearest_w1_arcsec=dist,w1_nearest_id=str(m['unwise_detid']),w1_flux=float(m['flux']),w1_sigma=float(m['dflux']),
                w1_fracflux=float(m['fracflux']),w1_rchi2=float(m['rchi2']),w1_flags=int(m['flags_unwise']),w1_qf=float(m['qf']),
                n_local=len(ds),local_sigma_p90=float(noise)))
    print(tile,'tracks',len(tracks),'pairs',len(result),flush=True)
    return result,receipts

def main():
    ap=argparse.ArgumentParser();ap.add_argument('tracks',type=Path);ap.add_argument('cache',type=Path);ap.add_argument('out',type=Path);a=ap.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    name='untimely_index-neo7.fits.gz';payload=fetch(BASE+name,a.cache/name);idx=table(payload)
    bytile=defaultdict(list)
    for t in csv.DictReader(a.tracks.open()):bytile[t['tile']].append(t)
    jobs=sorted(bytile,key=lambda tile:-max(int(t['n']) for t in bytile[tile]))
    outputs=[];receipts=[dict(band=0,tile='index',file=name,sha256=hashlib.sha256(payload).hexdigest())]
    with ThreadPoolExecutor(5) as pool:
        for result,sources in pool.map(run_tile,((n,bytile[n],idx,a.cache) for n in jobs)):
            outputs.extend(result);receipts.extend(sources)
            write(a.out/'paired_epochs.csv',outputs);write(a.out/'input_hashes.csv',receipts)
if __name__=="__main__":main()

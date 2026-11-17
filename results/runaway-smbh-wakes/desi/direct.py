"""Reconcile sky cones to the DESI DR1 iron HEALPix manifest; inspect covered redrock files."""
import argparse
import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import astropy.units as u
from astropy.io import fits
from astropy_healpix import HEALPix
import numpy as np

HERE=Path(__file__).resolve().parent
URL='https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/'

def cached(url,cache):
    dest=cache/url.rsplit('/',1)[-1]
    for n in range(3):
        try:
            if not dest.exists():
                with urlopen(url,timeout=75) as response: raw=response.read()
                if not raw.startswith(b'SIMPLE  '):raise OSError('not FITS')
                dest.write_bytes(raw)
            b=dest.read_bytes()
            return dest,dict(url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
        except OSError:
            if n==2:raise
            time.sleep(2**n)

def write(file,rows,fields):
    with file.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def run(cache,candidates,review):
    cache.mkdir(exist_ok=True,parents=True)
    hp=HEALPix(nside=64,order='nested',frame='icrs')
    file,index=cached(URL+'tilepix.fits',cache)
    with fits.open(file) as hd:
        t=hd['TILEPIX'].data
        index['manifest_rows']=len(t)
        mapping={}
        for pix,s,program in zip(t['HEALPIX'],t['SURVEY'],t['PROGRAM']):
            mapping.setdefault(int(pix),set()).add((str(s).strip(),str(program).strip(),int(pix)))
    selected={r['id'] for r in csv.DictReader(review.open())}
    targets=[]
    for r in csv.DictReader(candidates.open()):
        if r['id'] in selected:
            targets.append((r['id'],float(r['ra_deg']),float(r['dec_deg']),float(r['length_arcsec'])/2+5))
    if {r[0] for r in targets} != selected:
        raise ValueError('Endpoint catalogue does not cover every visual review ID')
    targets.append(('HST_XDF_west',53.18577725,-27.78513879,10.8))
    cones=[];products={}
    for name,ra,dec,radius in targets:
        hpix=sorted(set(map(int,hp.cone_search_lonlat(ra*u.deg,dec*u.deg,radius*u.arcsec))))
        groups=sorted({group for pix in hpix for group in mapping.get(pix,[])})
        for group in groups:products.setdefault(group,[]).append((name,ra,dec,radius))
        cones.append(dict(id=name,ra_deg=ra,dec_deg=dec,radius_arcsec=radius,healpix=' '.join(map(str,hpix)),
                          published_groups=len(groups),products_read=0,spectra=0,matches_on_sky=''))
    def download(item):
        (s,prog,pix),tg=item
        url=URL+f'{s}/{prog}/{pix//100}/{pix}/redrock-{s}-{prog}-{pix}.fits'
        matches=[]
        try:
            zfile,receipt=cached(url,cache)
            with fits.open(zfile) as h:
                z=h['REDSHIFTS'].data;fm=h['FIBERMAP'].data
                zlookup={int(k):i for i,k in enumerate(z['TARGETID'])}
                rak='MEAN_FIBER_RA' if 'MEAN_FIBER_RA' in fm.names else 'TARGET_RA'
                dk='MEAN_FIBER_DEC' if 'MEAN_FIBER_DEC' in fm.names else 'TARGET_DEC'
                a=fm[rak];d=fm[dk]
                receipt.update(rows=len(z),coord_columns=[rak,dk])
                for name,ra,dec,r in tg:
                    hav=np.sin(np.radians(d-dec)/2)**2+np.cos(np.radians(dec))*np.cos(np.radians(d))*np.sin(np.radians(a-ra)/2)**2
                    sep=np.degrees(2*np.arcsin(np.sqrt(np.clip(hav,0,1))))*3600
                    for ix in np.flatnonzero(sep<=r):
                        tid=int(fm['TARGETID'][ix]);q=z[zlookup[tid]]
                        matches.append(dict(id=name,targetid=tid,survey=s,program=prog,pix=pix,ra=float(a[ix]),dec=float(d[ix]),sep_arcsec=float(sep[ix]),
                                            z=float(q['Z']),zerr=float(q['ZERR']),zwarn=int(q['ZWARN']),spectype=str(q['SPECTYPE']).strip(),deltachi2=float(q['DELTACHI2'])))
            return receipt,[r[0] for r in tg],matches
        except Exception as e:return dict(url=url,error=repr(e)),[],[]
    with ThreadPoolExecutor(5) as pool:res=list(pool.map(download,sorted(products.items())))
    done={r['id']:r for r in cones}
    allmatches=[]
    for receipt,names,mm in res:
        for name in names:done[name]['products_read']+=1
        allmatches.extend(mm)
    for r in cones:
        mm=[m for m in allmatches if m['id']==r['id']]
        r['spectra']=len(mm)
        r['matches_on_sky']=' '.join(str(m['targetid']) for m in mm)
        r['status']='no DR1 HEALPix coverage' if not r['published_groups'] else 'complete' if r['products_read']==r['published_groups'] else 'incomplete'
    write(HERE/'direct_cones.csv',cones,list(cones[0]))
    write(HERE/'direct_matches.csv',allmatches,['id','targetid','survey','program','pix','ra','dec','sep_arcsec','z','zerr','zwarn','spectype','deltachi2'])
    (HERE/'direct_inputs.json').write_text(json.dumps({'tilepix':index,'redrock':[r[0] for r in res]},indent=2)+'\n')
    print('cones',len(cones),'products',len(res),'matched rows',len(allmatches),'errors',sum('error' in r[0] for r in res),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--candidates',type=Path,default=HERE.parent/'euclid/endpoint_candidates.csv')
    p.add_argument('--review',type=Path,default=HERE.parent/'detector/review_rows.csv')
    a=p.parse_args();run(a.cache,a.candidates,a.review)

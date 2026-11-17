"""Seeing-limited Legacy Surveys DR10/DR11 grz imaging at Euclid ridge positions."""
import argparse
import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

from acquire import sha,write_csv


def photometry(a,xy,scale):
    dp=xy[1]-xy[0];v=dp/np.linalg.norm(dp);perp=np.array([-v[1],v[0]])
    center=xy[0]+np.linspace(.15,.85,45)[:,None]*dp
    widths=np.r_[np.linspace(-.5,.5,5),np.linspace(-2.25,-1.25,5),np.linspace(1.25,2.25,5)]/scale
    on=center[:,None,:]+widths[None,:,None]*perp
    offsets=np.array([(i,j) for i in [-7,-5,-3,0,3,5,7] for j in [-7,-5,-3,0,3,5,7]
                      if abs(np.dot((i,j),perp))>2.6])/scale
    pix=on[None,:,:,:]+np.vstack([[0,0],offsets])[:,None,None,:]
    # A zero-valued survey hole is absent data even when finite.
    arr=np.where(np.isfinite(a)&(a!=0),a,np.nan)
    z=map_coordinates(arr,[pix[...,1].ravel(),pix[...,0].ravel()],order=1,mode='constant',cval=np.nan).reshape(pix.shape[:-1])
    dif=z[:,:,:5].mean(axis=2)-z[:,:,5:].mean(axis=2)
    flux=np.median(dif.reshape(len(dif),5,9).mean(axis=2),axis=1)
    blank=flux[1:][np.isfinite(flux[1:])]
    if len(blank)<12 or not np.isfinite(flux[0]):return dict(usable=0,n_blank=len(blank),contrast='',blank_mad='',sn='')
    med=np.median(blank);sig=1.4826*np.median(abs(blank-med))
    if sig<=0:return dict(usable=0,n_blank=len(blank),contrast='',blank_mad='',sn='')
    return dict(usable=1,n_blank=len(blank),contrast=float(flux[0]-med),blank_mad=float(sig),sn=float((flux[0]-med)/sig))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('candidates','metrics','cache','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--jobs',type=int,default=3);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True);a.cache.mkdir(exist_ok=True,parents=True)
    rows={r['id']:r for r in csv.DictReader(a.candidates.open())}
    ids=[r['id'] for r in csv.DictReader(a.metrics.open()) if r['sign']=='positive' and r['band']=='VIS' and
         r['usable']=='1' and int(r['positive_bins'])>=4 and float(r['central_ridge_contrast_sigma'])>=3.]
    def run(task):
        ident,release=task;r=rows[ident];layer='ls-'+release
        url='https://www.legacysurvey.org/viewer/fits-cutout?'+urlencode(dict(ra=f"{float(r['ra_deg']):.8f}",dec=f"{float(r['dec_deg']):.8f}",
            pixscale=.262,width=144,height=144,layer=layer,bands='grz'))
        fp=a.cache/(sha(url.encode())+'.fits');rec=dict(id=ident,release=release,url=url)
        try:
            if not fp.exists():
                for j in range(3):
                    try:
                        with urlopen(url,timeout=45) as res:b=res.read()
                        if not b.startswith(b'SIMPLE'):raise ValueError(b[:160].decode(errors='replace'))
                        fp.write_bytes(b);break
                    except (OSError,ValueError) as exc:
                        if j==2:raise
                        if getattr(exc,'code',None)==429:
                            h=exc.headers.get('Retry-After','30');time.sleep(int(h) if h.isdigit() else 30)
                        else: time.sleep(2**j)
            raw=fp.read_bytes();rec.update(bytes=len(raw),sha256=sha(raw))
            with fits.open(BytesIO(raw)) as hdus:
                h=hdus[0].header;cube=hdus[0].data;w=WCS(h,naxis=2);rec.update(survey=h['SURVEY'],version=h['VERSION'])
                xy=w.all_world2pix([[float(r[f'ra{i}_deg']),float(r[f'dec{i}_deg'])] for i in (0,1)],0)
                result=[dict(id=ident,release=release,band=h[f'BAND{k}'],**photometry(a,xy,.262)) for k,a in enumerate(cube)]
            return rec,result
        except (OSError,ValueError,KeyError) as e:
            rec['error']=repr(e)
            return rec,[]
    with ThreadPoolExecutor(a.jobs) as pool:
        groups=list(pool.map(run,[(i,rel) for i in ids for rel in ('dr11','dr10')]))
    (a.out/'legacy_queries.json').write_text(json.dumps([q for q,_ in groups],indent=1)+'\n')
    vals=[s for _,v in groups for s in v]
    if vals:write_csv(a.out/'legacy_photometry.csv',vals)


if __name__=='__main__':main()

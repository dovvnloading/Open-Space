"""Independent detector-frame aperture, peak-quality and centroid checks for NEXUS."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

S3='https://stpubdata.s3.amazonaws.com/'
OBS={7:2,8:4,9:6}
BAD=1|2|4|1024|2048
R=12
Y,X=np.mgrid[-R:R+1,-R:R+1]

def get_products(obs):
    paths=[]
    for visit in range(1,7):
        prefix=f'jwst/public/jw09263/jw09263{obs:03}{visit:03}/'
        url=S3+'?list-type=2&prefix='+prefix
        root=ET.fromstring(urlopen(url,timeout=60).read())
        if root.findtext('{*}IsTruncated')=='true': raise ValueError('paged listing')
        paths.extend(v.text for v in root.findall('.//{*}Key') if v.text.endswith('_cal.fits') and '_nrc' in v.text)
    return paths


def stamp(a,x,y):
    x,y=int(np.rint(x)),int(np.rint(y))
    return a[y-R:y+R+1,x-R:x+R+1]

def measure(s,e,q,fx,fy,rad,unit):
    rr=np.hypot(X-fx,Y-fy)
    core=rr<=rad
    ann=(rr>=7)&(rr<=11)&np.isfinite(s)&((q&3)==0)
    if ann.sum()<65: return None
    sky=np.median(s[ann])
    mask=core&np.isfinite(s)&np.isfinite(e)
    if mask.sum()!=core.sum(): return None
    f=s-sky
    peak=int(np.argmax(np.where(core,f,-np.inf)))
    # Fixed normalization, signed flux and DQ are retained even for flagged rows.
    sigma=1.4826*np.median(abs(f[ann]))
    ferr=max(np.sqrt(np.sum(e[core]**2)),sigma*np.sqrt(core.sum()))*unit
    w=np.maximum(f,0)*(rr<=3.)
    dx,dy=(np.nan,np.nan) if w.sum()==0 else ((X*w).sum()/w.sum()-fx,(Y*w).sum()/w.sum()-fy)
    return dict(flux_njy=f[core].sum()*unit,err_njy=ferr,dx_pix=dx,dy_pix=dy,
        core_dq_or=int(np.bitwise_or.reduce(q[core])),peak_dq=int(q.flat[peak]),peak_above_sky=f.flat[peak],
        bad_core_pixels=int(np.count_nonzero((q[core]&BAD)!=0)))

def angular_position(w,pos,h):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        x,y=w.all_world2pix(pos.ra_deg,pos.dec_deg,0,quiet=True)
    ok=np.isfinite(x)&np.isfinite(y)&(x>=R+2)&(y>=R+2)&(x<h['NAXIS1']-R-3)&(y<h['NAXIS2']-R-3)
    return x,y,ok

def run_one(arg):
    path,ep,band,pos,anchors,out=arg
    name=path.split('/')[-1]; base=out/name.removesuffix('.fits')
    if base.with_suffix('.json').exists(): return
    url=S3+path
    from urllib.request import Request
    request=Request(url,headers={'Range':'bytes=0-65535'})
    with urlopen(request,timeout=100) as f:
        payload=f.read(); assert f.status==206
    file=BytesIO(payload);pr=fits.Header.fromfile(file);h=fits.Header.fromfile(file)
    assert h['EXTNAME']=='SCI' and pr['FILTER'].lower()==band
    w=WCS(h)
    px,py,ok=angular_position(w,pos,h)
    meta=dict(epoch=ep,filter=band,product=name,url=url,header_sha256=hashlib.sha256(payload).hexdigest(),
              header_bytes=len(payload),positions=int(ok.sum()),mjd=pr['EXPMID'],exptime=pr['EFFEXPTM'])
    if ok.any():
        with urlopen(url,timeout=150) as f: raw=f.read()
        meta.update(sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
        with fits.open(BytesIO(raw),memmap=False) as hd:
            a=hd['SCI'].data;e=hd['ERR'].data;q=hd['DQ'].data
            unit=h['PIXAR_SR']*1e15
            ax,ay,aok=angular_position(w,anchors,h)
            offsets=[]
            for j in np.flatnonzero(aok):
                x,y=ax[j],ay[j];s=stamp(a,x,y)
                r=np.hypot(X-(x-np.rint(x)),Y-(y-np.rint(y)))
                bg=np.nanmedian(s[(r>7)&(r<11)])
                weights=np.where((r<3.5)&np.isfinite(s),np.maximum(s-bg,0),0)
                if weights.sum()>0:
                    ox=(weights*X).sum()/weights.sum()-(x-np.rint(x))
                    oy=(weights*Y).sum()/weights.sum()-(y-np.rint(y))
                    if np.hypot(ox,oy)<2: offsets.append((ox,oy))
            sx,sy=np.median(offsets,axis=0) if len(offsets)>=3 else (0.,0.)
            meta.update(anchors=len(offsets),offset_x=sx,offset_y=sy)
            px+=sx;py+=sy
            rows=[]
            for j in np.flatnonzero(ok):
                x,y=px[j],py[j]
                r=measure(*(stamp(im,x,y) for im in (a,e,q)),x-np.rint(x),y-np.rint(y),2.8 if band=='f200w' else 2.2,unit)
                if r:rows.append(dict(id=pos.iloc[j].id,epoch=ep,band=band,product=name,mjd=pr['EXPMID'],
                                      x_pix=x,y_pix=y,anchor_count=len(offsets),**r))
            if rows: pd.DataFrame(rows).to_csv(base.with_suffix('.csv'),index=False)
    base.with_suffix('.json').write_text(json.dumps(meta)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('positions',type=Path);p.add_argument('calibrators',type=Path)
    p.add_argument('reference',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--epochs',type=int,nargs='+',default=[9,8]);p.add_argument('--bands',nargs='+',default=['f444w','f200w'])
    p.add_argument('--jobs',type=int,default=6);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    pos=pd.read_csv(a.positions);anc=pd.read_csv(a.calibrators);w=WCS(fits.getheader(a.reference))
    anc['ra_deg'],anc['dec_deg']=w.all_pix2world(anc.x,anc.y,0)
    with ThreadPoolExecutor(a.jobs) as ex:
        for e in a.epochs:
            items=[]
            for path in get_products(OBS[e]):
                band='f444w' if 'long_cal' in path else 'f200w'
                if band in a.bands:items.append((path,e,band,pos,anc,a.out))
            list(ex.map(run_one,items))
    rows=[pd.read_csv(f) for f in sorted(a.out.glob('jw*.csv'))]
    if rows: pd.concat(rows).to_csv(a.out/'native.csv',index=False)
    (a.out/'receipts.json').write_text(json.dumps([json.loads(f.read_text()) for f in sorted(a.out.glob('jw*.json'))],
                                            separators=(',',':'))+'\n')

if __name__=='__main__':main()

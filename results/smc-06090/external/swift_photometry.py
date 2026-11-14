"""Aperture photometry with a common nearby astrometric reference and local exposure maps."""
from pathlib import Path
import csv,json,warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.time import Time
from common import D,SKY,hjd,shape
warnings.filterwarnings('ignore',module='astropy')
def centroid(img,pos,scale):
    y,x=np.indices(img.shape)
    for _ in range(3):
        rr=np.hypot(x-pos[0],y-pos[1])*scale
        mask=rr<3.3
        if not np.any(mask):return None
        flux=np.maximum(img[mask]-np.median(img[(rr>4)&(rr<7)]),0)
        if flux.sum()<50:return None
        pos=np.array([np.sum(flux*x[mask]),np.sum(flux*y[mask])])/flux.sum()
    return pos
def aper(img,emap,pos,scale,r):
    y,x=np.indices(img.shape)
    dist=np.hypot(x-pos[0],y-pos[1])*scale
    pop=dist<r;sky=(dist>=6.5)&(dist<10.5)
    if pop.sum()<.85*np.pi*(r/scale)**2 or sky.sum()<.85*np.pi*(10.5**2-6.5**2)/scale**2:return None
    exp=np.mean(emap[pop])
    if exp<=0 or np.min(emap[pop])<.85*exp or np.min(emap[sky])<.75*exp:return None
    apcounts=img[pop].sum();skyc=img[sky].sum()
    ratio=emap[pop].sum()/emap[sky].sum()
    return float((apcounts-ratio*skyc)/exp),float((max(apcounts,1)+ratio**2*max(skyc,1))**.5/exp),float(exp)
def run():
    refs=Table.read(D/'inputs/gaia-references.csv'); sky=[[SKY.ra.deg,SKY.dec.deg]]
    for s in ['4686449235706531584','4686449029548388352']:
        r=refs[np.array(refs['source_id']).astype(str)==s][0]
        sky.append([r['ra'],r['dec']])
    rows=[]
    for path in sorted((D/'swift_stamps').glob('*.json')):
        obs=json.loads(path.read_text())
        if not path.with_suffix('.exp.npz').exists():continue
        with np.load(path.with_suffix('.npz')) as images,np.load(path.with_suffix('.exp.npz')) as exps:
          for z in obs['planes']:
            key=str(z['idx']);name='exp'+key
            if name not in exps:continue
            image=images['image'+key];exp=exps[name]
            w=WCS(fits.Header.fromstring(z['wcs']));scale=abs(w.pixel_scale_matrix[0,0])*3600
            xy=w.all_world2pix(sky,0)-z['ci']+[48,48]
            center=centroid(image,xy[1],scale)
            if center is None:continue
            off=center-xy[1]
            if np.linalg.norm(off)*scale>3.4:continue
            data={}
            for lab,pt in zip(['target','ref1','ref2'],xy+off):
                for rad in [2.7,3.5,4.5]:
                    out=aper(image,exp,pt,scale,rad)
                    if out:data.update({f'{lab}_{rad}_{k}':v for k,v in zip(['rate','err','exposure'],out)})
            if 'target_3.5_rate' not in data or 'ref1_3.5_rate' not in data:continue
            mid=np.mean(Time([z['date_obs'],z['date_end']],scale='utc').mjd)
            rows.append(dict(obsid=obs['obsid'],hdu=int(key),start_utc=z['date_obs'],date_end_utc=z['date_end'],mjd_mid=mid,pixel_scale_arcsec=scale,shift_arcsec=float(np.linalg.norm(off)*scale),**data))
    if rows:
        h=hjd(np.array([r['mjd_mid'] for r in rows]));ph,d1,d2=shape(h)
        for i,r in enumerate(rows):
            r.update(hjd_mid=h[i],phase=ph[i],delta_primary_cycles=d1[i],delta_secondary_cycles=d2[i])
        with (D/'swift-photometry.csv').open('w') as f:
            keys=list(dict.fromkeys(k for r in rows for k in r))
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    print('valid',len(rows))
if __name__=='__main__':run()

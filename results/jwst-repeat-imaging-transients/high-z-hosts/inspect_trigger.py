"""Extract calibrated detector photometry and DQ masks for the one F200W high-redshift trigger."""
import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

HOST=25213
S3='https://stpubdata.s3.amazonaws.com/'
TAP='https://jwst.esac.esa.int/server/tap/sync?'


def measure(im,er,dq,x,y,scale,unit):
    r=round(.65/scale);ix,iy=round(x),round(y)
    if ix-r<0 or iy-r<0 or ix+r>=im.shape[1] or iy+r>=im.shape[0]:return {}
    sl=np.s_[iy-r:iy+r+1,ix-r:ix+r+1]; data=im[sl];err=er[sl];mask=dq[sl]
    yy,xx=np.mgrid[-r:r+1,-r:r+1];rr=scale*np.hypot(xx-(x-ix),yy-(y-iy))
    ok=np.isfinite(data)&np.isfinite(err)&(err>0)&((mask&3)==0)
    ann=(rr>=.36)&(rr<=.50)&ok
    if ann.sum()<.65*np.count_nonzero((rr>=.36)&(rr<=.50)):return {}
    sky=np.median(data[ann]);result={}
    for size in (12,18,24):
        inside=rr<=size/100;good=inside&ok
        result[f'coverage_{size}']=good.sum()/inside.sum()
        result[f'dq_or_{size}']=int(np.bitwise_or.reduce(mask[inside]))
        if np.all(ok[inside]):
            result[f'ap{size}_njy']=(data[inside]-sky).sum()*unit
            result[f'err{size}_njy']=np.sqrt((err[inside]**2).sum())*unit
    return result


def one(args):
    ob,targets=args
    prefix=ob.obs_id.split('_')[0];name=ob.obs_id+'_cal.fits'
    url=S3+f'jwst/public/{prefix[:7]}/{prefix}/{name}'
    with urlopen(url,timeout=150) as r:data=r.read()
    with fits.open(BytesIO(data)) as hd:
        pr=hd[0].header; sh=hd['SCI'].header;w=WCS(sh)
        im=hd['SCI'].data;er=hd['ERR'].data;dq=hd['DQ'].data
        n,m=im.shape;scale=np.sqrt(sh['PIXAR_A2']);unit=sh['PIXAR_SR']*1e15
        x,y=w.all_world2pix(targets.ra,targets.dec,0,quiet=True)
        xh,yh=x[targets.id==HOST][0],y[targets.id==HOST][0]
        # Relative to survey MSA coordinates, translate the WCS using high-S/N field objects.
        anchors=(targets['sample']=='control')&(targets.f444>.5)&np.isfinite(x)&np.isfinite(y)&(x>17)&(y>17)&(x<m-17)&(y<n-17)
        off=[];yy,xx=np.mgrid[-8:9,-8:9]
        for j in np.flatnonzero(anchors):
            ix,iy=int(round(x[j])),int(round(y[j]));p=im[iy-8:iy+9,ix-8:ix+9];rr=np.hypot(xx-(x[j]-ix),yy-(y[j]-iy))
            sky=np.nanmedian(p[rr>=6]); vv=np.where((rr<.15/scale)&np.isfinite(p),np.maximum(p-sky,0),0)
            if vv.sum()<=0:continue
            dx=(xx*vv).sum()/vv.sum()-(x[j]-ix);dy=(yy*vv).sum()/vv.sum()-(y[j]-iy)
            if np.hypot(dx,dy)*scale<.12:off.append([dx,dy])
        dx,dy=np.median(off,axis=0) if len(off)>=3 else (0,0)
        xh+=dx;yh+=dy
        row=dict(product=name,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),filter=pr['FILTER'],mjd=pr['EXPMID'],
             anchors=len(off),dx=dx,dy=dy,edge_arcsec=min(xh,yh,m-1-xh,n-1-yh)*scale)
        if np.isfinite(xh+yh):row.update(measure(im,er,dq,xh,yh,scale,unit))
        return row


def main(targetfile,out):
    targets=pd.read_csv(targetfile);o=targets[targets.id==HOST].iloc[0]
    q=f"SELECT TOP 2000 obs_id,t_min,t_max,calib_level,access_url,em_min,em_max,instrument_name FROM ivoa.ObsCore WHERE calib_level=2 AND instrument_name='NIRCAM/IMAGE' AND dataproduct_type='image' AND t_min>60700 AND 1=CONTAINS(POINT('ICRS',{o.ra},{o.dec}),s_region) ORDER BY t_max DESC"
    raw=urlopen(TAP+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q)),timeout=90).read()
    obs=pd.read_csv(BytesIO(raw))
    # Cover both matched filters plus other independent spectral constraints from the Wide epoch.
    with ThreadPoolExecutor(max_workers=3) as pool: rows=list(pool.map(one,[(o,targets) for o in obs.itertuples()]))
    out.mkdir(exist_ok=True,parents=True)
    pd.DataFrame(rows).sort_values('mjd').to_csv(out/'detector_photometry.csv',index=False,float_format='%.10g')
    (out/'archive_inventory.csv').write_bytes(raw)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('targets',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args();main(args.targets,args.out)

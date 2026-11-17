"""Compare QLP and independent FFI extractions at disputed events, in 30-minute bins."""
import csv
import io
import json
import hashlib
import sys
import zipfile
from urllib.parse import urlencode

import numpy as np
from acquire import fetch,product_url
from pixels import read_stamp,target_pixel


def reduce(row):
    tic=int(row['tic'])
    sector=int(row['sector'])
    peak=float(row['peak_btjd'])
    url=product_url(sector,tic)
    lc=fetch(url)
    _,_,tab=read_stamp(lc)
    pxurl='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(dict(
        ra=row['ra_obj'],dec=row['dec_obj'],sector=sector,x=9,y=9,units='px'))
    with zipfile.ZipFile(io.BytesIO(fetch(pxurl))) as z:
        source=z.read(next(n for n in z.namelist() if n.endswith('.fits')))
    _,e,px=read_stamp(source)
    xy=target_pixel(e,float(row['ra_obj']),float(row['dec_obj']))
    yy,xx=(v.ravel() for v in np.indices((9,9)))
    border=(yy==0)|(xx==0)|(yy==8)|(xx==8)
    centered=(abs(xx-round(xy[0]))<=1)&(abs(yy-round(xy[1]))<=1)
    f=px['FLUX'].astype(float)
    sky=np.nanmedian(f[:,border],axis=1)
    f-=sky[:,None]
    ds=[(tab['TIME'],tab['QUALITY'],{k:tab[k] for k in ('SAP_FLUX','DET_FLUX','DET_FLUX_SML','DET_FLUX_LAG','SAP_BKG')}),
        (px['TIME'],px['QUALITY'],{'pix1':f[:,int(np.argmin((xx-xy[0])**2+(yy-xy[1])**2))],
                                   'pix3':np.sum(f[:,centered],axis=1),'perimeter':sky})]
    records=[]
    for j in range(-96,97):
        center=peak+j/48
        rec={'btjd':center}
        for t,q,cols in ds:
            good=(q==0)&(abs(t-center)<1/96)
            for k,y in cols.items():
                finite=good&np.isfinite(y)
                rec[k]=float(np.median(y[finite])) if finite.sum()>=4 else None
        records.append(rec)
    return dict(tic=tic,peak_btjd=peak,qlp_url=url,qlp_sha256=hashlib.sha256(lc).hexdigest(),
                pixel_url=pxurl,pixel_flux_sha256=hashlib.sha256(px['FLUX'].tobytes()).hexdigest(),
                aperture_xy=xy.tolist(),bins=records)


if __name__=='__main__':
    targets={316937197:4217.632,441246218:4208.340,2027133756:4227.965,116658378:4207.195}
    rows=csv.DictReader(open(sys.argv[1]))
    print(json.dumps([reduce(r) for r in rows if abs(float(r['peak_btjd'])-targets.get(int(r['tic']),0))<.01],
                     separators=(',',':'),allow_nan=False))

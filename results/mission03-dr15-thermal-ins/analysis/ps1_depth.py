"""Measure the J180019 field's r depth from a linear Pan-STARRS1 stack cutout."""
import csv
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.signal import fftconvolve

RA, DEC, R99 = 270.081101, 8.084437, 3.151
FX = 1.711016839692758e-14
def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    base='https://ps1images.stsci.edu/cgi-bin/'
    index=base+'ps1filenames.py?'+urlencode(dict(ra=RA,dec=DEC,filters='r',type='stack'))
    with urlopen(index,timeout=45) as r:rows=r.read().decode().splitlines()
    name=rows[1].split()[rows[0].split().index('filename')]
    url=base+'fitscut.cgi?'+urlencode(dict(ra=RA,dec=DEC,size=480,format='fits',red=name))
    path=cache/'j180019-ps1.fits'
    if not path.exists():
        with urlopen(url,timeout=45) as r:path.write_bytes(r.read())
    query='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv?'+urlencode(dict(ra=RA,dec=DEC,radius=.025,columns='[objID,raStack,decStack,primaryDetection,nDetections,rPSFMag,rPSFMagErr,rKronMag]'))
    with urlopen(query,timeout=45) as r: cat=list(csv.DictReader(io.StringIO(r.read().decode())))
    dat,h=fits.getdata(path,header=True);w=WCS(h).celestial
    x,y=w.world_to_pixel_values(RA,DEC)
    scale=np.mean([v.value for v in w.proj_plane_pixel_scales()])*3600
    yy,xx=np.mgrid[-22:23,-22:23];radius=np.hypot(yy,xx)*scale
    def kernel(lo,hi):return ((radius>=lo)&(radius<hi)).astype(float)
    aperture,annulus=kernel(0,1.25),kernel(2.75,4.5)
    weights=aperture-annulus*aperture.sum()/annulus.sum()
    sm=fftconvolve(np.nan_to_num(dat),weights,mode='same')
    mask=fftconvolve(np.isfinite(dat).astype(float),(radius<4.5).astype(float),mode='same')>(radius<4.5).sum()-.5
    Y,X=np.indices(dat.shape);dist=np.hypot(X-x,Y-y)*scale
    grid=(X%12==0)&(Y%12==0)&(dist>7)&(dist<51)&mask
    calibrators=[]
    for q in cat:
        if q['primaryDetection']!='1':continue
        a,b=w.world_to_pixel_values(float(q['raStack']),float(q['decStack']))
        if not np.isfinite([a,b]).all():continue
        rmag=float(q['rPSFMag']);err=float(q['rPSFMagErr']);k=float(q['rKronMag'])
        if 0<rmag<23 and 0<err<0.2:
            grid &= np.hypot(X-a,Y-b)*scale>3
        i,j=round(float(a)),round(float(b))
        if (16.5<rmag<20 and 0<err<.05 and abs(rmag-k)<.12 and int(q['nDetections'])>=3 and 24<i<456 and 24<j<456 and mask[j,i] and sm[j,i]>0):
            calibrators.append(rmag+2.5*np.log10(sm[j,i]))
    blank=sm[grid];med=np.median(blank);sig=1.4826*np.median(abs(blank-med));zp=np.median(calibrators)
    region=(dist<=R99)
    peak=(sm[region].max()-med)/sig
    lim=zp-2.5*np.log10(5*sig)
    return dict(ra=RA,dec=DEC,r99_arcsec=R99,image_url=url,cat_query=query,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),pixel_scale_arcsec=scale,all_apertures_in_circle_valid=bool(mask[region].all()),n_calibrators=len(calibrators),zeropoint_median=zp,zeropoint_mad=np.median(abs(np.array(calibrators)-zp)),blank_apertures=len(blank),sky_ap_sigma=sig,unshifted_central_sn=(sm[round(float(y)),round(float(x))]-med)/sig,max_sn_circle=peak,five_sigma_depth_ab=lim,fx_over_nufnu=FX/(2.99792458e18/6231*3.631e-20*10**(-.4*lim)))
if __name__=='__main__':
    print(json.dumps(run(Path(sys.argv[1])),indent=2))

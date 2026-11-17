"""Metadata-only SO/PHI geometry prerequisites; no science array is read."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from sunpy.coordinates import HeliographicStonyhurst, Helioprojective, propagate_with_solar_surface

REQUIRED=('DATE-BEG','DATE-AVG','DATE-END','TIMESYS','BUNIT','DSUN_OBS','HGLN_OBS','HGLT_OBS','RSUN_REF','CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','DSUN_AU','PXBEG1','PXBEG2','NAXIS1','NAXIS2','FILENAME')

def read_header(path):
    """Read primary FITS header only; never touch primary image data."""
    h=fits.getheader(path,0)
    absent=[k for k in REQUIRED if k not in h]
    if absent:raise ValueError('missing FITS metadata: '+', '.join(absent))
    if h['TIMESYS']!='UTC' or h['BUNIT'].strip().lower()!='gauss':raise ValueError('unexpected time or field unit')
    if (h['CTYPE1'],h['CTYPE2'])!=('HPLN-TAN','HPLT-TAN'):raise ValueError('unsupported projection')
    if not (Time(h['DATE-BEG'])<=Time(h['DATE-AVG'])<=Time(h['DATE-END'])):raise ValueError('invalid timing')
    return h

def observer(h):
    return SkyCoord(lon=h['HGLN_OBS']*u.deg,lat=h['HGLT_OBS']*u.deg,radius=h['DSUN_OBS']*u.m,frame=HeliographicStonyhurst,obstime=h['DATE-AVG'])

def to_surface(h,x,y):
    """Pixels are zero-based; return near-side solar surface or NaNs off disk."""
    w=WCS(h,naxis=2);tx,ty=w.all_pix2world(x,y,0)
    obs=observer(h);frame=Helioprojective(observer=obs,obstime=h['DATE-AVG'],rsun=h['RSUN_REF']*u.m)
    return SkyCoord(Tx=tx*u.deg,Ty=ty*u.deg,frame=frame).make_3d()

def mu(h,coordinate):
    """Local surface-normal cosine toward actual spacecraft, not B_r."""
    hgs=coordinate.transform_to(HeliographicStonyhurst(obstime=h['DATE-AVG']));point=hgs.cartesian.xyz.to_value(u.m)
    vec=observer(h).cartesian.xyz.to_value(u.m)-point
    return np.sum(point*vec,axis=0)/(np.linalg.norm(point,axis=0)*np.linalg.norm(vec,axis=0))

def transport(h1,h2,x,y,rotation_model='howard'):
    """Project a surface point between observers with differential rotation."""
    c=to_surface(h1,x,y);target=Helioprojective(observer=observer(h2),obstime=h2['DATE-AVG'],rsun=h2['RSUN_REF']*u.m)
    with propagate_with_solar_surface(rotation_model=rotation_model):new=c.transform_to(target)
    w=WCS(h2,naxis=2)
    x2,y2=w.all_world2pix(new.Tx.to_value(u.deg),new.Ty.to_value(u.deg),0)
    return x2,y2,new.is_visible()

def dust_center_zero_based(h):
    """Release team formula returns FITS one-based pixel coordinates."""
    d=h['DSUN_AU'];x=(-.23*d*d-1.19*d+1305.08)-h['PXBEG1']+2
    y=(5.96*d*d-10.74*d+975.62)-h['PXBEG2']+2
    return x-1,y-1

def metadata(path):
    h=read_header(path); digest=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    keys=REQUIRED+('OBS_ID','PARENT','VERSION','ACCRAWIT','PC1_1','PC1_2','PC2_1','PC2_2')
    return {'sha256':digest,'primary_header':{k:h.get(k) for k in keys},'dust_center_zero_based':dust_center_zero_based(h)}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('fits');a=p.parse_args();print(json.dumps(metadata(a.fits),indent=2))

"""Download two ASKAP and one GMRT survey products and repeat the radio-image checks.
Requires astropy, numpy, scipy. Files are retrieved by explicit public product URL.
"""
from pathlib import Path
from urllib.request import urlopen
import hashlib
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares

from io import BytesIO
import csv

ROOT=Path(__file__).resolve().parent
def get(row):
    body=urlopen(row['url'],timeout=90).read()
    if hashlib.sha256(body).hexdigest()!=row['sha256']:
        raise ValueError(f"Product changed: {row['survey']}")
    return fits.open(BytesIO(body))[0]

def askap(hdus):
    a,d=273.4162741379,-21.2195075583
    ims=np.array([h.data.squeeze() for h in hdus])
    w=WCS(hdus[0].header).celestial
    ox,oy=w.world_to_pixel_values(a,d)
    sx,sy=w.world_to_pixel_values(273.416125,-21.2167777778)
    def sample(x,y):
        coords=[np.atleast_1d(y),np.atleast_1d(x)]
        return np.array([map_coordinates(im,coords,order=1) for im in ims])
    y,x=np.indices(ims[0].shape)
    mask=(x>12)&(x<500)&(y>12)&(y<500)&(np.hypot(x-ox,y-oy)>20)
    bg=ims[:,mask]
    rms=1.4826*np.nanmedian(abs(bg-np.nanmedian(bg,axis=1)[:,None]),axis=1)
    a0,b0=sample(ox,oy)[:,0],sample(sx,sy)[:,0]
    rng=np.random.default_rng(2046)
    cx,cy=rng.uniform(13,499,130000),rng.uniform(13,499,130000)
    keep=np.hypot(cx-ox,cy-oy)>26
    cx,cy=cx[keep],cy[keep]
    ctrl=sample(cx,cy)
    displaced=sample(cx+sx-ox,cy+sy-oy)
    print('ASKAP latest and alternate positions, mJy/beam',1000*a0,1000*b0)
    print('Local robust image RMS, mJy/beam',1000*rms)
    print('Fraction exceeding both at same sky offset',np.mean(np.all(ctrl>=a0[:,None],axis=0)))
    print('Fraction exceeding the two-position difference',np.mean(np.sum((ctrl-displaced)/rms[:,None],axis=0)>=np.sum((a0-b0)/rms)))

def tgss(hdu):
    im=hdu.data.squeeze()*1000
    w=WCS(hdu.header).celestial
    yy,xx=np.indices(im.shape)
    a,b=w.pixel_to_world_values(xx,yy)
    ra,dec=296.021834,-17.836415
    u=(a-ra)*np.cos(np.deg2rad(dec))*3600
    v=(b-dec)*3600
    near=np.hypot(u,v)<30
    def residual(p):
        x=u[near]-p[0]
        y=v[near]-p[1]
        return p[2]+p[3]*np.exp(-.5*((x/11)**2+(y/16)**2))-im[near]
    p=least_squares(residual,[0,0,0,25]).x
    print('TGSS sensitivity fit RA, Dec (deg), peak (mJy/beam)',ra+p[0]/(3600*np.cos(np.deg2rad(dec))),dec+p[1]/3600,p[3])

if __name__=='__main__':
    with (ROOT/'image_products.csv').open() as stream:
        images=[get(r) for r in csv.DictReader(stream)]
    askap(images[:2])
    tgss(images[2])

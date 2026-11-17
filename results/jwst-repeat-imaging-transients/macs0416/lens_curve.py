"""Nearest zero of the tangential Jacobian eigenvalue in the HFF CATS v4.1 map."""
from io import BytesIO
from urllib.request import urlopen
from hashlib import sha256
from pathlib import Path
import json

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from contourpy import contour_generator
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
RA,DEC=64.03667859265462,-24.06729999780779
ZL,ZS=0.396,0.9397
ROOT="https://archive.stsci.edu/pub/hlsp/frontier/macs0416/models/cats/v4.1/"
STEM="hlsp_frontier_model_macs0416_cats_v4.1_"

def fetch_map(kind):
    url=ROOT+STEM+kind+".fits"
    with urlopen(url,timeout=40) as r: raw=r.read()
    with fits.open(BytesIO(raw)) as h:
        return h[0].data.copy(),WCS(h[0].header,naxis=2),{"url":url,"sha256":sha256(raw).hexdigest()}

def separation(sky):
    vertices=(sky-[RA,DEC])*[3600*np.cos(np.deg2rad(DEC)),3600.]
    p=vertices[:-1]
    d=np.diff(vertices,axis=0)
    lensq=np.sum(d*d,axis=1)
    valid=lensq>0
    if not valid.any(): return np.inf
    p,d,lensq=p[valid],d[valid],lensq[valid]
    f=np.clip(-np.sum(p*d,axis=1)/lensq,0,1)
    return np.hypot(*(p+f[:,None]*d).T).min()

def main():
    kappa,w,r1=fetch_map('kappa')
    gamma,w2,r2=fetch_map('gamma')
    if kappa.shape!=gamma.shape or not w.wcs.compare(w2.wcs): raise ValueError("Map WCS differs")
    def einv(z): return 1/np.sqrt(.3*(1+z)**3+.7)
    beta=quad(einv,ZL,ZS)[0]/quad(einv,0,ZS)[0]
    eig=1-beta*(kappa+gamma)
    curves=contour_generator(z=eig).lines(0)
    min_sep=min(separation(w.all_pix2world(line,0)) for line in curves if len(line)>1)
    res=dict(model="HFF CATS v4.1",ra=RA,dec=DEC,lens_z=ZL,source_z=ZS,omega_m=.3,beta=beta,
             model_pixel_arcsec=(3600*proj_plane_pixel_scales(w)).tolist(),dist_to_tangential_curve_arcsec=float(min_sep),
             inputs=[r1,r2])
    (HERE/"warhol_lens.json").write_text(json.dumps(res,indent=2)+"\n")
    print(round(min_sep,3))
if __name__=="__main__":main()

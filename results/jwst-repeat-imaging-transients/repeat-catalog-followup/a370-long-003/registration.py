"""Test positional motion relative to nearby stationary field sources, 2022 to Jan 2026."""
from pathlib import Path
from io import BytesIO
from urllib.request import urlopen
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from lightcurve import RA,DEC,product_url

def read(product):
 raw=urlopen(product_url(product),timeout=90).read()
 a=pd.read_csv(BytesIO(raw),sep=r'\s+',comment='#',usecols=['sky_centroid.ra','sky_centroid.dec','aper_total_flux','aper_total_flux_err'])
 a.columns=['ra','dec','flux','err']
 xy=(a[['ra','dec']].to_numpy()-[RA,DEC])*[3600*np.cos(np.radians(DEC)),3600]
 keep=np.isfinite(xy).all(axis=1)&(a.flux/a.err>12)&(np.hypot(*xy.T)<60)
 return xy[keep],a.flux.to_numpy()[keep]
def run(b):
 a,fa=read('jw01208-o002_t001_nircam_clear-'+b);bb,fb=read('jw05324-o014_t026_nircam_clear-'+b)
 dist,ii=cKDTree(bb).query(a);_,ij=cKDTree(a).query(bb)
 ok=(dist<.18)&(ij[ii]==np.arange(len(a)))&(np.hypot(*a.T)>1)&(fa>1e-7)
 off=bb[ii]-a; anc=off[ok]; drift=np.median(anc,axis=0);scat=1.4826*np.median(abs(anc-drift),axis=0)
 j=np.argmin(np.hypot(*a.T));d=off[j]-drift
 return dict(filter=b,controls=len(anc),control_dx=drift[0],control_dy=drift[1],scatter_x=scat[0],scatter_y=scat[1],rel_dx=d[0],rel_dy=d[1],motion_arcsec=np.hypot(*d))
if __name__=='__main__':
 from concurrent.futures import ThreadPoolExecutor
 with ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(run,['f115w','f150w','f200w','f277w','f356w','f444w']))
 pd.DataFrame(r).to_csv('registration.csv',index=False,float_format='%.7g')

"""Measure field-star centroid offsets for the same integrated-Gaussian model."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import ndtr
from measure import load_cube

def measure(path,catalog):
    t,cube,wcs=load_cube(path)
    im=np.median(cube[(t>3288)&(t<3290)],axis=0)
    tab=pd.read_csv(catalog,dtype={"source_id":str})
    ra=tab.ra+7.95*tab.pmra.fillna(0)/3.6e6/np.cos(np.deg2rad(tab.dec))
    dec=tab.dec+7.95*tab.pmdec.fillna(0)/3.6e6
    tab["x"],tab["y"]=wcs.world_to_pixel_values(ra,dec)
    results=[]
    for s in tab.itertuples():
        if not (1.5<s.x<14.5 and 1.5<s.y<14.5 and s.phot_g_mean_mag<16.5):
            continue
        dist=np.hypot(s.x-tab.x,s.y-tab.y)
        near=tab[(dist<3.5)&(dist>.05)].phot_g_mean_mag.min()
        if near<s.phot_g_mean_mag+1.5:continue
        xx,yy=np.meshgrid(np.arange(round(s.x)-2,round(s.x)+3),np.arange(round(s.y)-2,round(s.y)+3))
        image=im[yy,xx]
        def resid(v):
            a,x,y,sx,sy,b,bx,by=v
            px=ndtr((xx+.5-x)/sx)-ndtr((xx-.5-x)/sx)
            py=ndtr((yy+.5-y)/sy)-ndtr((yy-.5-y)/sy)
            return (a*px*py+b+bx*(xx-s.x)+by*(yy-s.y)-image).ravel()
        v=least_squares(resid,[np.ptp(image)*4,s.x,s.y,.65,.65,np.min(image),0,0],
            bounds=([0,s.x-1,s.y-1,.2,.2,-np.inf,-np.inf,-np.inf],
                    [np.inf,s.x+1,s.y+1,2,2,np.inf,np.inf,np.inf])).x
        results.append(dict(source_id=s.source_id,mag=s.phot_g_mean_mag,dx_pixel=v[1]-s.x,dy_pixel=v[2]-s.y))
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fits');p.add_argument('gaia')
    p.add_argument('--output',default='field-calibration.json');a=p.parse_args()
    Path(a.output).write_text(json.dumps(measure(a.fits,a.gaia),indent=2)+'\n')

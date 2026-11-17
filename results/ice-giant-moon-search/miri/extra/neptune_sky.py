"""Measure native-image centroids and sky coordinates of the single-visit F1000W points."""
import argparse,csv
from pathlib import Path
import numpy as np
from astropy.time import Time
from scipy.optimize import least_squares
from coverage import OUT,frame,get_states
from fetch import read_product,products
from miri_audit import write_csv


def measure(data,x,y):
    half=7;ix,iy=round(x),round(y);gy,gx=np.mgrid[-half:half+1,-half:half+1];gx=gx+ix-x;gy=gy+iy-y
    cut=np.s_[iy-half:iy+half+1,ix-half:ix+half+1];image=data['SCI'][cut].astype(float);err=data['ERR'][cut];dq=data['DQ'][cut]
    ok=np.isfinite(image)&np.isfinite(err)&(err>0)&((dq==0)|(dq==4))&(np.hypot(gx,gy)<=7)
    anchor=(np.hypot(gx,gy)>4.5)&ok;med=np.median(image[anchor]);noise=max(1e-6,1.4826*np.median(abs(image[anchor]-med)))
    # Fit a fixed F1000W Gaussian plus sky plane. Soft-L1 guards isolated high pixels.
    valid=ok;gx=gx[valid];gy=gy[valid];values=image[valid];sig=1.3
    def resid(p):
        dx,dy,amp,b,sx,sy=p
        return (amp*np.exp(-.5*((gx-dx)**2+(gy-dy)**2)/sig**2)+b+sx*gx+sy*gy-values)/noise
    p0=[0.,0.,max(noise,np.max(values)-med),med,0.,0.]
    fit=least_squares(resid,p0,loss='soft_l1',bounds=([-3.5,-3.5,0.,-np.inf,-np.inf,-np.inf],[3.5,3.5,np.inf,np.inf,np.inf,np.inf]))
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1.,2*fit.cost/(len(values)-6));unc=np.sqrt(np.maximum(0,cov.diagonal()))
    return x+fit.x[0],y+fit.x[1],unc[0],unc[1],fit.x[2],np.hypot(*fit.x[:2])


def run(cache,head):
    obs=list(csv.DictReader((OUT/'neptune_alias_representatives.csv').open()));states=get_states()
    fs={f['file']:(f,p) for f,p in products() if f['visit']=='01604005001'}
    result=[]
    for family,representative in [('inner',151978.),('outer',271086.)]:
        rows=sorted([r for r in obs if r['visit']=='01604005001' and float(r['a_km'])==representative],key=lambda r:r['file'])
        for row in rows:
            f,p=fs[row['file']];g=frame(f,p,head,states[f['file']],{});dat,_=read_product(p,cache,True)
            x,y,sx,sy,amp,disp=measure(dat,float(row['x']),float(row['y']))
            ra,dec=g['w'].all_pix2world([[x,y]],0)[0];t=Time(float(f['exp_mid_mjd_utc']),format='mjd',scale='utc')
            result.append(dict(feature=family,file=row['file'],mjd_mid_utc=t.mjd,mid_utc=t.isot+'Z',native_x_zero_based=x,native_y_zero_based=y,
                               ra_icrs_deg=ra,dec_icrs_deg=dec,formal_fit_x_pix=sx,formal_fit_y_pix=sy,peak_mjy_sr=amp,
                               offset_from_trial_pix=disp,aperture_mJy=float(row['aperture_mJy']),wcs_frame=g['w'].wcs.radesys))
    write_csv(OUT/'neptune_single_visit_astrometry.csv',list(result[0]),result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('headers',type=Path);a=p.parse_args();run(a.cache,a.headers)

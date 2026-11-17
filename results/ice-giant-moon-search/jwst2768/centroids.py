"""Fit individual native CAL tiles near a specified phase-shifted Perdita trajectory."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import ARCSEC,sky_xy
from prepare import inverse_sky
from search import load_npz
from perdita import rotate_tracks


def fit_blob(xs,ys,data,err,dq,sigma):
    x0,y0=xs[len(xs)//2],ys[len(ys)//2];height,width=data.shape
    if x0<6 or y0<6 or x0>width-7 or y0>height-7:
        return None
    x1,y1=int(np.rint(x0))-6,int(np.rint(y0))-6
    crop=np.s_[y1:y1+13,x1:x1+13];dy,dx=np.indices((13,13))
    dx=dx+x1-x0;dy=dy+y1-y0
    raw=data[crop];er=err[crop];q=dq[crop]
    valid=np.isfinite(raw)&np.isfinite(er)&(er>0)&((q==0)|(q==4))
    if valid.sum()<80 or np.mean(valid[np.hypot(dx,dy)<2.8])<.65:
        return None
    x,y=dx[valid],dy[valid];f=raw[valid];unc=er[valid]
    tx,ty=xs-x0,ys-y0;wt=np.ones(len(tx));wt[[0,-1]]=.5;wt/=wt.sum()
    def profile(a,b):
        arr=np.exp(-((x[:,None]-a-tx)**2+(y[:,None]-b-ty)**2)/(2*sigma**2))
        return arr@wt
    def residual(par):
        a,b,amplitude,offset,sx,sy=par
        return (amplitude*profile(a,b)+offset+sx*x+sy*y-f)/unc
    b0=np.median(f);amp=max(float(np.max(f[np.hypot(x,y)<3])-b0),float(np.median(unc)))
    fit=least_squares(residual,[0,0,amp,b0,0,0],bounds=([-2.5,-2.5,0,-np.inf,-np.inf,-np.inf],[2.5,2.5,np.inf,np.inf,np.inf,np.inf]),
                      loss='soft_l1',f_scale=2.)
    if not fit.success or np.max(np.abs(fit.x[:2]))>=2.4:
        return None
    scatter=max(1.,1.4826*np.median(np.abs(fit.fun-np.median(fit.fun))))
    try:
        cov=np.linalg.inv(fit.jac.T@fit.jac)*scatter**2
    except np.linalg.LinAlgError:
        return None
    errors=np.sqrt(np.diag(cov));sn=fit.x[2]/errors[2]
    if sn<5 or np.max(errors[:2])>.8:
        return None
    return x0+fit.x[0],y0+fit.x[1],float(sn),errors[:2],float(scatter),int(valid.sum())


def run(file,receipt,tiles,out,offset,body=725):
    cube=load_npz(file); band=str(cube['band']).split('/')[1]; info=json.loads(Path(receipt).read_text())
    centers={r['exposure_id']:r for r in info['registrations']}
    paths=sorted(Path(tiles).glob('*.fits'));ids=set(r['exposure_id'] for r in info['tile_inputs'])
    xx,yy=rotate_tracks(cube,np.array([offset]),body);xx=xx[0];yy=yy[0]
    results=[]
    for i,g in enumerate(cube['exposure_id']):
        c=centers[g];mid=xx.shape[1]//2; sc=float(cube['scale_arcsec'])
        e=-(xx[i]-cube['center_xy'][0])*sc/ARCSEC;n=(yy[i]-cube['center_xy'][1])*sc/ARCSEC
        for path in paths:
            if path.stem.rsplit('_',2)[0]!=g or path.stem not in ids:continue
            with fits.open(path,memmap=False) as hdul:
                sci=hdul['SCI'];hdr=sci.header;w=WCS(hdr,naxis=2)
                cra,cdec=inverse_sky(hdr['MT_RA'],hdr['MT_DEC'],c['east_shift_arcsec']/ARCSEC,c['north_shift_arcsec']/ARCSEC)
                wr,wd=inverse_sky(cra,cdec,e,n);xs,ys=w.all_world2pix(wr,wd,0,quiet=True)
                res=fit_blob(xs,ys,sci.data,hdul['ERR'].data,hdul['DQ'].data,float(cube['psf_sigma_px']))
                if res is None:continue
                x,y,snr,xy_unc,scatter,count=res
                tr,td=w.all_pix2world(x,y,0)
                sky=sky_xy(tr,td,cra,cdec,sc,cube['center_xy'])
                er=-(sky[0]-cube['center_xy'][0])*sc/ARCSEC
                no=(sky[1]-cube['center_xy'][1])*sc/ARCSEC
                ra,dec=inverse_sky(cube['path_ra'][i,mid],cube['path_dec'][i,mid],er,no)
                nominal=cube['known_xy'][i,list(cube['known_id']).index(body)]
                sep=np.linalg.norm(cube['known_xy'][i]-sky,axis=1);sep[list(cube['known_id']).index(body)]=np.inf;k=int(np.argmin(sep))
                # Approximate native pixel size is 31 mas. Includes local formal
                # fit uncertainty only. Limb/WCS registration uncertainty is separate.
                results.append(dict(frame=path.stem,band=band,visit=str(cube['visit'][i]),mjd_utc=float(cube['mjd'][i]),
                    ra_deg=float(ra),dec_deg=float(dec),native_x=float(x),native_y=float(y),sky_x=float(sky[0]),sky_y=float(sky[1]),
                    nominal_offset_x_px=float(sky[0]-nominal[0]),nominal_offset_y_px=float(sky[1]-nominal[1]),
                    fitted_snr=snr,formal_sigma_native_x_px=float(xy_unc[0]),formal_sigma_native_y_px=float(xy_unc[1]),
                    residual_error_scale=scatter,n_valid_native=count,nearest_other_id=int(cube['known_id'][k]),nearest_other_px=float(sep[k])))
    if not results:
        raise ValueError('no independent native centroid above S/N 5')
    with open(out,'w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=results[0]);wr.writeheader();wr.writerows(results)
    print(band,len(results),{v:sum(r['visit']==v for r in results) for v in set(q['visit'] for q in results)})


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','receipt','tiles','output']:a.add_argument(k)
    a.add_argument('--offset',type=float,required=True)
    a.add_argument('--body',type=int,default=725)
    p=a.parse_args();run(p.cube,p.receipt,p.tiles,p.output,p.offset,p.body)

"""Fit complete 20-pixel F140M trails around a predetermined longitude offset."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
base=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(base/'jwst2768'));sys.path.insert(0,str(base))
from perdita import rotate_tracks
from prepare import inverse_sky
from geometry import ARCSEC,sky_xy
from search import load_npz


def fit(x,y,img,err,dq):
    mid=len(x)//2;cx,cy=x[mid],y[mid]
    x0,y0=int(np.floor(min(x)))-6,int(np.floor(min(y)))-6
    x1,y1=int(np.ceil(max(x)))+7,int(np.ceil(max(y)))+7
    if x0<0 or y0<0 or x1>img.shape[1] or y1>img.shape[0]:return None
    gy,gx=np.mgrid[y0:y1,x0:x1];sl=np.s_[y0:y1,x0:x1]
    mask=np.isin(dq[sl],[0,4])&np.isfinite(img[sl]+err[sl])&(err[sl]>0)
    gx,gy,flux,unc=gx[mask],gy[mask],img[sl][mask],err[sl][mask]
    dx=gx-cx;dy=gy-cy
    weights=np.ones(len(x));weights[[0,-1]]=.5;weights/=weights.sum()
    def model(par):
        xx,yy,amplitude,sigma,off,sx,sy=par
        profile=(np.exp(-((gx[:,None]-x-xx)**2+(gy[:,None]-y-yy)**2)/(2*sigma*sigma))*weights).sum(1)
        return amplitude*profile+off+sx*dx+sy*dy
    init=[0,0,max(.1,10*np.percentile(flux-np.median(flux),95)),.8,np.median(flux),0,0]
    sol=least_squares(lambda p:(model(p)-flux)/unc,init,loss='soft_l1',f_scale=2.,
      bounds=([-4,-4,0,.55,-np.inf,-np.inf,-np.inf],[4,4,np.inf,1.5,np.inf,np.inf,np.inf]))
    p=sol.x;scatter=max(1.,1.4826*np.median(abs(sol.fun-np.median(sol.fun))))
    cov=np.linalg.pinv(sol.jac.T@sol.jac)*scatter**2
    return cx+p[0],cy+p[1],p[0],p[1],p[3],p[2]/np.sqrt(cov[2,2]),np.sqrt(cov[0,0]),np.sqrt(cov[1,1]),len(gx)


def run(path,files,out):
    cube=load_npz(path);sc=float(cube['scale_arcsec']);ref=cube['center_xy'];mid=cube['path_mjd'].shape[1]//2
    xx,yy=rotate_tracks(cube,np.array([-19.8]));angles=np.arange(-25,-15,.01);rx,ry=rotate_tracks(cube,angles)
    rows=[]
    for i,r in enumerate(json.loads(Path(str(path)+'.inputs.json').read_text())):
        with fits.open(files/r['filename']) as h:
            w=WCS(h['SCI'].header,naxis=2);cr,cd=w.all_pix2world(r['center_native_x'],r['center_native_y'],0)
            ra,dec=inverse_sky(cr,cd,-(xx[0,i]-ref[0])*sc/ARCSEC,(yy[0,i]-ref[1])*sc/ARCSEC)
            x,y=w.all_world2pix(ra,dec,0);q=fit(x,y,h['SCI'].data,h['ERR'].data,h['DQ'].data)
            if q is None:continue
            cx,cy,dx,dy,sigma,snr,ex,ey,n=q
            ra,dec=w.all_pix2world(cx,cy,0);sx,sy=sky_xy(ra,dec,cr,cd,sc,ref)
            ll=angles[np.argmin(np.hypot(rx[:,i,mid]-sx,ry[:,i,mid]-sy))]
            ra,dec=inverse_sky(cube['path_ra'][i,mid],cube['path_dec'][i,mid],-(sx-ref[0])*sc/ARCSEC,(sy-ref[1])*sc/ARCSEC)
            ok=snr>=5 and abs(dx)<3.9 and abs(dy)<3.9 and .56<sigma<1.49
            rows.append(dict(frame=r['exposure_id'],visit=cube['visit'][i],mjd_utc=cube['mjd'][i],accepted=ok,
                ra_deg=float(ra),dec_deg=float(dec),x_sky=float(sx),y_sky=float(sy),native_x=cx,native_y=cy,dx=dx,dy=dy,
                fitted_sigma_px=sigma,amplitude_formal_snr=snr,sigma_x_px=ex,sigma_y_px=ey,valid=n,longitude_deg=ll))
    with out.open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
    for r in rows:print(r['visit'],r['accepted'],*[round(r[k],2) for k in ['dx','dy','fitted_sigma_px','amplitude_formal_snr','longitude_deg']],flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for n in ['cube','files','out']:a.add_argument(n,type=Path)
    p=a.parse_args();run(p.cube,p.files,p.out)

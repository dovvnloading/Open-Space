"""Robustly measure the compact feature near the fixed October trial on each RATEINTS plane."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from link_visits import search
import positive_control as pc
from independent_vet import save
ROOT=Path(__file__).resolve().parent


def fit_patch(frame,x0,y0,guess,vx,vy):
    """Integrate the PSF over a linear short track, with a local background plane."""
    ix,iy=np.rint(np.array(guess)-[x0,y0]).astype(int)
    yy,xx=np.mgrid[-7:8,-7:8]
    cut=frame[iy-7:iy+8,ix-7:ix+8].astype(float);valid=np.isfinite(cut)
    back=valid&(xx**2+yy**2>30)
    if back.sum()<20 or valid.sum()<80:return None
    m=np.median(cut[back]);noise=max(1e-5,np.median(abs(cut[back]-m))*1.4826)
    weights=np.array([1,2,2,2,1])/8;tf=np.linspace(-.5,.5,5)
    def residual(p):
        x,y,amp,sig,b,gx,gy=p
        profile=sum(w*np.exp(-((xx-x-vx*t)**2+(yy-y-vy*t)**2)/(2*sig**2)) for w,t in zip(weights,tf))
        return ((b+gx*xx+gy*yy+amp*profile-cut)/noise)[valid]
    start=[0,0,max(noise,np.nanpercentile(cut,95)-m),1.,m,0,0]
    sol=least_squares(residual,start,bounds=([-4.5,-4.5,0,.55,-np.inf,-np.inf,-np.inf],[4.5,4.5,np.inf,2.5,np.inf,np.inf,np.inf]),loss='soft_l1',f_scale=2.5,max_nfev=150)
    x,y,amp,sig,*_=sol.x
    return dict(x=ix+x+x0,y=iy+y+y0,psf_sigma=sig,amplitude_DN_s=amp,background_mad=noise,
                model_peak_contrast=amp/noise,valid_pixels=int(valid.sum()),edge_parameter=int(abs(x)>4.4 or abs(y)>4.4 or sig>2.45 or sig<.57))


def run(args):
    rec=list(csv.DictReader((ROOT/'integration_check/integration_positions.csv').open()))
    c=search.load_npz(args.cube);rows=[]
    for j,f in enumerate(c['exposure_id']):
        m=pc.metadata(args.cache,str(f).removesuffix('_cal'));im,x0,y0=pc.wide_images(args.cache,m)
        selected=[r for r in rec if r['file']==m['file'] and r['body']=='trial']
        xy=np.array([[float(r['x']),float(r['y'])] for r in selected]);patches=np.array(pc.patches_for_positions(im,x0,y0,{'trial':xy})['trial'])
        co=np.nanmedian(patches,axis=0);yy,xx=np.mgrid[-pc.RADIUS:pc.RADIUS+1,-pc.RADIUS:pc.RADIUS+1]
        use=(np.isfinite(patches).sum(0)>=5)&(np.hypot(xx,yy)<20);iy,ix=np.unravel_index(np.nanargmax(np.where(use,co,np.nan)),co.shape); shift=np.array([xx[iy,ix],yy[iy,ix]])
        speed=np.gradient(xy,axis=0)*m['sec']/np.mean(np.diff(m['times'])*86400)
        for k,(guess,v) in enumerate(zip(xy+shift,speed)):
            p=fit_patch(im[k],x0,y0,guess,*v)
            if p is not None:rows.append(dict(j=j,k=k,t=m['times'][k],**p))
    save(args.output,rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('cube','cache','output'):p.add_argument(k,type=Path)
    run(p.parse_args())

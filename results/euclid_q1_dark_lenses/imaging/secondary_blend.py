"""Competing unlensed fit: an intrinsic curved feature plus one independent compact knot."""
import argparse
import csv
import numpy as np
from scipy.optimize import least_squares,nnls
from scipy.signal import fftconvolve
from cutouts import CutoutArchive,image_extensions
from measure import HERE,table,stamp
from secondary_fit import ID,profiles

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--index-cache");arg=p.parse_args()
    arc=CutoutArchive("lens",arg.index_cache);raw=image_extensions(arc.members(ID,[ID+".fits"])[ID+".fits"])
    im=raw['VIS_FLUX'][1][125:175,125:175];err=raw['VIS_RMS'][1][125:175,125:175]
    q=raw['VIS_PSF'][1].astype(float);q/=q.sum()
    yy,xx=np.indices(im.shape);yy=yy-24.5;xx=xx-24.5
    mask=(np.hypot(yy,xx)<21)&(err>0)&(err<1e15)&np.isfinite(im)
    dat=im[mask];w=1/err[mask];avg=np.average(dat,weights=w*w)
    init=np.fromstring(next(r['params'] for r in csv.DictReader((HERE/'secondary_models.csv').open()) if r['model']=='intrinsic_arc'),sep=' ')
    def render(t):
        fg,src,knot=0.,0.,0.
        for dx in (-.25,.25):
            for dy in (-.25,.25):
                x,y=xx+dx,yy+dy;f,s=profiles(t[:12],"intrinsic_arc",x,y);fg+=f/4;src+=s/4
                u=(x-t[12])*np.cos(t[16])+(y-t[13])*np.sin(t[16])
                v=-(x-t[12])*np.sin(t[16])+(y-t[13])*np.cos(t[16])
                knot+=np.exp(-.5*((u/t[14])**2+(v/t[15])**2))/4
        comps=np.stack([fftconvolve(z,q,mode='same') for z in (fg,src,knot)])
        mat=comps[:,mask].T;mean=np.average(mat,axis=0,weights=w*w)
        amps,_=nnls((mat-mean)*w[:,None],(dat-avg)*w)
        return np.einsum("i,ijk->jk",amps,comps)+(avg-mean@amps)
    def resid(t):return (render(t)[mask]-dat)*w
    lo=[-7,-7,1,.4,.18,-3.15,-8,-8,2,-6.3,.6,.08,-18,-18,.4,.4,-3.15]
    hi=[7,7,25,5,1,3.15,8,8,20,6.3,6,2.2,18,18,7,7,3.15]
    fits=[]
    for x,y in [(0,3),(-8,-6),(-8,5)]:
        fits.append(least_squares(resid,np.r_[init,x,y,1.3,1.3,0],bounds=(lo,hi),x_scale='jac',max_nfev=130))
    f=min(fits,key=lambda a:a.cost);mid=render(f.x)
    table(HERE/'secondary_blend.csv',[dict(model='intrinsic arc + independent Gaussian knot',diagonal_chi2=2*f.cost,
         usable_pixels=mask.sum(),n_free=len(f.x)+4,converged=bool(f.success),nfev=f.nfev,params=" ".join(f"{v:.8g}" for v in f.x))])
    stamp([("VIS",im),("arc + knot",mid),("residual",im-mid)],ID+"_blend",HERE)


if __name__=='__main__':main()

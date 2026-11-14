"""Check all phase-displaced templates, alternative sky extraction, and orbital time groups."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import median_filter
from astropy.io import fits
from astropy.wcs import WCS
import test_tess as tt

def simple(s):
    with fits.open(f's{s:03}.fits') as h:
        d=h[1].data
        t=d['TIME'].astype(float);img=d['FLUX'].astype(float)
        xy=WCS(h[2].header,naxis=2).world_to_pixel_values(tt.RA,tt.DEC)
        good=(d['QUALITY']==0)&np.isfinite(t)&np.all(np.isfinite(img),axis=(1,2));t,img=t[good],img[good]
    yy,xx=np.indices((11,11));x,y=map(float,xy)
    r=np.hypot(xx-x,yy-y);sky=(r>3.1)&(np.median(img,axis=0)<np.median(np.median(img,axis=0)[r>3.1]))
    bg=np.median(img[:,sky],axis=1)
    out=[]
    for ap in ((abs(xx-round(x))<=1)&(abs(yy-round(y))<=1),r<=np.sort(r.ravel())[3],r==r.min()):
        f=img[:,ap].sum(axis=1)-ap.sum()*bg
        v=np.zeros_like(f)
        for seg in np.split(np.arange(len(t)),np.flatnonzero(np.diff(t)>.1)+1):
            if len(seg)<12:continue
            width=max(5,int(.62/np.median(np.diff(t[seg])))//2*2+1)
            v[seg]=f[seg]-median_filter(f[seg],size=min(width,len(seg)//2*2+1),mode='nearest')
        out.append(v)
    flux=np.array(out).T
    keep=np.all(abs(flux-np.median(flux,axis=0))<8*tt.md(flux),axis=1)
    t=t[keep];flux=flux[keep];p=((t+2457000-70/86400-tt.T0)/tt.P)%1
    cycle=np.floor((t+2457000-70/86400-tt.T0)/tt.P+.17).astype(int)
    ifs=[]
    for n in np.unique(cycle):
        ix=cycle==n
        ids=np.floor(p[ix]*144).astype(int)%144;b=np.bincount(ids,minlength=144)
        if (b>0).sum()<.78*144 or (b@tt.tmpl((np.arange(144)+.5)/144)).min()<7: continue
        a=np.zeros((144,3));np.add.at(a,ids,flux[ix]);ifs.append((a,b))
    ar=[np.stack(w) for w in zip(*ifs)]; v=tt.evalf(ar);return {'n':len(ifs),'centers':v[5:7].tolist(),'sep':v[6]-v[5],'depths':v[3:5].tolist()}

def phase_control(ifs):
    ar=[np.stack(w) for w in zip(*ifs)]
    a,b=ar[0].sum(axis=0),ar[1].sum(axis=0)
    p=(np.arange(144)+.5)/144;w=np.sqrt(b)
    dep=[]
    for d in np.arange(144)/144:
        D=np.column_stack([np.ones(144),np.cos(2*np.pi*p),np.sin(2*np.pi*p),-tt.tmpl(p,tt.C+d)])
        dep.append(np.linalg.lstsq(D*w[:,None],a/w[:,None],rcond=None)[0][-2:])
    ds=np.array(dep)
    score=np.min(ds,axis=1)
    return {'rank_of_undisplaced_by_min_depth':(1+(score>score[0]).sum(axis=0)).tolist(),'phase_of_best_min_depth':(np.argmax(score,axis=0)/144).tolist(),'n_shifts':144}

if __name__=='__main__':
    out={}
    for s in map(int,sys.argv[1:]):
        r=np.load(f'folds-{s}.npz');ifs=list(zip(r['flux'],r['count']))
        out[s]={'phase_control':phase_control(ifs),'alternative_sky':simple(s)}
    print(json.dumps(out,indent=2))

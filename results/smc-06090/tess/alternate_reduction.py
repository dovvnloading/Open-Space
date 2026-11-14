"""Independent local-baseline photometry without principal components or block regression."""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter
from tess_timing import GROUPS,fit

HERE=Path(__file__).resolve().parent
def reduce(path):
    with fits.open(path) as h:
        t=h[1].data['TIME'].astype(float);quality=h[1].data['QUALITY']
        cube=h[1].data['FLUX'].astype(float)
        x,y=WCS(h[2].header,naxis=2).world_to_pixel_values(21.06275,-73.32594444)
    ok=np.isfinite(t)&(quality==0)&np.isfinite(cube).all(axis=(1,2))
    t=t[ok];cube=cube[ok];dt=np.median(np.diff(t))
    yy,xx=np.indices(cube.shape[1:]);d=np.hypot(xx-x,yy-y)
    median=np.median(cube,axis=0)
    candidates=(d>=3.5);sky=candidates&(median<=np.quantile(median[candidates],.6))
    aper=(np.abs(xx-round(float(x)))<=1)&(np.abs(yy-round(float(y)))<=1)
    y=cube[:,aper].sum(axis=1)-aper.sum()*np.median(cube[:,sky],axis=1)
    segments=np.r_[0,1+np.where(np.diff(t)>.15)[0],len(t)]
    flux=y.copy()
    for a,b in zip(segments[:-1],segments[1:]):
        window=2*int(.4/dt)+1
        flux[a:b]-=median_filter(y[a:b],size=window,mode='nearest')
    med=np.median(flux);scatter=1.4826*np.median(np.abs(flux-med))
    sel=np.abs(flux-med)<7*scatter;flux=flux[sel];t=t[sel]
    phase=(t+2457000.-70/86400.-2460000.0369)/1.1696719
    cyc=np.floor(phase+.2);binid=np.floor(120*(phase%1)).astype(int)
    ss=[];nn=[]
    for c in np.unique(cyc):
        mask=cyc==c;counts=np.bincount(binid[mask],minlength=120)
        if counts.sum()*dt<=.65*1.1696719 or (counts>0).sum()<min(55,int(.65*1.1696719/dt)):continue
        ss.append(np.bincount(binid[mask],weights=flux[mask], minlength=120));nn.append(counts)
    return np.array(ss),np.array(nn)
if __name__=='__main__':
    r={}
    for group,sects in GROUPS.items():
        sets=[reduce(HERE/'raw'/f's{s:03}.fits') for s in sects];S,N=map(list,zip(*sets))
        q=fit(S,N)
        r[group]={'usable_cycles':sum(map(len,S)),'centers_and_spacing':q[:3].tolist(),'depths_e_s':q[3:5].tolist()}
        print(group,q[:3].round(4),flush=True)
    (HERE/'alternate_reduction.json').write_text(json.dumps(r,indent=2)+'\n')

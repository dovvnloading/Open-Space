"""Contrast the flux-growth curve of MSAID 25213 with same-epoch stars and nearby sky."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from acquire_nexus import R
ID=25213
Y,X=np.mgrid[-R:R+1,-R:R+1]
UNIT=84.6159499407523

def aper(a,x,y,r,ann):
    rad=np.hypot(X-x,Y-y)
    w=np.zeros_like(rad)
    for dx in (np.arange(5)+.5)/5-.5:
        for dy in (np.arange(5)+.5)/5-.5:w+=np.hypot(X+dx-x,Y+dy-y)<=r
    keep=w>0;bg=a[(rad>=ann[0])&(rad<=ann[1])]
    if not np.isfinite(a[keep]).all() or np.isfinite(bg).sum()<.8*len(bg):return np.nan
    return UNIT*np.sum(w[keep]*(a[keep]-np.nanmedian(bg)))/25

def run(cache,out):
    t=pd.read_csv(cache/'targets.csv');j=t.index[t.id==ID][0]
    epcal=pd.read_csv(out.parent/'f200w/epoch_calibration.csv')
    # At least 0.6 arcsec from the target, sufficiently far from the cutout edge, 0.36 arcsec sampling.
    offsets=[(0,0)]+[(x,y) for x in range(-18,19,6) for y in range(-18,19,6) if np.hypot(x,y)>=10]
    rows=[];starrows=[];stars=np.flatnonzero((t.z==0)&(t.f444>.45)&(t['sample']=='control'))
    epochs=[5,6]
    for k in range(8):
        inp=np.load(cache/f'ep{k+1:02d}_data.npz');a=inp['stamps']
        x=inp['px']-np.rint(inp['px'])+epcal.dx[k];y=inp['py']-np.rint(inp['py'])+epcal.dy[k]
        for st in stars:
            v=[aper(a[st],x[st],y[st],r,(8,11)) for r in (2,6)]
            if np.all(np.isfinite(v)) and min(v)>0:starrows.append((k+1,t.id[st],*v,v[0]/v[1]))
        for ox,oy in offsets:
            vals=[aper(a[j],x[j]+ox,y[j]+oy,r,(8,11))/epcal.scale[k] for r in (2,6)]
            rows.append((k+1,ox,oy,*vals))
    s=pd.DataFrame(starrows,columns=['epoch','id','small','large','growth_ratio'])
    ratio=float(s[s.epoch.isin(epochs)].growth_ratio.median())
    r=pd.DataFrame(rows,columns=['epoch','ox','oy','small','large'])
    sel=r.dropna().groupby(['ox','oy']).filter(lambda x:len(x)==8)
    g=sel.groupby(['ox','oy','epoch']).mean().reset_index()
    a=g[g.epoch.isin(epochs)].groupby(['ox','oy'])[['small','large']].mean()-g[~g.epoch.isin(epochs)].groupby(['ox','oy'])[['small','large']].mean()
    a['excess_small_minus_psf']=a.small-ratio*a.large
    a['star_ratio']=ratio
    ctrl=a.drop((0,0)).excess_small_minus_psf
    spread=1.4826*np.median(abs(ctrl-ctrl.median()))
    a['centered_local_score']=(a.excess_small_minus_psf-ctrl.median())/spread
    a.to_csv(out/'radial_contrast.csv',float_format='%.7g');s.to_csv(out/'star_growth.csv',index=False,float_format='%.7g')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args();run(args.cache,args.out)

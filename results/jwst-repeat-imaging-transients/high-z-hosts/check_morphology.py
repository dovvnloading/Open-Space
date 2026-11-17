"""Repeat the triggered host aperture measurements across radii and display two-band image residuals."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from acquire_nexus import R, DATES

ID=25213
Y,X=np.mgrid[-R:R+1,-R:R+1]
UNIT=84.6159499407523


def run(root,tfile,out):
    t=pd.read_csv(tfile);idx=t.index[t.id==ID][0]; all_rows=[];fig,axes=plt.subplots(4,8,figsize=(13,6.5))
    for b,(band,cache) in enumerate(root):
        imgs=[];fxs=[];fys=[]
        for k in range(8):
            sci=np.load(cache/f'ep{k+1:02d}_data.npz');img=sci['stamps'][idx]
            er=np.load(cache/f'ep{k+1:02d}_error.npz')['stamps'][idx]
            fx=sci['px'][idx]-np.rint(sci['px'][idx]);fy=sci['py'][idx]-np.rint(sci['py'][idx])
            imgs.append(img);fxs.append(fx);fys.append(fy)
            rr=np.hypot(X-fx,Y-fy)
            for r in (2.,2.5,3.,3.5,4.,5.,6.):
                for lo,hi in [(6,8.3),(8,11),(10,14)]:
                    sky=np.nanmedian(img[(rr>=lo)&(rr<=hi)]);w=np.zeros_like(rr)
                    for xx in (np.arange(5)+.5)/5-.5:
                        for yy in (np.arange(5)+.5)/5-.5:w+=np.hypot(X-fx+xx,Y-fy+yy)<=r
                    w=w/25;good=(w>0)&np.isfinite(img)&(er>0)
                    if np.any((w>0)&(~good)):continue
                    f=UNIT*np.sum((img[good]-sky)*w[good]);e=UNIT*np.sqrt(np.sum((er[good]*w[good])**2))
                    all_rows.append(dict(epoch=k+1,band=band,radius_arcsec=.06*r,sky_inner_arcsec=.06*lo,sky_outer_arcsec=.06*hi,flux_njy=f,pipeline_sigma_njy=e))
        imgs=np.array(imgs);med=np.nanmedian(imgs,axis=0);di=imgs-med
        core=np.s_[R-12:R+13,R-12:R+13]
        vm=np.nanpercentile(np.abs(di[:,core[0],core[1]]),98.5);v1,v2=np.nanpercentile(imgs[:,core[0],core[1]],[5,99.3])
        for k in range(8):
            for row,plane,cmap,lo,hi in [(0,imgs[k],'gray',v1,v2),(1,di[k],'RdBu_r',-vm,vm)]:
                ax=axes[b*2+row,k];ax.imshow(plane[core],origin='lower',cmap=cmap,vmin=lo,vmax=hi,extent=[-.75,.75,-.75,.75])
                ax.add_patch(plt.Circle((fxs[k]*.06,fys[k]*.06),.18,fill=False,lw=.55,color='lime'))
                ax.set_xticks([]);ax.set_yticks([])
                if b==0 and row==0: ax.set_title(DATES[k],fontsize=8)
                if k==0: ax.set_ylabel(band.upper()+('\nSCI' if row==0 else '\nSCI – median'),fontsize=9)
    pd.DataFrame(all_rows).to_csv(out/'aperture_sensitivity.csv',index=False,float_format='%.8g')
    fig.suptitle('MSAID 25213 | 1.5 arcsec stamps, 0.18 arcsec apertures',fontsize=12)
    fig.tight_layout();fig.savefig(out/'stamps.png',dpi=120);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache200',type=Path);p.add_argument('cache444',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();run([('f200w',a.cache200),('f444w',a.cache444)],a.cache200/'targets.csv',a.out)

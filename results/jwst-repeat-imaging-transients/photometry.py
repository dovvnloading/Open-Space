"""Independent fixed-aperture light curves and forced corroboration in paired bands."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as nd, signal

EPOCH_DATE={1:'2025-06-01',2:'2025-07-30',3:'2025-09-28',4:'2025-11-28',5:'2026-01-28',
            6:'2026-03-28',7:'2026-05-28',8:'2026-08-07',9:'2026-09-28'}

def aperture(stamp,err,fx,fy,unit,radius):
    cy,cx=np.array(stamp.shape)//2
    y,x=np.indices(stamp.shape)
    rr=np.hypot(x-cx-fx,y-cy-fy)
    aperture=rr<=radius
    ann=(rr>=radius*2)&(rr<radius*3.5)
    ok=np.isfinite(stamp)&np.isfinite(err)&(err>0)
    if not ok[aperture].all() or ok[ann].sum()<.75*ann.sum():
        return np.nan,np.nan,np.nan
    sky=np.median(stamp[ann&ok])
    flux=((stamp[aperture]-sky).sum())*unit
    formal=np.sqrt((err[aperture]**2).sum())*unit
    # Identical apertures placed on a local ring are a check on both correlated
    # noise and structured host background. No independence assumption is needed.
    comparison=[]
    for theta in np.linspace(0,2*np.pi,16,endpoint=False):
        off=4*radius
        r=np.hypot(x-cx-fx-off*np.cos(theta),y-cy-fy-off*np.sin(theta))
        ap=r<=radius
        if ok[ap].all(): comparison.append((np.mean(stamp[ap])-sky)*aperture.sum()*unit)
    mad=(1.4826*np.median(abs(np.asarray(comparison)-np.median(comparison)))) if len(comparison)>=8 else np.nan
    return flux,formal,mad


def photometry(folder,out):
    rows=[]
    for file in sorted(folder.glob('ep*_data.npz')):
        stem=file.name.split('_')
        ep=int(stem[0][2:]);band=stem[1]
        errfile=file.with_name(file.name.replace('_data','_error'))
        if not errfile.exists(): continue
        with np.load(file) as a, np.load(errfile) as err:
            for i,name in enumerate(a['id']):
                flux,formal,emp=aperture(a['stamps'][i],err['stamps'][i],a['fx'][i],a['fy'][i],a['njy_per_unit'],
                                         4. if band=='f200w' else 3.)
                rows.append(dict(id=name,epoch=ep,nominal_date=EPOCH_DATE[ep],band=band,
                    host_plus_source_aperture_njy=flux,uncorr_formal_njy=formal,local_blank_scatter_njy=emp))
    pd.DataFrame(rows).to_csv(out,index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stamps',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();photometry(a.stamps,a.output)

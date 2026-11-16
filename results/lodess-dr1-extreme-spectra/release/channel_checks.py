"""Inspect J112022 and nearby controls in the six LoLSS channel catalogues.

Requires the six channel source-list FITS files named chan0.fits ... chan5.fits
from the LoLSS DR1 v1.1 public release. Acquires DR1 source list if absent.
Usage: python channel_checks.py DATA_DIRECTORY CHANNEL_DIRECTORY
"""
import argparse
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
from acquire_lodess import retrieve

def check(root, channels):
    out=root/'j112022'
    out.mkdir(parents=True,exist_ok=True)
    base='https://lofar-surveys.org/public/lolss/catalogues/LoLSS_DR1_v1.1'
    channel_files=[channels/f'chan{k}.fits' for k in range(6)]
    channels.mkdir(parents=True,exist_ok=True)
    for k,f in enumerate(channel_files):
        if not f.exists():retrieve(base+f'-chan{k:02d}.srl.fits',f)
    full=channels/'lolss.fits'
    if not full.exists():retrieve(base+'.srl.fits',full)
    s=fits.getdata(full)
    refs=SkyCoord(s['RA']*u.deg,s['DEC']*u.deg)
    target=SkyCoord(170.092166*u.deg,50.493238*u.deg)
    sep=target.separation(refs).deg
    _,near,_=refs.match_to_catalog_sky(refs,nthneighbor=2)
    ix=np.where((sep<2)&(near.arcsec>60)&(s['S_Code']=='S')&(s['Peak_flux']/s['Isl_rms']>30)&
                (s['DC_Maj']<12)&(s['Total_flux']>0))[0]
    freq=np.arange(44,65,4)
    matrices={k:np.empty((len(ix),6)) for k in ['total','peak','e_total','e_peak','maj','min','ra','dec']}
    exists=np.ones(len(ix),dtype=bool)
    scode=[]
    for k,file in enumerate(channel_files):
        a=fits.getdata(file)
        nn,d,_=refs[ix].match_to_catalog_sky(SkyCoord(a['RA']*u.deg,a['DEC']*u.deg))
        pick=a[nn]
        exists &= (d.arcsec<5)&(pick['Peak_flux']/pick['Isl_rms']>10)
        for name,col in [('total','Total_flux'),('peak','Peak_flux'),('e_total','E_Total_flux'),('e_peak','E_Peak_flux'),('maj','Maj'),('min','Min'),('ra','RA'),('dec','DEC')]:
            matrices[name][:,k]=pick[col]
        scode.append(pick['S_Code'])
    it=np.argmin(sep[ix])
    if not exists[it] or sep[ix[it]]>1/3600:raise ValueError('Target did not meet channel quality verification')
    d=pd.DataFrame(dict(mhz=freq,**{k:v[it] for k,v in matrices.items()},s_code=[a[it] for a in scode]))
    d.to_csv(out/'channel_target.csv',index=False)
    X=np.column_stack([np.ones(6),np.log(freq/54)])
    reg=np.linalg.pinv(X)[1]
    data={'lolss_row':ix,'ra_deg':s['RA'][ix].astype(float),'dec_deg':s['DEC'][ix].astype(float),'separation_deg':sep[ix]}
    for f in ['total','peak']:
        fl=matrices[f]
        data['alpha_'+f]=np.log(fl)@reg
    data['delta_alpha']=data['alpha_peak']-data['alpha_total']
    df=pd.DataFrame(data)[exists]
    df.to_csv(out/'channel_local_slopes.csv',index=False)
    targetdata=df[df.lolss_row==ix[it]].iloc[0]
    ctrl=df[df.lolss_row!=ix[it]]
    print('log-linear unweighted slopes',targetdata.to_dict())
    print('controls',len(ctrl))
    print(ctrl[['alpha_total','alpha_peak','delta_alpha']].quantile([.16,.5,.84]).to_string())
    print('total slope above/equal target',int((ctrl.alpha_total>=targetdata.alpha_total).sum()))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('channels',type=Path)
    args=p.parse_args();check(args.root,args.channels)

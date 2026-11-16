"""Conditional steady-coma spatial test: PSF-convolved rho^-p, not chemistry."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import log
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
from scipy.optimize import lsq_linear
log.setLevel('ERROR')
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--cache',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=[]
for row in json.load(open(a.input)):
 if not(row['start'].startswith('2025-12-25') or row['start'].startswith('2026-01-02 10:')):continue
 path=Path(a.cache)/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')
 with fits.open(path) as h:
  im=h['IMAGE'].data;var=h['VARIANCE'].data;fl=h['FLAGS'].data;x,y=row['x'],row['y'];yy,xx=np.indices(im.shape);r=np.hypot(xx-x,yy-y)
  good=np.isfinite(im)&np.isfinite(var)&(var>0)&((fl&~(2097152|524288))==0)
  if ((r<=2.5)&~good).any():continue
  w=WCS(h['IMAGE'].header,key='A');w.sip=None;ax,ay=w.wcs_pix2world(x,y,0);tb=h['EPSF'].data;k=np.argmin((tb['XCENTER']-ax)**2+(tb['YCENTER']-ay)**2);ps=tb['EPSF'][k].astype(float);ps/=ps.sum()
  # 40 native pixels radius, 5x sampling. Core radius is regularization, not nucleus size.
  gy,gx=np.indices((401,401));rr=np.hypot(gx-200,gy-200)/5
  point=25*map_coordinates(ps,[(yy-y)*5+16,(xx-x)*5+16],order=3,mode='constant')
  entry=dict(start=row['start'],wavelength=row['wavelength_um'],url=row['url'],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),fits=[])
  for radius in [3.5,5.,8.]:
   mask=good&(r<=radius);z=im[mask];err=np.sqrt(var[mask]);pt=point[mask]
   def fit(columns):
    A=np.array(columns).T;fit=lsq_linear(A/err[:,None],z/err,bounds=([0]*(len(columns)-1)+[-np.inf],[np.inf]*len(columns)),tol=1e-12)
    return dict(parameters=fit.x.tolist(),chi2=float(np.sum(((A@fit.x-z)/err)**2)),n=int(len(z)),parameter_count=len(columns))
   base=fit([pt,np.ones(len(z))])
   for power in [.5,1.,1.5]:
    for core in [.1,.5]:
     coma=(rr**2+core**2)**(-power/2);conv=fftconvolve(coma,ps,mode='same');model=map_coordinates(conv,[(yy-y)*5+200,(xx-x)*5+200],order=3,mode='constant');model/=model[r<=5].sum()
     cm=model[mask];pure=fit([cm,np.ones(len(z))]);mixed=fit([pt,cm,np.ones(len(z))])
     entry['fits'].append(dict(radius=radius,power=power,core_pixels=core,point=base,coma=pure,mixed=mixed))
  out.append(entry)
Path(a.output).write_text(json.dumps(dict(caveats=['Fixed ephemeris center; no centroid optimization','rho^-1 represents optically thin stationary isotropic constant-speed coma without significant destruction inside aperture; dust can also follow it','Free power/core are sensitivity grid, not measured physical parameters','OUTLIER included; other bad pixels excluded; strict core fails for these signals','Formal independent-pixel noise only; no significance or gas identification','EPSF finite wings and neighboring sources unresolved; background constant jointly fitted','Coma normalization is unit sum inside native radius5; nucleus coefficient uses truncated PSF normalization'],results=out),indent=2)+'\n')
for e in out:
 best=min([f for f in e['fits'] if f['radius']==5],key=lambda f:f['mixed']['chi2']);print(e['start'],e['wavelength'],'point',best['point']['chi2'],'coma',best['coma']['chi2'],'mixed',best['mixed']['chi2'],'p',best['power'],'core',best['core_pixels'],'coef',best['mixed']['parameters'])

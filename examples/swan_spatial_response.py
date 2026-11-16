"""Separable light-profile weighted spectral response, not full radiometric retrieval."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve
from scipy.integrate import trapezoid
from scipy.optimize import lsq_linear
from astropy.io import fits
from astropy.wcs import WCS
from astropy import log
log.setLevel('ERROR')
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--response-input',required=True);p.add_argument('--calibration-root',required=True);p.add_argument('--cache',required=True);p.add_argument('--spectrum',required=True);p.add_argument('--output',required=True);a=p.parse_args();cal=Path(a.calibration_root)
bc=fits.getdata(cal/'ParameterFiles/20260604_SSDC_BC_Band5.fits',1);ch=fits.getdata(cal/'ChannelMaps/spectral_channels_17ch_2040.fits',1);su=np.genfromtxt(next((cal/'FidSpecSummary/SubChannel17/Band5').glob('*.csv')),delimiter=',',names=True);spectra={}
for k in range(1,18):spectra[k]=np.genfromtxt(cal/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{k:02d}_17_Cent.csv',delimiter=',',names=True)
s=np.loadtxt(a.spectrum);wave=s[:,0];mol=s[:,-1];mol/=trapezoid(mol,wave)
responses={r['observationid']:r for r in json.load(open(a.response_input))};out=[]
for r in json.load(open(a.manifest)):
 if r.get('error') or r['start'][:10] not in ['2025-12-22','2025-12-25','2026-01-02'] or r['observationid'] not in responses:continue
 z=responses[r['observationid']]
 if not z['diagnostic_outlier_inclusive'].get('2.5') or z['map_bc_um']<4.1:continue
 with fits.open(Path(a.cache)/(hashlib.sha256(r['url'].encode()).hexdigest()+'.fits')) as h:
  yy,xx=np.indices(h['IMAGE'].data.shape);rad=np.hypot(xx-r['x'],yy-r['y']);m=rad<=12;w=WCS(h['IMAGE'].header,key='A');w.sip=None;ax,ay=w.wcs_pix2world(xx[m],yy[m],0);bcs=map_coordinates(bc,[ay,2039-ax],order=1);chs=ch[np.rint(ay).astype(int),np.rint(ax).astype(int)];q=[]
  def response(c,k):
   sp=spectra[int(k)];wv=sp['Wavelength_um']+c-su['BandCenter_um'][int(k)-1];R=sp['NormTransmission'];return trapezoid(np.interp(wave,wv,R,left=0,right=0)*mol*wave,wave)/trapezoid(R/wv,wv)
  for b,k in zip(bcs,chs):q.append(response(b,k))
  q=np.array(q);tx,ty=w.wcs_pix2world(r['x'],r['y'],0);tb=h['EPSF'].data;k=np.argmin((tb['XCENTER']-tx)**2+(tb['YCENTER']-ty)**2);ps=tb['EPSF'][k].astype(float);ps/=ps.sum();gy,gx=np.indices((401,401));rr=np.hypot(gx-200,gy-200)/5;weights={}
  for power in [0.,1.,1.5]:
   if power==0:profile=25*map_coordinates(ps,[(yy-r['y'])*5+16,(xx-r['x'])*5+16],order=3,mode='constant')
   else:profile=map_coordinates(fftconvolve((rr**2+.1**2)**(-power/2),ps,mode='same'),[(yy-r['y'])*5+200,(xx-r['x'])*5+200],order=3,mode='constant')
   for ap in [2.5,3.5]:
    inside=rad[m]<=ap;sky=(rad[m]>=8)&(rad[m]<=12);wt=profile[m]*inside
    for bg in [False,True]:
     weight=wt-profile[m]*sky*inside.sum()/sky.sum() if bg else wt
     weights[f'p{power}_r{ap}_meansky{bg}']=float(np.sum(weight*q)/np.sum(weight))
  out.append(dict(start=r['start'],id=r['observationid'],center_response=response(z['map_bc_um'],z['channel']),weighted_responses=weights,photometry=z['diagnostic_outlier_inclusive']))
fitsout=[]
for key in out[0]['weighted_responses']:
 ap='3.5' if '_r3.5_' in key else '2.5';rows=[]
 for epoch in ['2025-12-22','2025-12-25','2026-01-02']:
  es=[e for e in out if e['start'].startswith(epoch)];X=np.array([[1,e['weighted_responses'][key]] for e in es]);y=np.array([e['photometry'][ap]['sum_mjysr'] for e in es]);err=np.array([e['photometry'][ap]['error_mjysr'] for e in es]);beta=lsq_linear(X/err[:,None],y/err,bounds=([-np.inf,0],[np.inf,np.inf])).x;rows.append(dict(epoch=epoch,chi2=float(np.sum(((X@beta-y)/err)**2)),dof=len(y)-2))
 fitsout.append(dict(model=key,epochs=rows))
Path(a.output).write_text(json.dumps(dict(caveats=['100K CO2 fixed from earlier December selection','Assumes spatially invariant spectrum and separable observed PSF/profile; not response-before-optical-convolution simulation','Mean-sky transfer approximates actual nonlinear median sky','Negative sky weights included only in meanskyTrue; truncated ePSF wings unresolved','Formal independent noise only; no gas identification'],observations=out,fits=fitsout),indent=2)+'\n')
for f in fitsout:print(f['model'],[round(e['chi2'],3) for e in f['epochs']])

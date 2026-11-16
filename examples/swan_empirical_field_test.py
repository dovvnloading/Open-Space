"""Translated-track field controls with local calibrated response; diagnostic only."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.integrate import trapezoid
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u,log
from spherex_moving_screen import aperture
log.setLevel('ERROR')
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--training',required=True);p.add_argument('--calibration',required=True);p.add_argument('--cache',required=True);p.add_argument('--output',required=True);a=p.parse_args()
rows=json.load(open(a.input));train=json.load(open(a.training))['fits'][0];mu=train['center_um'];sigma=train['sigma_um'];cal=Path(a.calibration)
bc=fits.getdata(cal/'ParameterFiles/20260604_SSDC_BC_Band5.fits',1);ch=fits.getdata(cal/'ChannelMaps/spectral_channels_17ch_2040.fits',1)
summ=np.genfromtxt(next((cal/'FidSpecSummary/SubChannel17/Band5').glob('*.csv')),delimiter=',',names=True)
curves={n:np.genfromtxt(cal/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{n:02d}_17_Cent.csv',delimiter=',',names=True) for n in range(1,18)}
def response(h,x,y):
 w=WCS(h['IMAGE'].header,key='A');w.sip=None;ax,ay=w.wcs_pix2world(x,y,0);center=float(map_coordinates(bc,[[ay],[2039-ax]],order=1)[0]);c=int(ch[round(float(ay)),round(float(ax))]);s=curves[c];wave=s['Wavelength_um']+center-summ['BandCenter_um'][c-1];r=s['NormTransmission'];g=np.exp(-.5*((wave-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi));t=trapezoid(r*g*wave,wave)/trapezoid(r/wave,wave);return center,float(t)
def position(seq,offset):
 ys=[];es=[];bs=[];ts=[]
 for row,h in seq:
  coord=SkyCoord(row['ra']*u.deg,row['dec']*u.deg).spherical_offsets_by(offset[0]*u.arcsec,offset[1]*u.arcsec)
  x,y=WCS(h['IMAGE'].header,h).all_world2pix([[coord.ra.deg,coord.dec.deg]],0)[0]
  m=aperture(h['IMAGE'].data,h['VARIANCE'].data,h['FLAGS'].data,x,y,2.5,True)
  if m is None:return {'offset_arcsec':offset,'valid':False}
  b,t=response(h,x,y);ys.append(m['sum_mjysr']);es.append(m['error_mjysr']);bs.append(b);ts.append(t)
 y=np.array(ys);e=np.array(es);b=np.array(bs);t=np.array(ts);X=np.column_stack([np.ones(len(y)),b-4.26,t]);normal=X.T@(X/e[:,None]**2)
 if np.linalg.cond(normal)>1e12:return {'offset_arcsec':offset,'valid':False,'reason':'ill_conditioned_response'}
 cov=np.linalg.inv(normal);beta=cov@X.T@(y/e**2);pred=X@beta;se=np.sqrt(cov[2,2]);X0=X[:,:2];beta0=np.linalg.lstsq(X0/e[:,None],y/e,rcond=None)[0]
 return dict(offset_arcsec=offset,valid=True,flux=ys,error=es,centers=bs,response=ts,feature_amplitude=float(beta[2]),formal_amplitude_error=float(se),studentized_amplitude=float(beta[2]/se),continuum_chi2=float(np.sum(((y-X0@beta0)/e)**2)),feature_chi2=float(np.sum(((y-pred)/e)**2)),dof=len(y)-3)
report={'warning':'Local translated-track controls overlap and are not independent trials. Finite ranks are descriptive, not discovery p-values. Same OUTLIER-inclusive extraction for target and controls; no catalog-based postselection. Response shape frozen from December training.', 'frozen_center_um':mu,'frozen_sigma_um':sigma,'epochs':[]}
for date,lo,hi in [('2026-01-02','10:18','10:25'),('2026-09-11','17:39','21:02')]:
 selected=[r for r in rows if r['start'][:10]==date and lo<=r['start'][11:16]<=hi and 'error' not in r];seq=[(r,fits.open(Path(a.cache)/(hashlib.sha256(r['url'].encode()).hexdigest()+'.fits'))) for r in selected]
 target=position(seq,[0,0]);controls=[position(seq,[dx,dy]) for dx in range(-120,121,20) for dy in range(-120,121,20) if 60<=np.hypot(dx,dy)<=140]
 valid=[c for c in controls if c['valid']];v=np.array([c['studentized_amplitude'] for c in valid]);median=float(np.median(v));mad=float(1.4826*np.median(abs(v-median)))
 z=target.get('studentized_amplitude');summary=dict(n_valid=len(valid),median_studentized=median,robust_scale_studentized=mad,largest_studentized=float(max(v)),target_above_count=int(np.sum(v>=z)) if z is not None else None,descriptive_upper_rank=float((1+np.sum(v>=z))/(1+len(v))) if z is not None else None)
 report['epochs'].append(dict(date=date,observations=[r['observationid'] for r in selected],target=target,controls=controls,summary=summary))
 for _,h in seq:h.close()
Path(a.output).write_text(json.dumps(report,indent=2)+'\n')
for ep in report['epochs']:print(ep['date'],ep['target'],ep['summary'])

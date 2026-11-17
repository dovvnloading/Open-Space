"""Frozen development-only mask, S/N and narrow-line control; never reads holdout."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
src=Path(sys.argv[1]);rec=json.loads(Path(sys.argv[2]).read_text());out=Path(sys.argv[3]);records=[]
for r in rec:
 if r['kind']=='hold':continue
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p,memmap=False) as h:
  co=h['COADD'].data;z=float(next(x['z'] for x in rec if x['sdss_id']==r['sdss_id'] and x['kind']=='baseline'))
  wave=10**co['loglam']/(1+z);flux=np.asarray(co['flux'],float);iv=np.asarray(co['ivar'],float);mask=co['and_mask'];valid=(iv>0)&np.isfinite(flux)&np.isfinite(iv)&(mask==0)
  def sel(a,b):return (wave>=a)&(wave<b)
  def frac(a,b):
   w=sel(a,b);return float(np.sum(valid&w)/np.sum(w)) if np.sum(w) else 0
  cont=(sel(4750,4770)|sel(5080,5100))&valid;sn=float(np.median(flux[cont]*np.sqrt(iv[cont]))) if cont.any() else float('nan')
  a=sel(4975,4988)&valid;b=sel(5025,5038)&valid;core=sel(4998,5018)&valid
  if min(frac(4975,4988),frac(5025,5038),frac(4998,5018))>=.8 and a.any() and b.any() and core.any():
   xa=np.median(wave[a]);xb=np.median(wave[b]);fa=np.median(flux[a]);fb=np.median(flux[b]);x=wave[core];fw=flux[core];dx=np.gradient(wave)[core]; continuum=fa+(fb-fa)*(x-xa)/(xb-xa)
   line=float(np.sum((fw-continuum)*dx));err=float(np.sqrt(np.sum(dx**2/iv[core])));n=int(core.sum())
  else:line=err=float('nan');n=0
  records.append(dict(sdss_id=r['sdss_id'],kind=r['kind'],domain_coverage=round(frac(4750,5100),6),continuum_snr=round(sn,5),oiii_flux=round(line,5),oiii_formal_core_err=round(err,5),oiii_pixels=n))
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
print('development spectra',len(records),'passing coverage and SNR',sum(x['domain_coverage']>=.8 and x['continuum_snr']>=5 for x in records))

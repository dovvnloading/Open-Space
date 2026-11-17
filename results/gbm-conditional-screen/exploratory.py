"""Explicitly authorized descriptive joint characterization; no p-values or rankings."""
import argparse,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from manifest import RAW,sha,verify
from fits_reader import metadata
from observational import SPECS

def authorize(path,digest,manifest_digest):
 if sha(path)!=digest:raise ValueError('authorization digest')
 a=json.loads(Path(path).read_text())
 if a.get('authorized') is not True or a.get('phase')!='exploratory_characterization' or a.get('manifest_sha256')!=manifest_digest:raise ValueError('authorization content')
def describe(times,pha,selected,start,n):
 end=start+4*n;mask=(times>=start)&(times<end);times=times[mask];pha=pha[mask];edges=start+np.arange(n*256+1)/64;total=np.histogram(times,bins=edges)[0];overflow=np.histogram(times[pha==127],bins=edges)[0];dead=(total-overflow)*2.6e-6+overflow*1e-5
 if np.any(dead>=1/64):raise ValueError('deadtime exceeds bin duration')
 # Full-band diagnostics are computed before selecting any energy channels.
 histogram=np.bincount(pha.astype(int),minlength=128);band=np.isin(pha,selected);light=np.histogram(times[band],bins=edges)[0];windows=[];power=[]
 for w in range(n):
  q=light[w*256:(w+1)*256];norm=int(q.sum());x=q-q.mean();p=np.zeros(129) if norm==0 else 2*np.abs(np.fft.rfft(x))**2/norm;p[-1]*=.5;power.append(p[20:129].tolist());sl=slice(w*256,(w+1)*256);windows.append({'fullband_count':int(total[sl].sum()),'overflow_count':int(overflow[sl].sum()),'nominal_deadtime_seconds':float(dead[sl].sum()),'nominal_livetime_seconds':float(4-dead[sl].sum()),'band_count':norm})
 return {'fullband_count':int(total.sum()),'overflow_count':int(overflow.sum()),'pha_histogram':histogram.tolist(),'fullband_lightcurve':total.tolist(),'band_lightcurve':light.tolist(),'windows':windows,'fourier_power':power,'limitations':['nominal event-deadtime accounting only; boundary crossing, unrecorded events, pileup, saturation and housekeeping remain unverified']}
def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();m=verify(a.manifest,a.manifest_sha256);authorize(a.authorization,a.authorization_sha256,a.manifest_sha256)
 if a.output.exists() or a.output.resolve().drive.upper()!='D:':raise ValueError('output collision or not D drive')
 handles=[];out={}
 try:
  for bid,det,name,trig,start,n in SPECS:
   filename=f'glg_tte_{det}_bn{bid}_v00.fit'
   if filename not in m['tte_sha256']:raise ValueError('input not in manifest')
   h=fits.open(RAW/filename,memmap=True);handles.append(h);metadata(h,name,bid,trig,start,start+n*4)
  # All four headers, GTIs and EBOUNDS are verified before any EVENTS decoding.
  for h,(bid,det,name,trig,start,n) in zip(handles,SPECS):
   ch,selected=metadata(h,name,bid,trig,start,start+n*4)
   if not np.array_equal(ch,np.arange(128)):raise ValueError('overflow channel mapping unverified')
   ev=h['EVENTS'].data;times=np.asarray(ev['TIME']);pha=np.asarray(ev['PHA']);
   if not np.isfinite(times).all() or np.any(np.diff(times)<0) or not np.isin(pha,ch).all():raise ValueError('invalid photon rows')
   out.setdefault(bid,{})[det]=describe(times-trig,pha,selected,start,n)
  for bid,da,db in [('261002867','n5','n2'),('261007236','nb','n8')]:
   out[bid]['joint_power']=np.minimum(out[bid][da]['fourier_power'],out[bid][db]['fourier_power']).tolist()
  with a.output.open('x',encoding='utf-8') as f:json.dump({'manifest_sha256':a.manifest_sha256,'phase':'exploratory_characterization','dt_seconds':1/64,'frequencies_hz':np.fft.rfftfreq(256,1/64)[20:129].tolist(),'bursts':out,'scope':'descriptive powers; no p-values, ranking, discovery or independent-confirmation claim'},f,allow_nan=False);f.write('\n')
 finally:
  for h in handles:h.close()
if __name__=='__main__':main()

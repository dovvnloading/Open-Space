"""Fit OGLE-predicted dip amplitudes to the Gaia DR3 BP/RP transit fluxes."""
import csv,json
import erfa
from common import D,np,SKY,shape,PARAMS
rows=list(csv.DictReader((D/'inputs/gaia-epochs.csv').open()))
line=np.array([np.cos(SKY.ra.radian)*np.cos(SKY.dec.radian),np.sin(SKY.ra.radian)*np.cos(SKY.dec.radian),np.sin(SKY.dec.radian)])
def hjd_utc(t):
 a,b=erfa.tcbtdb(2451545.,t+2455197.5-2451545.)
 hel,bar=erfa.epv00(a,b);proj=(bar['p']-hel['p'])@line*499.004783836/86400
 a,b=erfa.tttai(a,b);a,b=erfa.taiutc(a,b)
 return a+b-proj
def run(b):
 tname=b+'_obs_time';fn=b+'_flux';en=b+'_flux_error';m=b+'_mag'
 clean=[r for r in rows if r[tname] and r[fn] and r[en] and r[m] and r['variability_flag_'+b+'_reject']=='false' and r['rejected_by_photometry']=='false' and int(r[b+'_other_flags'])==0 and float(r[en])>0]
 t=hjd_utc(np.array([float(r[tname]) for r in clean]));f=np.array([float(r[fn]) for r in clean]);se=np.array([float(r[en]) for r in clean])
 phi,d1,d2=shape(t);scale=np.median(f[(abs(d1)>.09)&(abs(d2)>.09)])
 f=f/scale;e=np.hypot(se/scale,.05)
 X=np.column_stack((np.ones(len(t)),-np.exp(-.5*(d1/PARAMS[4])**2),-np.exp(-.5*(d2/PARAMS[5])**2)))
 p=np.linalg.lstsq(X/e[:,None],f/e,rcond=None)[0]
 nights=np.floor(t-.5);S=X/e[:,None];res=(f-X@p)/e;Ai=np.linalg.pinv(S.T@S)
 sc=np.array([S[nights==d].T@res[nights==d] for d in np.unique(nights)]);cov=Ai@sc.T@sc@Ai*len(t)/(len(t)-3)
 return dict(count=len(f),hard_half_depth=[int((X[:,i]<-.5).sum()) for i in [1,2]],depths=p[1:].tolist(),cluster_se=np.sqrt(np.diag(cov))[1:].tolist(),off_eclipse_rate_e_s=scale,flags='variability and photometry accepted, per-band other_flags=0')
out={b:run(b) for b in ['bp','rp']};(D/'gaia-colors.json').write_text(json.dumps(out,indent=2));print(out)

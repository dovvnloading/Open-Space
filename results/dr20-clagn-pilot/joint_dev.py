"""Frozen joint development fits and deterministic real-spectrum injections; no held-out flux."""
import csv,hashlib,json,sys,collections
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import lsq_linear
src=Path(sys.argv[1]); rec=json.loads(Path(sys.argv[2]).read_text());qa=list(csv.DictReader(open(sys.argv[3],newline='')));q=collections.defaultdict(dict);pairs=collections.defaultdict(dict)
for x in qa:q[int(x['sdss_id'])][x['kind']]=x
for x in rec:
 if x['kind']!='hold':pairs[x['sdss_id']][x['kind']]=x

def gaussian(x,mu,sig):return np.exp(-.5*((x-mu)/sig)**2)/(np.sqrt(2*np.pi)*sig)
def load(r,z):
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p,memmap=False) as h:
  a=h['COADD'].data;w=10**a['loglam']/(1+z);fl=np.asarray(a['flux'],float);iv=np.asarray(a['ivar'],float);good=(w>=4750)&(w<5100)&(iv>0)&np.isfinite(fl)&np.isfinite(iv)&(a['and_mask']==0);return w[good],fl[good],iv[good]
def fit(w,y,iv,sigma=20,shift=0,degree=2):
 x=(w-4925)/175;cont=[np.ones(len(w)),x]
 if degree==2:cont.append(x*x)
 absorption=gaussian(w,4862.68+shift,7); narrow=gaussian(w,4862.68+shift,2.5); oiii=(gaussian(w,5008.239+shift,2.5)+gaussian(w,4960.295+shift,2.5)/2.98)/(1+1/2.98);broad=gaussian(w,4862.68+shift,sigma)
 A=np.column_stack(cont+[absorption,narrow,oiii,broad]); weight=np.sqrt(iv)/1.2;Aw=A*weight[:,None];yw=y*weight;nf=len(cont);lower=np.array([-np.inf]*nf+[-np.inf,0,0,0]);upper=np.array([np.inf]*nf+[0,np.inf,np.inf,np.inf]);res=lsq_linear(Aw,yw,bounds=(lower,upper),tol=1e-10,max_iter=300)
 cov=np.linalg.pinv(Aw.T@Aw);n,o,b=nf+1,nf+2,nf+3;of,bf=res.x[o],res.x[b];oe,be=np.sqrt(max(0,cov[o,o])),np.sqrt(max(0,cov[b,b]));ratio=bf/of if of>5*oe else np.nan;g=np.zeros(len(res.x));
 if of>5*oe:g[b]=1/of;g[o]=-bf/of**2
 err=np.sqrt(max(0,g@cov@g)) if of>5*oe else np.nan
 return dict(success=res.success,broad=float(bf),broad_err=float(be),oiii=float(of),oiii_err=float(oe),ratio=float(ratio),ratio_err=float(err),chi2=float(np.sum((Aw@res.x-yw)**2)),dof=len(y)-len(res.x))
rows=[];inject=[];spectra={}
for sid,pair in sorted(pairs.items()):
 if not all(float(q[sid][k]['domain_coverage'])>=.8 and float(q[sid][k]['continuum_snr'])>=5 for k in ('baseline','dev')):continue
 data={k:load(pair[k],pair['baseline']['z']) for k in ('baseline','dev')};spectra[sid]=data
 for sigma in (20,35):
  for shift in (-5,0,5):
   for degree in (1,2):
    a=fit(*data['baseline'],sigma,shift,degree);b=fit(*data['dev'],sigma,shift,degree)
    diff=b['ratio']-a['ratio'];err=np.hypot(a['ratio_err'],b['ratio_err']); valid=all(x['success'] and x['oiii']>5*x['oiii_err'] for x in [a,b]);high=max(a['ratio'],b['ratio']);low=min(a['ratio'],b['ratio']);decision=bool(valid and low>0 and high/low>=5 and abs(diff)>5*err and ((a['broad']>5*a['broad_err']) if a['ratio']>=b['ratio'] else (b['broad']>5*b['broad_err'])))
    rows.append(dict(sdss_id=sid,sigma=sigma,shift=shift,degree=degree,valid=int(valid),baseline_broad=a['broad'],baseline_broad_err=a['broad_err'],baseline_oiii=a['oiii'],baseline_oiii_err=a['oiii_err'],baseline_ratio=a['ratio'],baseline_ratio_err=a['ratio_err'],baseline_chi2=a['chi2'],baseline_dof=a['dof'],dev_broad=b['broad'],dev_broad_err=b['broad_err'],dev_oiii=b['oiii'],dev_oiii_err=b['oiii_err'],dev_ratio=b['ratio'],dev_ratio_err=b['ratio_err'],dev_chi2=b['chi2'],dev_dof=b['dof'],conditional_threshold=int(decision)))
 a=fit(*data['baseline']);b0=fit(*data['dev']);w,y,iv=data['dev']
 for factor in (0,.5,2):
  added=factor*a['broad'];b=fit(w,y+added*gaussian(w,4862.68,20),iv);inject.append(dict(sdss_id=sid,factor=factor,baseline_broad=a['broad'],injected=added,original_dev=b0['broad'],recovered_dev=b['broad'],response=b['broad']-b0['broad']))
def write(path,rs):
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rs[0].keys());w.writeheader();w.writerows(rs)
write(sys.argv[4],rows);write(sys.argv[5],inject)
primary=[r for r in rows if r['sigma']==20 and r['shift']==0 and r['degree']==2];print('pairs',len(primary),'oiii qualified',sum(r['valid'] for r in primary),'conditional triggers',sum(r['conditional_threshold'] for r in primary))

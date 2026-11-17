"""Post-inspection shared narrow-profile sensitivity, development spectra only."""
import csv,hashlib,json,sys,collections
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import lsq_linear
src=Path(sys.argv[1]);rec=json.loads(Path(sys.argv[2]).read_text());qa=list(csv.DictReader(open(sys.argv[3],newline='')));q=collections.defaultdict(dict);pairs=collections.defaultdict(dict)
for r in qa:q[int(r['sdss_id'])][r['kind']]=r
for r in rec:
 if r['kind']!='hold':pairs[r['sdss_id']][r['kind']]=r

def g(w,mu,s):return np.exp(-.5*((w-mu)/s)**2)/(np.sqrt(2*np.pi)*s)
def data(r,z):
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p) as h:
  c=h['COADD'].data;w=10**c['loglam']/(1+z);y=np.array(c['flux'],float);iv=np.array(c['ivar'],float);m=(w>=4750)&(w<5100)&(iv>0)&np.isfinite(y)&(c['and_mask']==0);return w[m],y[m],iv[m]
def fit(w,y,iv,shift,width,broad=20,degree=2):
 x=(w-4925)/175; cont=[np.ones(len(w)),x]+([x*x] if degree==2 else []);n=len(cont)
 def doublet(s):return (g(w,5008.239+shift,s)+g(w,4960.295+shift,s)/2.98)/(1+1/2.98)
 A=np.column_stack(cont+[g(w,4862.68+shift,7),g(w,4862.68+shift,width),doublet(width),doublet(2*width),g(w,4862.68+shift,broad)]);weight=np.sqrt(iv)/1.2;Aw=A*weight[:,None];yw=y*weight;lo=[-np.inf]*n+[-np.inf,0,0,0,0];hi=[np.inf]*n+[0,np.inf,np.inf,np.inf,np.inf];r=lsq_linear(Aw,yw,bounds=(lo,hi),tol=1e-9,max_iter=300);cov=np.linalg.pinv(Aw.T@Aw);oi=r.x[n+2]+r.x[n+3];bi=r.x[n+4];ov=cov[n+2,n+2]+cov[n+3,n+3]+2*cov[n+2,n+3];oe=np.sqrt(max(0,ov));be=np.sqrt(max(0,cov[n+4,n+4]));ratio=bi/oi if oi>5*oe else np.nan;v=np.zeros(len(r.x));
 if oi>5*oe:v[n+4]=1/oi;v[n+2]=v[n+3]=-bi/(oi*oi)
 re=np.sqrt(max(0,v@cov@v)) if oi>5*oe else np.nan
 return dict(success=r.success,chi2=float(np.sum((Aw@r.x-yw)**2)),dof=len(y)-len(r.x),broad=float(bi),broad_err=float(be),oiii=float(oi),oiii_err=float(oe),ratio=float(ratio),ratio_err=float(re))
rows=[]
for sid,pair in sorted(pairs.items()):
 if not all(float(q[sid][k]['domain_coverage'])>=.8 and float(q[sid][k]['continuum_snr'])>=5 for k in ('baseline','dev')):continue
 z=pair['baseline']['z'];d={k:data(pair[k],z) for k in ('baseline','dev')};grid=[]
 for shift in (-4,-2,0,2,4):
  for width in (2.5,4,6):
   a=fit(*d['baseline'],shift,width);b=fit(*d['dev'],shift,width);grid.append((a['chi2']+b['chi2'],shift,width))
 _,shift,width=min(grid)
 for broad in (20,35):
  for degree in (1,2):
   a=fit(*d['baseline'],shift,width,broad,degree);b=fit(*d['dev'],shift,width,broad,degree);row=dict(sdss_id=sid,shift=shift,narrow_sigma=width,broad_sigma=broad,degree=degree)
   for k,v in a.items():row['baseline_'+k]=v
   for k,v in b.items():row['dev_'+k]=v
   rows.append(row)
with open(sys.argv[4],'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('pairs',len(rows)//4)

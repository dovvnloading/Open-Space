"""Post-inspection RM17 residual-block sensitivity; development data only."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import lsq_linear
src=Path(sys.argv[1]);records=json.loads(Path(sys.argv[2]).read_text());items=[x for x in records if x['sdss_id']==62239511 and x['kind']!='hold'];z=next(x['z'] for x in items if x['kind']=='baseline');rng=np.random.default_rng(20261117)
def gauss(x,m,s):return np.exp(-.5*((x-m)/s)**2)/(np.sqrt(2*np.pi)*s)
def prepare(r):
 p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p) as h:
  c=h['COADD'].data;wave=10**c['loglam']/(1+z);flux=np.array(c['flux'],float);iv=np.array(c['ivar'],float);m=(wave>=4750)&(wave<5100)&(iv>0)&np.isfinite(flux)&(c['and_mask']==0);return wave[m],flux[m],iv[m]
def fit(w,y,iv):
 shift=-2;width=2.5;x=(w-4925)/175
 def db(s):return (gauss(w,5008.239+shift,s)+gauss(w,4960.295+shift,s)/2.98)/(1+1/2.98)
 A=np.column_stack([np.ones(len(w)),x,x*x,gauss(w,4862.68+shift,7),gauss(w,4862.68+shift,width),db(width),db(2*width),gauss(w,4862.68+shift,20)]);weight=np.sqrt(iv)/1.2;Aw=A*weight[:,None];r=lsq_linear(Aw,y*weight,bounds=([-np.inf]*4+[0]*4,[np.inf]*3+[0]+[np.inf]*4),tol=1e-9,max_iter=300);oi=r.x[5]+r.x[6];return (r.x[7]/oi if oi>0 else np.nan),A@r.x,r.success
series={}
for kind in ('baseline','dev'):
 w,y,iv=prepare(next(x for x in items if x['kind']==kind));ratio,model,ok=fit(w,y,iv);standard=(y-model)*np.sqrt(iv)/1.2;standard-=np.mean(standard);series[kind]=(w,y,iv,model,standard,ratio)
results=[]
for i in range(200):
 vals={};success=True
 for kind in ('baseline','dev'):
  w,y,iv,model,res,ratio=series[kind];starts=rng.integers(0,max(1,len(res)-4),size=(len(res)+4)//5);sample=np.concatenate([res[s:s+5] for s in starts])[:len(res)];v,_,ok=fit(w,model+sample*1.2/np.sqrt(iv),iv);vals[kind]=v;success &= ok and np.isfinite(v)
 if success:results.append([vals['baseline'],vals['dev'],vals['dev']-vals['baseline']])
a=np.asarray(results);o=dict(n_requested=200,n_success=len(results),baseline_original=series['baseline'][-1],dev_original=series['dev'][-1],percentiles_2p5_50_97p5=np.percentile(a,[2.5,50,97.5],axis=0).tolist());Path(sys.argv[3]).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))

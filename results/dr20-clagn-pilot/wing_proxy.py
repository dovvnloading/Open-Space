"""Exploratory fixed-wing/[O III] proxy; development data only, not state classification."""
import csv,hashlib,json,sys,collections
from pathlib import Path
import numpy as np
from astropy.io import fits
src=Path(sys.argv[1]);rec=json.loads(Path(sys.argv[2]).read_text());qa=list(csv.DictReader(open(sys.argv[3],newline='')));q=collections.defaultdict(dict)
for a in qa:q[int(a['sdss_id'])][a['kind']]=a
rby=collections.defaultdict(dict)
for r in rec:
 if r['kind']!='hold':rby[r['sdss_id']][r['kind']]=r
out=[]
for sid,pair in sorted(rby.items()):
 ok=True
 for k in ['baseline','dev']:
  a=q[sid][k];ok &= float(a['domain_coverage'])>=.8 and float(a['continuum_snr'])>=5 and float(a['oiii_flux'])>5*float(a['oiii_formal_core_err'])
 if not ok:continue
 z=pair['baseline']['z'];res={'sdss_id':sid}
 for k in ['baseline','dev']:
  r=pair[k];p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
  with fits.open(p,memmap=False) as h:
   co=h['COADD'].data;w=10**co['loglam']/(1+z);fl=np.asarray(co['flux'],float);iv=np.asarray(co['ivar'],float);good=(iv>0)&np.isfinite(fl)&(co['and_mask']==0)
   def window(a,b):return (w>=a)&(w<b)
   left=window(4750,4770)&good;right=window(5080,5100)&good;core=(window(4805,4840)|window(4885,4920))&good;allcore=(window(4805,4840)|window(4885,4920));cover=core.sum()/allcore.sum() if allcore.sum() else 0
   if cover<.8 or not left.any() or not right.any():res[k+'_proxy']='';continue
   xa=np.median(w[left]);xb=np.median(w[right]);fa=np.mean(fl[left]);fb=np.mean(fl[right]);continuum=fa+(fb-fa)*(w[core]-xa)/(xb-xa);area=float(np.sum((fl[core]-continuum)*np.gradient(w)[core]));res[k+'_proxy']=round(area/float(q[sid][k]['oiii_flux']),6)
 if res.get('baseline_proxy','')!='' and res.get('dev_proxy','')!='':res['difference']=round(res['dev_proxy']-res['baseline_proxy'],6)
 else:res['difference']=''
 out.append(res)
with open(sys.argv[4],'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['sdss_id','baseline_proxy','dev_proxy','difference']);w.writeheader();w.writerows(out)
print('diagnostic pairs',len(out),'complete',sum(x['difference']!='' for x in out))

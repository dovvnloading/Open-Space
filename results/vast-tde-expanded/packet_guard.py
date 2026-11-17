"""Metadata guard and synthetic-only conditional rise demonstration.
Never opens radio measurement parquet or real source flux.
"""
import csv,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def guard():
 p=json.loads((ROOT/'proposed-packet.json').read_text(encoding='utf-8'))
 if p['holdout_flux_access_authorized'] is not False:raise ValueError('holdout must remain sealed')
 for name,digest in p['inputs_sha256'].items():
  if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('input digest mismatch: '+name)
 rows=list(csv.DictReader((ROOT/'expanded-metadata.csv').open(newline='',encoding='utf-8')))
 if len(rows)!=134 or sum(r['status'] in ('candidate association; unverified','no catalog counterpart within 30 arcsec') for r in rows)!=103:raise ValueError('parent denominator mismatch')
 for sid,item in p['sources'].items():
  dev=item['development'];hold=item['holdout'];groups=dev+hold; flat=[x for group in groups for x in group]
  if len(dev)<2 or len(hold)!=2 or len(flat)!=len(set(flat)) or set(map(str,flat))!=set(item['sbids']):raise ValueError('invalid image split')
  if len(set(item['sbids'].values()))!=len(flat):raise ValueError('SBID overlap')
 controls=list(csv.DictReader((ROOT/'control-metadata.csv').open(newline='',encoding='utf-8')))
 for sid in p['sources']:
  r=[x for x in controls if x['target_source']==sid]
  if len(r)<20 or len(r)>200 or len({x['control_source'] for x in r})!=len(r):raise ValueError('invalid metadata controls')
 return p

def toy_statistics(early,hold,early_sigma,hold_sigma):
 e=np.asarray(early,float);h=np.asarray(hold,float);s=np.asarray(early_sigma,float);v=np.asarray(hold_sigma,float)
 if e.shape!=s.shape or len(e)<2 or h.shape!=(2,) or v.shape!=(2,) or np.any(s<=0) or np.any(v<=0) or not np.isfinite(np.r_[e,h,s,v]).all():raise ValueError('invalid input')
 w=1/s**2;mean=np.sum(w*e)/sum(w);return (h-mean)/np.sqrt(v**2+1/sum(w))

def synthetic():
 rng=np.random.default_rng(20261117);n=100000;early=rng.normal(size=(n,6));hold=rng.normal(size=(n,2));mean=early.mean(axis=1);z=(hold-mean[:,None])/np.sqrt(1+1/6)
 null=np.mean(np.all(z>=3,axis=1));injected=np.mean(np.all(z+4/np.sqrt(1+1/6)>=3,axis=1));print('synthetic null crossing',null,'injected crossing',injected,'n',n)
 assert null<0.002 and injected>0.4
 assert np.allclose(toy_statistics([1,1],[1,1],[1,1],[1,1]),0)

if __name__=='__main__':
 guard();synthetic();print('metadata guard and synthetic checks passed; no flux read')

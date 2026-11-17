"""Offline replay from post-unblinding row projections and frozen calculations."""
import csv,gzip,hashlib,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent;E=ROOT/'evidence';sys.path.insert(0,str(ROOT.parent/'vast-tde-expanded'))
from packet_guard import guard
from conditional_pipeline import phase_rows,develop,evaluate_holdout

def check(a,b,path='root'):
 if isinstance(a,dict) and isinstance(b,dict):
  if set(a)!=set(b):raise AssertionError(path+' keys')
  for k in a:check(a[k],b[k],path+'.'+str(k))
 elif isinstance(a,(list,tuple)) and isinstance(b,(list,tuple)):
  if len(a)!=len(b):raise AssertionError(path+' length')
  for i,(x,y) in enumerate(zip(a,b)):check(x,y,path+f'[{i}]')
 elif isinstance(a,bool) or isinstance(b,bool):
  if a!=b or type(a)!=type(b):raise AssertionError(path+' flag')
 elif isinstance(a,(int,float)) and isinstance(b,(int,float)):
  if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10):raise AssertionError(f'{path}: {a} != {b}')
 elif a!=b:raise AssertionError(f'{path}: {a!r} != {b!r}')

def rows(name,expected):
 p=E/name
 if hashlib.sha256(p.read_bytes()).hexdigest()!=expected['sha256']:raise ValueError('projection digest mismatch '+name)
 with gzip.open(p,'rt',encoding='utf-8') as f:r=[json.loads(x) for x in f]
 if len(r)!=expected['rows']:raise ValueError('projection row count mismatch')
 return r

def check_failed_preflight():
 raw=gzip.decompress((E/'failed-preflight-manifest.json.gz').read_bytes())
 if hashlib.sha256(raw).hexdigest()!='27f0c6f291bd4cdbd14461d72104fe03f81fafc7207fc620a8bfe6130bc0661a':raise ValueError('failed preflight digest mismatch')
 old=json.loads(raw)
 if old['git_head']!='4e32585258bd2419538e9698c8821f8910778443':raise ValueError('failed preflight head mismatch')
 blob=gzip.decompress((E/'failed-attributes-git-blob.txt.gz').read_bytes())
 reconstructed=gzip.decompress((E/'failed-attributes-crlf-reconstructed.txt.gz').read_bytes())
 if hashlib.sha256(blob).hexdigest()!='8c5d9c4b533c370055e534748b39c27cbe529de4769af068fafd695bdde3eede':raise ValueError('historical blob mismatch')
 if hashlib.sha256(reconstructed).hexdigest()!=old['worktree_sha256']['results/vast-tde-expanded/.gitattributes'] or reconstructed==blob or reconstructed.replace(b'\r\n',b'\n')!=blob:raise ValueError('historical checkout reconstruction mismatch')
 print('historical failed preflight and reconstructed EOL difference verified')

def main():
 check_failed_preflight()
 packet=guard();manifest=json.loads((E/'projection-manifest.json').read_text());observed=json.loads((E/'holdout-result.json').read_text());controls=list(csv.DictReader((ROOT.parent/'vast-tde-expanded/control-metadata.csv').open(newline='',encoding='utf-8')))
 for sid,cfg in packet['sources'].items():
  basefile=E/f'{sid}-development.json'
  if hashlib.sha256(basefile.read_bytes()).hexdigest()!=manifest['development_sha256'][sid]:raise ValueError('frozen baseline mismatch')
  expected=json.loads(basefile.read_text(encoding='utf-8'))['result'];pool={int(x['control_source']) for x in controls if x['target_source']==sid};ids={i for g in cfg['development'] for i in g}
  data=rows(f'{sid}-development.jsonl.gz',manifest['projections'][f'{sid}-development.jsonl.gz']);good,drops=phase_rows(data,ids,pool|{int(sid)});result=develop(int(sid),packet,good,pool);result['drops']=drops+result['drops'];check(expected,json.loads(json.dumps(result)),sid+'.development')
  ids={i for g in cfg['holdout'] for i in g};pool={int(k) for k in expected['control_baselines']};data=rows(f'{sid}-holdout.jsonl.gz',manifest['projections'][f'{sid}-holdout.jsonl.gz']);good,drops=phase_rows(data,ids,pool|{int(sid)});result=evaluate_holdout(expected,good);result['drops']=drops+result['drops'];actual=observed['sources'][sid];check({k:v for k,v in actual.items() if k!='status'},json.loads(json.dumps(result)),sid+'.holdout');print(sid,'development and holdout reproduced',result['result'])
 print('Offline replay passed; no network or native Parquet required.')
if __name__=='__main__':main()

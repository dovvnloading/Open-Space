"""Predeclared pseudo-test on last two DEVELOPMENT sessions only."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'vast-tde-expanded'))
from packet_guard import guard
from read_development import read_development
from conditional_pipeline import develop,evaluate_holdout

PSEUDO_SHA='aef29ad1b93ac5f960e61e46fb4ba67528c48d49e4bacfc48582641a72867cc9'

def run(pseudo_path,parquet,images,out):
 raw=Path(pseudo_path).read_bytes()
 if hashlib.sha256(raw).hexdigest()!=PSEUDO_SHA:raise ValueError('pseudo target list mismatch')
 config=json.loads(raw);packet=guard()
 if config['actual_holdout_authorized'] is not False:raise ValueError('actual holdout must be sealed')
 if Path(out).exists():raise FileExistsError('output must be new')
 controls=list(csv.DictReader((ROOT.parent/'vast-tde-expanded/control-metadata.csv').open()))
 results=[]
 for sid in (14043503,13894502,14481388):
  rows,drops=read_development(sid,parquet,images)
  pool={int(r['control_source']) for r in controls if int(r['target_source'])==sid}
  dev=packet['sources'][str(sid)]['development'];train=dev[:-2];test=dev[-2:]
  trainids={i for g in train for i in g};testids={i for g in test for i in g}
  if testids & {i for g in packet['sources'][str(sid)]['holdout'] for i in g}:raise ValueError('actual holdout overlap')
  for pseudo in config['pseudotargets'][str(sid)]:
   record={'parent_source':sid,'pseudo_source':pseudo,'status':'inconclusive','input_quality_drops':len(drops)}
   try:
    candidate_controls=pool-{pseudo};allowed=candidate_controls|{pseudo}
    training=[r for r in rows if r['image_id'] in trainids and r['source'] in allowed]
    held=[r for r in rows if r['image_id'] in testids and r['source'] in allowed]
    fake={'sources':{str(pseudo):{'development':train,'holdout':test}}}
    result=develop(pseudo,fake,training,candidate_controls)
    values={}
    for injection in (0.,.5,1.,2.):
     amended=[dict(r,flux_peak=float(r['flux_peak'])+injection) if r['source']==pseudo else r for r in held]
     values[str(injection)]=evaluate_holdout(result,amended)['result']
    record.update(status='measured',controls=len(result['control_baselines']),tests=values)
   except Exception as exc:record.update(error_type=type(exc).__name__,error=str(exc))
   results.append(record)
 with Path(out).open('x',encoding='utf-8') as f:json.dump(results,f,indent=2,allow_nan=False);f.write('\n')
 return results

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--pseudo',required=True);ap.add_argument('--parquet',required=True);ap.add_argument('--images',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();r=run(a.pseudo,a.parquet,a.images,a.output);print('completed',len(r),'measured',sum(x['status']=='measured' for x in r))

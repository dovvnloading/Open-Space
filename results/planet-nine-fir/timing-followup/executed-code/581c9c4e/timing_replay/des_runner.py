"""DES retrospective runner. Verify frozen inputs before opening image FITS."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
from preflight import verify_mapping
from scenario_states import states_for_rows, ROOT, p
from controls import trial_seed
from des_packet import verify as verify_packet

sys.path.insert(0,str(ROOT/'archives'))
import planet_injections as pi
import decam_difference as dd

LABELS=('DECam_605285','DECam_605571','mask_605285','mask_605571','template_605571')

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def require_d(path):
 if not path.is_absolute():raise ValueError('absolute paths required')
 path=path.resolve()
 if path.drive.upper()!='D:':raise ValueError('paths must be absolute on D:')
 return path

def configure(scenario,mapping_path,mapping_sha,cache,out):
 if scenario not in ('historical','dec1','dec2'):raise ValueError('unknown scenario')
 mapping_path,cache,out=map(require_d,(mapping_path,cache,out))
 if out.exists():raise FileExistsError(out)
 if out==cache or out in cache.parents or cache in out.parents:raise ValueError('cache and output must not overlap')
 if not cache.is_dir():raise ValueError('missing input cache')
 mapping=verify_mapping(mapping_path,mapping_sha)
 for path,expected in [('orbits/orbit_input.json',mapping['input_sha256']['orbits\\orbit_input.json']),('orbits/propagate.py',mapping['input_sha256']['orbits\\propagate.py'])]:
  if digest(ROOT/path)!=expected:raise ValueError('mapped source changed: '+path)
 receipts=json.loads((ROOT/'archives'/'decam_inputs.json').read_bytes())
 for label in LABELS:
  found=[r for r in receipts if r['label']==label]
  if len(found)!=1:raise ValueError('invalid source manifest')
  src=cache/(label+'.response');r=found[0]
  if not src.is_file() or src.stat().st_size!=r['size'] or digest(src)!=r['sha256']:raise ValueError('input digest mismatch: '+label)
 return mapping

def execute(scenario,mapping_path,mapping_sha,cache,out,manifest,manifest_sha,auth,auth_sha):
 verify_packet(require_d(manifest),manifest_sha,require_d(auth),auth_sha,scenario,mapping_path,mapping_sha,cache,out)
 mapping=configure(scenario,mapping_path,mapping_sha,cache,out)
 inputs=json.loads((ROOT/'orbits'/'orbit_input.json').read_bytes())
 rows=states_for_rows(mapping,inputs,scenario)
 outcomes={(r['experiment'],r['draw']):{'experiment':r['experiment'],'draw':r['draw'],'state_status':r['status'],'image_status':'not_tested','preplant_eligible':False,'models':[{'postplant_eligible':False,'recovered':None,'reason':'not_tested'} for _ in range(21)]} for r in rows}
 if len(outcomes)!=800:raise ValueError('duplicate identity')
 original_family=p.family; original_measure=pi.measure
 def scenario_family(inp,*args,**kwargs):
  copy=dict(inp,akari_midpoint_utc=__import__('timing_scenarios').SCENARIOS[scenario][0])
  return original_family(copy,*args,**kwargs)
 def population(stress,wcs,xs,ys,times):
  exp='centroid_stress' if stress else 'published_centroids'; selected=[r for r in rows if r['experiment']==exp]
  if [r['draw'] for r in selected]!=list(range(400)):raise ValueError('draw order invalid')
  table=np.zeros((400,13),dtype=float);table[:,4:6]=1.0
  bound=[i for i,r in enumerate(selected) if r['status']=='bound']
  if bound:
   state=np.array([selected[i]['state'] for i in bound]);st=(state[:,:3],state[:,3:6],state[:,6],{})
   dist=pi.distance_at(p,st,times[0]);positions=[]
   for t in times:
    ra,de=p.sky(st,t);x,y=wcs.world_to_pixel_values(ra,de);positions.extend([x-xs[0],y-ys[0]])
   table[bound]=np.column_stack((*positions,*dist,*state.T))
  return table
 def measure(*args,**kwargs):
  kwargs.update(population_override=population,draw_outcomes=outcomes,noise_seed_for=trial_seed)
  return original_measure(*args,**kwargs)
 # Exclusive creation occurs only after every digest guard; never reuses an old run.
 out.mkdir(parents=True,exist_ok=False)
 (out/'mapped_states.json').write_text(json.dumps(rows,indent=2)+'\n')
 (out/'run_manifest.json').write_text(json.dumps({'scenario':scenario,'mapping_sha256':mapping_sha,'input_cache':str(cache),'input_sha256':{label:digest(cache/(label+'.response')) for label in LABELS},'code_sha256':{str(path.relative_to(ROOT)):digest(path) for path in (Path(__file__),ROOT/'archives'/'planet_injections.py',ROOT/'archives'/'decam_difference.py',ROOT/'orbits'/'propagate.py')},'python':sys.version,'noise':'keyed by experiment/draw/model/exposure; shared keys do not imply identical Poisson realization','retrospective':True},indent=2)+'\n')
 status={'status':'running'}
 try:
  p.family=scenario_family;pi.measure=measure;dd.RNG=np.random.default_rng(192721)
  dd.main(cache,out)
  status={'status':'complete'}
 except BaseException as exc:
  status={'status':'failed','error_type':type(exc).__name__,'error':str(exc)}
  raise
 finally:
  p.family=original_family;pi.measure=original_measure
  (out/'draw_outcomes.json').write_text(json.dumps(list(outcomes.values()),indent=2)+'\n')
  (out/'run_status.json').write_text(json.dumps(status,indent=2)+'\n')

def main():
 a=argparse.ArgumentParser(description=__doc__)
 for name in ('scenario','mapping','mapping-sha256','cache','out','manifest','manifest-sha256','auth','auth-sha256'):a.add_argument('--'+name,required=True)
 q=a.parse_args();execute(q.scenario,Path(q.mapping),q.mapping_sha256,Path(q.cache),Path(q.out),Path(q.manifest),q.manifest_sha256,Path(q.auth),q.auth_sha256)
if __name__=='__main__':main()

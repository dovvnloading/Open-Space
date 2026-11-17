"""Compute scenario states by frozen pre-filter grid identity without pixel access."""
import json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'orbits'))
sys.path.insert(0,str(ROOT/'archives'))
import propagate as p
from planet_injections import offset
from timing_scenarios import SCENARIOS

# Use identical historical integrator while retaining unbound solutions.
source=(ROOT/'orbits'/'propagate.py').read_text(encoding='utf-8')
needle='    keep = energy < 0\n'
if source.count(needle)!=1:raise RuntimeError('unrecognized integrator source')
ns={'__name__':'all_states','__file__':str(ROOT/'orbits'/'propagate.py')}
exec(compile(source.replace(needle,'    bound_mask = energy < 0\n    keep = np.ones(len(energy), dtype=bool)\n').replace('    info = dict(\n','    info = dict(bound_mask=bound_mask.tolist(),\n'),ns['__file__'],'exec'),ns)

def states_for_rows(mapping,inputs,scenario):
 if scenario not in ('historical','dec1','dec2'):raise ValueError('unknown scenario')
 rows=mapping['rows'];out=[None]*len(rows)
 pa=np.deg2rad(79);iras=[(0,0),(18*np.sin(pa),18*np.cos(pa)),(-18*np.sin(pa),-18*np.cos(pa)),(7*np.cos(pa),-7*np.sin(pa)),(-7*np.cos(pa),7*np.sin(pa))]
 akari=[(0,0),(40,0),(-40,0),(0,40),(0,-40)]
 groups={}
 for n,row in enumerate(rows):
  if row['status']!='unique' or len(row['matches'])!=1:
   out[n]={'experiment':row['experiment'],'draw':row['draw'],'status':row['status']};continue
  match=row['matches'][0];i,j=match['iras_offset'],match['akari_offset'];ix=match['grid_index']
  if not 0<=i<5 or not 0<=j<5 or len(ix)!=3 or any(not 0<=x<n for x,n in zip(ix,(13,32,56))):raise ValueError('invalid draw identity')
  groups.setdefault((i,j),[]).append((n,ix))
 for (i,j),members in groups.items():
  inp=dict(inputs,akari_midpoint_utc=SCENARIOS[scenario][0],iras_radec_deg=offset(inputs['iras_radec_deg'],*iras[i],p),akari_radec_deg=offset(inputs['akari_radec_deg'],*akari[j],p))
  parameters=np.array([[ix[0]/12,ix[1]/31,ix[2]/55] for _,ix in members])
  r,v,t,info=ns['family'](inp,parameters=parameters)
  for k,(n,ix) in enumerate(members):
   out[n]={'experiment':rows[n]['experiment'],'draw':rows[n]['draw'],'status':'bound' if info['bound_mask'][k] else 'unbound','state':[*r[k].tolist(),*v[k].tolist(),float(t[k])]}
 return out

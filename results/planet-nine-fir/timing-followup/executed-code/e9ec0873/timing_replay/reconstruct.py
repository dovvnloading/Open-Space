"""Reconstruct historical draw identities, without reading image products."""
import argparse, csv, gzip, hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'orbits'))
sys.path.insert(0, str(ROOT/'archives'))
import propagate as original
from planet_injections import offset

# Instrument only the bound-mask return, preserving all historical calculations.
source = (ROOT/'orbits'/'propagate.py').read_text(encoding='utf-8')
needle = '    info = dict(\n'
assert source.count(needle) == 1
namespace = {'__name__': 'instrumented_propagate', '__file__': str(ROOT/'orbits'/'propagate.py')}
exec(compile(source.replace(needle, '    info = dict(grid_indices=np.flatnonzero(keep).tolist(),\n'), namespace['__file__'], 'exec'), namespace)
family = namespace['family']
TOL = np.array([1e-8]*3 + [1e-12]*3 + [1e-8])
KEYS = ['x_au','y_au','z_au','vx_au_d','vy_au_d','vz_au_d','reference_emission_tt_days_since_j2000']

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--out',required=True,type=Path);args=a.parse_args()
 args.out.mkdir(parents=True,exist_ok=False)
 paths=[ROOT/'archives'/'planet_draws.csv.gz',ROOT/'orbits'/'orbit_input.json',ROOT/'orbits'/'propagate.py',ROOT/'archives'/'planet_injections.py',Path(__file__)]
 rows=list(csv.DictReader(gzip.open(paths[0],'rt')))
 if len(rows)!=800 or len({(r['experiment'],r['draw']) for r in rows})!=800: raise ValueError('invalid frozen row membership')
 inputs=json.loads(paths[1].read_bytes());pa=np.deg2rad(79)
 iras=[(0,0),(18*np.sin(pa),18*np.cos(pa)),(-18*np.sin(pa),-18*np.cos(pa)),(7*np.cos(pa),-7*np.sin(pa)),(-7*np.cos(pa),7*np.sin(pa))]
 akari=[(0,0),(40,0),(-40,0),(0,40),(0,-40)]
 values=np.array([[float(r[k]) for k in KEYS] for r in rows]);matches=[[] for r in rows]
 for i,ii in enumerate(iras):
  for j,aa in enumerate(akari):
   inp=dict(inputs,iras_radec_deg=offset(inputs['iras_radec_deg'],*ii,original),akari_radec_deg=offset(inputs['akari_radec_deg'],*aa,original))
   states=family(inp,13,32,56);table=np.column_stack((*states[0].T,*states[1].T,states[2]))
   tree=cKDTree(table/TOL);nearest=tree.query_ball_point(values/TOL,r=np.sqrt(7),p=2)
   for rownum,found in enumerate(nearest):
    if rows[rownum]['experiment']=='published_centroids' and (i,j)!=(0,0):continue
    for x in found:
     delta=np.abs(values[rownum]-table[x]);
     if np.all(delta<=TOL):
      flat=states[3]['grid_indices'][x];t,r,q=np.unravel_index(flat,(13,32,56))
      matches[rownum].append({'iras_offset':i,'akari_offset':j,'grid_index':[int(t),int(r),int(q)],'max_normalized_residual':float(np.max(delta/TOL))})
 output=[]
 for row,hit in zip(rows,matches):
  output.append({'experiment':row['experiment'],'draw':int(row['draw']),'status':'unique' if len(hit)==1 else ('unmapped' if not hit else 'ambiguous'),'matches':hit})
 manifest={'code_revision':__import__('subprocess').check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'tolerance_position_au':1e-8,'tolerance_velocity_au_per_day':1e-12,'tolerance_time_day':1e-8,'status_counts':{s:sum(r['status']==s for r in output) for s in ['unique','ambiguous','unmapped']},'rows':output}
 (args.out/'mapping.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(manifest['status_counts']))
if __name__=='__main__':main()

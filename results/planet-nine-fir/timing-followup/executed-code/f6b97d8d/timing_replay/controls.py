"""Pure helpers for retrospective injection accounting and deterministic noise."""
import hashlib
import numpy as np

def trial_seed(experiment,draw,model,exposure,realization):
 # Deliberately independent of scenario for common random numbers; a changed
 # Poisson mean/shape still does not imply identical photon realizations.
 key='\0'.join(map(str,('pr143-retrospective-v1',experiment,draw,model,exposure,realization)))
 return int.from_bytes(hashlib.sha256(key.encode()).digest()[:16],'big')

def validate_coverage(x,y,width,height,margin,good_mask):
 if not all(isinstance(v,(int,np.integer)) for v in (x,y,width,height,margin)):raise ValueError('integer geometry required')
 if width<=0 or height<=0 or margin<0 or good_mask.shape!=(height,width):raise ValueError('invalid coverage geometry')
 if x-margin<0 or y-margin<0 or x+margin>=width or y+margin>=height:return 'outside_range'
 if not np.all(good_mask[y-margin:y+margin+1,x-margin:x+margin+1]):return 'masked'
 return 'covered'

def partition(rows):
 allowed={'unmapped','ambiguous','unbound','outside_range','masked','eligible','recovered','not_recovered'}
 out={k:0 for k in allowed}
 seen=set()
 for row in rows:
  key=(row['experiment'],row['draw']);s=row['status']
  if key in seen or s not in allowed:raise ValueError('duplicate draw or invalid status')
  seen.add(key);out[s]+=1
 return out

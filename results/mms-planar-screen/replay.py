"""Offline scalar gate recount on deposited authorized MMS1 projection."""
import json,hashlib,math
from pathlib import Path
import numpy as np
D=Path(__file__).resolve().parent;raw=(D/'evidence/mms1-projection.npz').read_bytes();receipt=json.loads((D/'evidence/projection-receipt.json').read_text());assert hashlib.sha256(raw).hexdigest()==receipt['sha256']
with np.load(D/'evidence/mms1-projection.npz',allow_pickle=False) as p:t=p['t'];b=p['field'];good=p['good']
counts={'tested_1s_windows':0,'quality_cadence_windows':0,'field_endpoint_windows':0,'rotation_windows':0};maximum=0.
for i in range(len(t)):
 j=int(np.searchsorted(t,t[i]+1.))
 if j>=len(t) or t[j]-t[i]>1.02:continue
 counts['tested_1s_windows']+=1
 lo=int(np.searchsorted(t,t[i]-5));hi=int(np.searchsorted(t,t[i]+5,side='right'))
 if t[i]-t[0]<5 or t[-1]-t[i]<5 or hi-lo<1000 or not good[lo:hi].all() or np.max(np.diff(t[lo:hi]))>.02 or not np.isfinite(b[lo:hi]).all():continue
 counts['quality_cadence_windows']+=1
 x=b[i];y=b[j];nx=float(np.linalg.norm(x));ny=float(np.linalg.norm(y))
 if min(nx,ny)<3:continue
 counts['field_endpoint_windows']+=1
 angle=math.degrees(math.acos(float(np.clip(x@y/(nx*ny),-1,1))));maximum=max(maximum,angle)
 if angle>=45:counts['rotation_windows']+=1
expected=json.loads((D/'evidence/development.json').read_text())['gates']
for k,v in counts.items():assert v==expected[k],(k,v,expected[k])
assert counts['rotation_windows']==0
print(json.dumps({'gate_counts':counts,'max_rotation_deg_descriptive':maximum,'status':'reproduced_from_projection'},indent=2))

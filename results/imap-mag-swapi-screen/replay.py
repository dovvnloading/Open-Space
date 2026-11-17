"""Independent offline recount on the first-12-hour MAG projection."""
from pathlib import Path
import json,hashlib,numpy as np
D=Path(__file__).resolve().parent;path=D/'evidence/mag-projection.npz';receipt=json.loads((D/'evidence/projection-receipt.json').read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==receipt['sha256']
with np.load(path,allow_pickle=False) as f:t=f['t'];b=f['b'];good=f['good']
s=[];audit=[]
for k in range(144):
 take=(t>=k*300)&(t<(k+1)*300);tt=t[take];ok=good[take]&np.isfinite(b[take]).all(axis=1);v=tt[ok];count=len(v)>=540;edge=bool(len(v) and v[0]<=k*300+2 and v[-1]>=(k+1)*300-2);gap=bool(len(v)>1 and np.max(np.diff(v))<=2);usable=count and edge and gap
 audit.append({'bin':k,'attempted':int(take.sum()),'valid':len(v),'count_ok':count,'edge_ok':edge,'gap_ok':gap,'usable':usable});s.append(float(np.median(np.linalg.norm(b[take][ok],axis=1))) if usable else None)
gates={'triplets_tested':0,'quality_triplets':0,'both_controls_ge_3nT':0,'absolute_drop_ge_1nT':0,'both_controls_fraction_gate':0,'stopped_at_first_selection':False};selection=None
for k in range(1,143):
 gates['triplets_tested']+=1;l,m,r=s[k-1:k+2]
 if None in (l,m,r):continue
 gates['quality_triplets']+=1
 if min(l,r)<3:continue
 gates['both_controls_ge_3nT']+=1
 if (l+r)/2-m<1:continue
 gates['absolute_drop_ge_1nT']+=1
 if m>.8*min(l,r):continue
 gates['both_controls_fraction_gate']+=1;gates['stopped_at_first_selection']=True;selection=k;break
saved=json.loads((D/'evidence/development.json').read_text());assert gates==saved['gates'] and audit==saved['bin_audit'] and selection is None and saved['status']=='no_candidate';print(json.dumps({'status':'reproduced_from_projection','valid_bins':sum(v is not None for v in s),'gates':gates},indent=2))

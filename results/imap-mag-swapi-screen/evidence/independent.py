"""Separate native-CDF recount of the authorized IMAP MAG-only screen."""
import hashlib
import json
from pathlib import Path
import cdflib
import numpy as np

BASE = Path('D:/Open-Space/campaigns')
SOURCE = BASE/'heliophysics-schema-20261117/raw/imap_mag_l2_norm-gse_20260731_v002.0003.cdf'
PRIMARY = BASE/'imap-mag-dev-0695ac72-96cc21b-20261117/development.json'
OUTPUT = BASE/'imap-independent-dev-0695ac72-20261117.json'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
assert sha(SOURCE) == 'f587e5ea2969418a54738580f8289ba2ef32a1d755b6af377c6ad3a8dc774727'
c = cdflib.CDF(str(SOURCE))
epoch = c.varget('epoch')
assert c.varinq('epoch').Data_Type_Description == 'CDF_TIME_TT2000'
assert np.all(np.diff(epoch)>0)
start = int(cdflib.cdfepoch.compute_tt2000([2026,7,31,0,0,0,0,0,0]))
indices = np.flatnonzero((epoch >= start) & (epoch < start+43_200_000_000_000))
i,j = int(indices[0]),int(indices[-1])
assert np.array_equal(indices,np.arange(i,j+1))
b = c.varget('b_gse',startrec=i,endrec=j)
attr = c.varattsget('b_gse')
assert attr['UNITS']=='nT' and list(c.varinq('b_gse').Dim_Sizes)==[3]
flags = c.varget('quality_flags',startrec=i,endrec=j)
bits = c.varget('quality_bitmask',startrec=i,endrec=j)
valid = ((flags==0)&(bits==0)&np.isfinite(b).all(axis=1)&
         (b!=attr['FILLVAL']).all(axis=1)&(b>=attr['VALIDMIN']).all(axis=1)&
         (b<=attr['VALIDMAX']).all(axis=1))
elapsed = epoch[indices]-start
strength = np.sqrt(np.sum(b*b,axis=1))
audit,medians = [],[]
for k in range(144):
    lo,hi = k*300_000_000_000,(k+1)*300_000_000_000
    take = (elapsed>=lo)&(elapsed<hi)
    keep = take&valid
    stamps = elapsed[keep]
    count_ok = int(keep.sum())>=540
    edge_ok = bool(len(stamps) and stamps[0]<=lo+2_000_000_000 and stamps[-1]>=hi-2_000_000_000)
    gap_ok = bool(len(stamps)>1 and np.max(np.diff(stamps))<=2_000_000_000)
    usable = bool(count_ok and edge_ok and gap_ok)
    audit.append(dict(bin=k,attempted=int(take.sum()),valid=int(keep.sum()),
                      count_ok=count_ok,edge_ok=edge_ok,gap_ok=gap_ok,usable=usable))
    medians.append(float(np.median(strength[keep])) if usable else None)
gates = dict(triplets_tested=0,quality_triplets=0,both_controls_ge_3nT=0,
             absolute_drop_ge_1nT=0,both_controls_fraction_gate=0,stopped_at_first_selection=False)
selected = None
for k in range(1,143):
    gates['triplets_tested']+=1
    l,m,r=medians[k-1:k+2]
    if any(x is None for x in (l,m,r)): continue
    gates['quality_triplets']+=1
    if min(l,r)<3: continue
    gates['both_controls_ge_3nT']+=1
    if (l+r)/2-m<1: continue
    gates['absolute_drop_ge_1nT']+=1
    if m>.8*min(l,r): continue
    gates['both_controls_fraction_gate']+=1
    gates['stopped_at_first_selection']=True
    selected=k
    break
coverage=dict(first_record=i,last_record=j,records=len(indices),good=int(valid.sum()))
primary=json.loads(PRIMARY.read_text())
assert coverage==primary['coverage']
assert audit==primary['bin_audit']
assert gates==primary['gates']
assert selected is None and primary['status']=='no_candidate'
result=dict(status='independent_native_null_reproduced',source_sha256=sha(SOURCE),
            primary_sha256=sha(PRIMARY),code_sha256=sha(Path(__file__)),
            coverage=coverage,valid_bins=sum(x['usable'] for x in audit),gates=gates,
            bin_audit=audit,scope='Same first12h MAG observations, separate calculation. No SWAPI or later MAG physical values read.')
with OUTPUT.open('x',encoding='utf-8') as out:
    json.dump(result,out,indent=2,allow_nan=False)
    out.write('\n')
print(json.dumps({k:v for k,v in result.items() if k!='bin_audit'},indent=2))

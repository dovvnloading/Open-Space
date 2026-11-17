"""Apply the frozen same-band observer-frame thresholds to ALeRCE positive difference PSF flux."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

QUALITY = ["isdipole","glint_trail","psfflux_flag","pixelflags_edge","pixelflags_crcenter",
           "pixelflags_saturatedcenter","pixelflags_streakcenter","pixelflags_injectedcenter",
           "pixelflags_injected_templatecenter"]
def screen(directory,out):
    # Validate all frozen raw replies before parsing or creating outputs.
    for key in ['parent'] + [f'day_{d}' for d in range(61345, 61361)]:
        receipt = json.loads((directory / (key + '.json')).read_text())
        raw = gzip.decompress((directory / (key + '.csv.gz')).read_bytes())
        if receipt.get('key') != key or receipt.get('status') != 200 or len(raw) != receipt.get('bytes') or hashlib.sha256(raw).hexdigest() != receipt.get('sha256'):
            raise ValueError('Frozen input mismatch: ' + key)
    if out.resolve().is_relative_to(Path(__file__).resolve().parents[1]):
        raise ValueError('Output must be outside the frozen result')
    if out.exists():
        raise FileExistsError(out)
    parent={}
    with gzip.open(directory/'parent.csv.gz','rt') as f:
        parent={r['oid']:r for r in csv.DictReader(f)}
    groups=defaultdict(dict)
    days={}
    for day in range(61345,61361):
        key=f'day_{day}'
        with gzip.open(directory/(key+'.csv.gz'),'rt') as f: rows=list(csv.DictReader(f))
        receipt=json.loads((directory/(key+'.json')).read_text())
        if receipt.get('rows')!=len(rows): raise RuntimeError(f'Missing rows {key}')
        days[day]=len(rows)
        for r in rows:
            prior=groups[r['oid']].setdefault(r['measurement_id'],r)
            if r!=prior: raise RuntimeError('Conflicting measurements')
    if set(groups)-set(parent): raise RuntimeError('Joined object outside parent')
    pairs=[]; counts=Counter(); positive=set(); sampled=set(); rapid=set(); clean=set()
    for oid,table in groups.items():
        bands=defaultdict(list)
        for r in table.values():
            try: flux,err=float(r['psfflux']),float(r['psffluxerr'])
            except ValueError: continue
            if not (math.isfinite(flux) and math.isfinite(err) and err>0 and flux>=5*err): continue
            if r['isnegative']=='True': continue
            positive.add(oid)
            bands[r['band']].append((float(r['mjd']),flux,err,r))
        for band,rows in bands.items():
            rows.sort(key=lambda item:item[0])
            for i,(t1,f1,e1,a) in enumerate(rows):
                for t2,f2,e2,b in rows[i+1:]:
                    dt=t2-t1
                    if dt>1: break
                    if dt<0.02 or a['visit']==b['visit']:continue
                    sampled.add(oid); counts['eligible_pairs']+=1
                    dm=2.5*math.log10(f1/f2)
                    sd=2.5/math.log(10)*math.hypot(e1/f1,e2/f2)
                    if abs(dm)<0.30 or abs(dm)/dt<=1 or abs(dm)<3*sd: continue
                    rapid.add(oid)
                    flags=[q for q in QUALITY if a[q]=='True' or b[q]=='True']
                    if any(r['ssobjectid'] not in ['', '0'] or r['diaobjectid'] in ['','0'] for r in (a,b)):flags.append('ss_association')
                    if not flags:clean.add(oid)
                    pairs.append(dict(oid=oid,ra=a['meanra'],dec=a['meandec'],band=band,
                        visit1=a['visit'],visit2=b['visit'],measurement1=a['measurement_id'],measurement2=b['measurement_id'],
                        t1=t1,t2=t2,dt_d=dt,dm_mag=dm,rate_mag_d=dm/dt,sigma_mag=sd,
                        snr1=f1/e1,snr2=f2/e2,rel1=a['reliability'],rel2=b['reliability'],flags='|'.join(flags)))
    out.mkdir(parents=True,exist_ok=False)
    summary=dict(parent_objects=len(parent),joined_objects=len(groups),joined_unique_measurements=sum(map(len,groups.values())),
            absent_parent_objects=len(set(parent)-set(groups)),parents_with_one_positive=len(positive),
            eligible_pair_objects=len(sampled),eligible_pairs=counts['eligible_pairs'],
            rapid_pair_objects=len(rapid),rapid_pairs=len(pairs),objects_with_unflagged_pair=len(clean),
            unflagged_pairs=sum(not r['flags'] for r in pairs),joined_rows_by_first_mjd_partition=days)
    print(json.dumps(summary,indent=2)); (out/'funnel.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (out/'pairs.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]) if pairs else ['oid']);w.writeheader();w.writerows(sorted(pairs,key=lambda r:(r['oid'],r['t1'],r['t2'])))
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('data',type=Path);ap.add_argument('out',type=Path);ap.add_argument('--flux-unit', choices=['nJy'], default='nJy');a=ap.parse_args()
    screen(a.data,a.out)

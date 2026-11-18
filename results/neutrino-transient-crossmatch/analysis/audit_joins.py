import csv,gzip,json,sys
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'analysis'))
import run
for archival in (False,True):
 reg=run.rows('ice/gcn_alerts.csv'); active={r['run_event'] for r in reg if r['status']=='active'}
 maps=[json.loads(x) for x in (root/'ice/gcn_hpd90.jsonl').read_text().splitlines() if json.loads(x)['run_event'] in active]
 with gzip.open(root/'ice/archival_maps.jsonl.gz','rt') as f:
  maps += [m for m in map(json.loads,f) if archival or run.days(m['utc'])>=run.days('2023-07-01T00:00:00Z')]
 sources,_=run.sources(archival); dirname='archival' if archival else 'modern'
 expected=list(csv.DictReader((root/'analysis'/dirname/'native_pairs.csv').open(newline='')))
 found=[]; support={}
 for fam, cat in sources.items():
  lo,hi=run.WINDOWS.get(fam,run.WINDOWS['mir']); eligible_pairs=0; ts_alerts=set(); ts_sources=set()
  for m in maps:
   mask=run.prepared_region(m)
   for j,s in enumerate(cat):
    lag=run.days(m['utc'])-s['day']
    if lo <= lag <= hi:
     eligible_pairs+=1; ts_alerts.add(m['run_event']); ts_sources.add(j)
     if bool(mask(np.array([s['ra']]),np.array([s['dec']]))[0]): found.append((fam,m['run_event'],s['source'],s['epoch']))
  support[fam]=(len(cat),len(ts_alerts),len(ts_sources),eligible_pairs)
 exp=[(r['family'],r['run_event'],r['source'],r['epoch']) for r in expected]
 assert sorted(found)==sorted(exp),(dirname,found,exp)
 coverage=list(csv.DictReader((root/'analysis'/dirname/'coverage.csv').open(newline='')))
 for row in coverage:
  actual=support.get(row['family'],(0,0,0,0))
  want=tuple(int(row[k]) for k in ('epochs','alerts_with_temporal_support','epochs_with_temporal_support','time_eligible_pairs'))
  assert actual==want,(dirname,row['family'],actual,want)
 print(dirname,'independent nested alert/source joins match',len(found),'pairs; support',support,flush=True)

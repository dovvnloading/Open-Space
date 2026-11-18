import csv,gzip,json,sys
from collections import defaultdict,Counter
from datetime import datetime
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(root/'analysis')); import run
for archival in (False,True):
 registry={r['run_event']:r for r in run.rows('ice/gcn_alerts.csv') if r['status']=='active'}
 maps=[]
 for line in (root/'ice/gcn_hpd90.jsonl').read_text().splitlines():
  m=json.loads(line)
  if m['run_event'] in registry: m['stream']=registry[m['run_event']]['stream']; maps.append(m)
 with gzip.open(root/'ice/archival_maps.jsonl.gz','rt') as f:
  maps.extend(m for m in map(json.loads,f) if archival or run.days(m['utc'])>=run.days('2023-07-01T00:00:00Z'))
 maps.sort(key=lambda m:m['utc']); n=len(maps); groups=defaultdict(list)
 for i,m in enumerate(maps):
  d=datetime.fromisoformat(m['utc'].replace('Z','+00:00')); groups[(d.year,(d.month-1)//3+1,m['stream'].upper().split('-')[-1])].append(i)
 rng=np.random.default_rng(20261120 if archival else 20261118); perms=np.tile(np.arange(n),(9999,1))
 for ids in groups.values(): perms[:,ids]=rng.permuted(perms[:,ids],axis=1)
 areas=np.array([m.get('published_area_deg2') or m['area_deg2'] for m in maps]); sources,_=run.sources(archival)
 saved=list(csv.DictReader((root/'analysis'/('archival' if archival else 'modern')/'null_histograms.csv').open(newline='')))
 for family in (('gamma','mir') if archival else ('optical','xray','gamma')):
  cat=sources[family]; lo,hi=run.WINDOWS[family]; times=np.array([s['day'] for s in cat]); ras=np.array([s['ra'] for s in cat]); decs=np.array([s['dec'] for s in cat]); stamps=np.array([run.days(m['utc']) for m in maps]); eligible=(stamps[:,None]-times>=lo)&(stamps[:,None]-times<=hi)
  # Explicit alert-geometry x source-time matrix: each cell is a physical contour/date cross-join.
  cross=np.zeros((n,n),bool)
  for i,m in enumerate(maps):
   sky=run.prepared_region(m)(ras,decs)
   for j in np.flatnonzero(sky): cross[i] |= eligible[:,j]
  hits=cross[np.arange(n)[None,:],perms]
  for label,values in [('distinct_alerts',hits.sum(axis=1)),('maximum_inverse_area_deg-2',np.max(np.where(hits,1/areas,0),axis=1))]:
   actual=Counter(map(float,values)); expected={float(r['value']):int(r['frequency']) for r in saved if r['family']==family and r['null']=='time' and r['statistic']==label}
   assert dict(actual)==expected,(archival,family,label)
  print(('archival' if archival else 'modern'),family,'9999 time trial histograms match',flush=True)

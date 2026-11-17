import csv,json,math,statistics,hashlib,platform
from pathlib import Path
import pyarrow.parquet as pq
P=Path(r'D:/Open-Space/worktrees/vast-holdout-review-20261117/results/vast-tde-expanded')
D=Path(r'D:/Open-Space/campaigns/vast-tde-dev-f323ab5-21ccdd9-20261117/run')
H=Path(r'D:/Open-Space/campaigns/vast-tde-holdout-2ca26d3-28ce15c-20261117')
O=Path(r'D:/Open-Space/campaigns/vast-tde-independent-holdout-2ca26d3')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(H/'holdout-authorization.json')=='c1d8be3f585f6666d0e04baff0d8ebac309d9671545de5b10ef6186821b51b0a'
assert sha(H/'run/holdout-result.json')=='a66cf11f370c5d6258e21cf2d85ea36877bb1f05dbaa46caaecd5642e8b08fe2'
assert sha(D/'summary.json')=='a851263393a0f59609fa14da6443e9c9bb6b751b94c4821b96da9af4f6e4bb60'
auth=json.loads((H/'holdout-authorization.json').read_text());summary=json.loads((D/'summary.json').read_text());packet=json.loads((P/'proposed-packet.json').read_text())
images=Path(r'D:/Open-Space/campaigns/vast-tde-expanded-20261117/image_df.csv');parquet=Path(r'D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_measurement_df.parquet')
assert sha(images)=='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195';assert sha(parquet)==packet['measurement_parquet_sha256']
meta={int(r['image_id']):r for r in csv.DictReader(images.open())};med=statistics.median;out={}
for sid in auth['fixed_sources']:
 entry=summary['sources'][str(sid)]; dpath=D/entry['file'];assert sha(dpath)==entry['sha256'];frozen=json.loads(dpath.read_text())['result'];controls={int(k):v for k,v in frozen['control_baselines'].items()};groups=frozen['holdout_images'];ids={int(i) for g in groups for i in g};assert ids=={int(i) for g in packet['sources'][str(sid)]['holdout'] for i in g};assert ids.isdisjoint({int(i) for g in packet['sources'][str(sid)]['development'] for i in g})
 sources=set(controls)|{sid};t=pq.read_table(parquet,columns=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings'],filters=[('source','in',list(sources)),('image_id','in',list(ids))]);assert set(t['source'].to_pylist())<=sources and set(t['image_id'].to_pylist())<=ids
 rows=t.to_pylist();seen=set();good=[];drops=[]
 for r in rows:
  k=(int(r['source']),int(r['image_id']));assert k not in seen;seen.add(k)
  try:
   v=[float(r[n]) for n in ('flux_peak','flux_peak_err','local_rms')]+[float(meta[k[1]][n]) for n in ('beam_bmaj','beam_bmin')]
   if not all(map(math.isfinite,v)) or min(v[1:])<=0:why='invalid flux/error/rms/beam'
   elif v[2]>1.82:why='local RMS exceeds 1.82 mJy/beam'
   elif r['forced'] not in (True,False):why='missing forced flag'
   elif r['flag_c4'] is True or r['has_siblings'] is True:why='fit flag'
   elif not r['forced'] and (r['flag_c4'] is not False or r['has_siblings'] is not False):why='missing detection flags'
   else:why=''
  except (TypeError,ValueError):why='malformed value'
  if why:drops.append({'source':k[0],'image_id':k[1],'reason':why})
  else:good.append(r)
 gains={};survivors={}
 for i in sorted(ids):
  rs=[r for r in good if r['image_id']==i and r['source'] in controls and not r['forced']];survivors[i]=[r['source'] for r in rs]; assert len(rs)>=20
  ratios=[float(r['flux_peak'])/controls[r['source']] for r in rs];g=med(ratios);scatter=1.4826*med([abs(q-g) for q in ratios]); frac=max(.03,1.2533*scatter/(g*math.sqrt(len(ratios))));gains[i]={'gain':g,'fractional_sigma':frac,'n':len(rs)}
 sessions=[];targetrows=[]
 for group in groups:
  rs=[r for r in good if r['source']==sid and r['image_id'] in group];assert rs;targetrows.extend(rs)
  val=sum(float(r['flux_peak'])/gains[r['image_id']]['gain'] for r in rs)/len(rs);var=max(float(r['flux_peak_err'])/gains[r['image_id']]['gain'] for r in rs)**2+(max(gains[r['image_id']]['fractional_sigma'] for r in rs)*abs(val))**2;sessions.append([val,var])
 pred=frozen['prediction'];z=[(v-pred['mean'])/math.sqrt(var+pred['variance']) for v,var in sessions]
 out[str(sid)]={'input_rows':len(rows),'frozen_control_count':len(controls),'control_survivors':survivors,'target_rows':targetrows,'drops':drops,'gains':gains,'sessions':sessions,'z':z,'both_at_least_3':all(v>=3 for v in z)}
p=O/'result.json';p.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print('result',sha(p),'python',platform.python_version())
primary=json.loads((H/'run/holdout-result.json').read_text())['sources']
for sid,x in out.items():
 y=primary[sid];assert y['status']=='measured';assert x['both_at_least_3']==y['result']['both_at_least_3'];assert sorted(x['drops'],key=lambda q:(q['source'],q['image_id']))==sorted(y['drops'],key=lambda q:(q['source'],q['image_id']))
 def compare(a,b):
  if isinstance(a,dict):assert set(map(str,a))==set(map(str,b));[compare(v,b[str(k)]) for k,v in a.items()]
  elif isinstance(a,list):assert len(a)==len(b);[compare(v,w) for v,w in zip(a,b)]
  else:assert math.isclose(a,b,abs_tol=1e-12,rel_tol=1e-12),(a,b)
 compare(x['gains'],y['holdout_gains']);compare(x['sessions'],y['holdout_sessions']);compare(x['z'],y['result']['z']);print(sid,'rows',x['input_rows'],'drops',len(x['drops']),'frozen',x['frozen_control_count'],'surviving',[(k,len(v)) for k,v in x['control_survivors'].items()],'z',x['z'],'all match')

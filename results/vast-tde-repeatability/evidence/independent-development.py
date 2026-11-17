import csv,json,math,statistics,hashlib,sys,platform
from pathlib import Path
import pyarrow.parquet as pq
P=Path(r'D:/Open-Space/worktrees/vast-tde-expanded-20261117/results/vast-tde-expanded')
packet=json.loads((P/'proposed-packet.json').read_text())
controls=list(csv.DictReader((P/'control-metadata.csv').open()))
images=Path(r'D:/Open-Space/campaigns/vast-tde-expanded-20261117/image_df.csv')
parquet=Path(r'D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_measurement_df.parquet')
assert hashlib.sha256(images.read_bytes()).hexdigest()=='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195'
assert hashlib.sha256(parquet.read_bytes()).hexdigest()==packet['measurement_parquet_sha256']
meta={int(r['image_id']):r for r in csv.DictReader(images.open())}
med=statistics.median
out={}
for sid in (14043503,13894502,14481388):
 cfg=packet['sources'][str(sid)]; groups=cfg['development']; allowed={int(i) for g in groups for i in g}; held={int(i) for g in cfg['holdout'] for i in g}; assert not allowed&held
 pool={int(r['control_source']) for r in controls if int(r['target_source'])==sid}; sources=pool|{sid}
 tab=pq.read_table(parquet,columns=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings'],filters=[('source','in',list(sources)),('image_id','in',list(allowed))])
 assert set(tab['source'].to_pylist())<=sources and set(tab['image_id'].to_pylist())<=allowed
 rows=tab.to_pylist(); keys=set(); good=[]; drops=[]
 for r in rows:
  key=(int(r['source']),int(r['image_id']));assert key not in keys;keys.add(key)
  try:
   vals=[float(r[k]) for k in ('flux_peak','flux_peak_err','local_rms')]+[float(meta[key[1]][k]) for k in ('beam_bmaj','beam_bmin')]
   if not all(map(math.isfinite,vals)) or min(vals[1:])<=0: why='invalid flux/error/rms/beam'
   elif vals[2]>1.82:why='local RMS exceeds 1.82 mJy/beam'
   elif r['forced'] not in (True,False):why='missing forced flag'
   elif r['flag_c4'] is True or r['has_siblings'] is True:why='fit flag'
   elif not r['forced'] and (r['flag_c4'] is not False or r['has_siblings'] is not False):why='missing detection flags'
   else:why=''
  except (TypeError,ValueError):why='malformed value'
  if why:drops.append({'source':key[0],'image_id':key[1],'reason':why})
  else:good.append(r)
 selected={}
 for c in pool:
  x=[]
  for g in groups:
   v=[float(r['flux_peak']) for r in good if r['source']==c and r['image_id'] in g and not r['forced'] and float(r['flux_peak'])/float(r['flux_peak_err'])>=5]
   if v:x.append(med(v))
  if len(x)>=4:
   m=med(x);mad=1.4826*med([abs(z-m) for z in x])
   if m>0 and mad/m<=.10:selected[c]=m
 result={'input_rows':len(rows),'drops':drops,'pool_count':len(pool),'selected':selected}
 try:
  if len(selected)<20:raise ValueError('too few stable development controls')
  gains={}
  for image in allowed:
   q=[float(r['flux_peak'])/selected[r['source']] for r in good if r['image_id']==image and r['source'] in selected and not r['forced']]
   if len(q)<20:raise ValueError('insufficient frozen controls for image')
   g=med(q);m=1.4826*med([abs(v-g) for v in q]);gains[image]={'gain':g,'fractional_sigma':max(.03,1.2533*m/(g*math.sqrt(len(q)))),'n':len(q)}
  sessions=[]
  for group in groups:
   rs=[r for r in good if r['source']==sid and r['image_id'] in group]
   if not rs:raise ValueError('missing source session')
   val=sum(float(r['flux_peak'])/gains[r['image_id']]['gain'] for r in rs)/len(rs)
   var=max(float(r['flux_peak_err'])/gains[r['image_id']]['gain'] for r in rs)**2+(max(gains[r['image_id']]['fractional_sigma'] for r in rs)*abs(val))**2
   sessions.append([val,var])
  w=[1/v for _,v in sessions];mean=sum(x*wi for (x,v),wi in zip(sessions,w))/sum(w);chi=sum((x-mean)**2/v for x,v in sessions);dof=len(sessions)-1
  result.update(status='qualified-for-review',gains=gains,sessions=sessions,prediction={'mean':mean,'variance':max(1,chi/dof)/sum(w),'chi2':chi,'dof':dof})
 except ValueError as e:result.update(status='inconclusive',error=str(e))
 out[str(sid)]=result
path=Path(r'D:/Open-Space/campaigns/vast-tde-independent-dev-f323ab5/result.json');path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print('output',path,'sha256',hashlib.sha256(path.read_bytes()).hexdigest(),'python',platform.python_version())
for k,v in out.items():print(k,v['status'],'rows',v['input_rows'],'drops',len(v['drops']),'pool',v['pool_count'],'selected',len(v['selected']),v.get('prediction',v.get('error')))

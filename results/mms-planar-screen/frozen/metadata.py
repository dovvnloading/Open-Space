import cdflib,numpy as np,json,hashlib
from pathlib import Path
ROOT=Path('D:/Open-Space/campaigns/heliophysics-schema-20261117/raw')
OUT=Path('D:/Open-Space/campaigns/mms-planar-metadata-20261117');OUT.mkdir(parents=True,exist_ok=True)
records={};streams=[];positions=[]
for n in range(1,5):
 files=[next(ROOT.glob(f'mms{n}_fgm_brst_l2_20260930132803_*.cdf')),next(ROOT.glob(f'mms{n}_fpi_brst_l2_dis-moms_20260930132803_*.cdf'))]
 for j,p in enumerate(files):
  c=cdflib.CDF(str(p));epoch=c.varget('Epoch');flag=c.varget(f'mms{n}_fgm_flag_brst_l2' if not j else f'mms{n}_dis_errorflags_brst');good=(flag==0);streams.append((epoch,good,'FGM' if not j else 'DIS',n))
  records[p.name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n':len(epoch),'good':int(good.sum()),'first_ns':int(epoch[0]),'last_ns':int(epoch[-1]),'largest_gap_seconds':float(np.diff(epoch).max()/1e9)}
  if not j:positions.append((c.varget('Epoch_state'),c.varget(f'mms{n}_fgm_r_gse_brst_l2')))
start=max(t[0] for t,g,_,_ in streams);end=min(t[-1] for t,g,_,_ in streams);t=np.arange(start,end,int(1e8),dtype=np.int64);accepted=np.ones(len(t),dtype=bool);details=[]
for epoch,good,typ,n in streams:
 ix=np.searchsorted(epoch,t);right=np.clip(ix,0,len(epoch)-1);left=np.clip(ix-1,0,len(epoch)-1);pick=np.where(abs(epoch[left]-t)<=abs(epoch[right]-t),left,right);dist=np.abs(epoch[pick]-t)/1e9;limit=.005 if typ=='FGM' else .08;pass_=good[pick]&(dist<=limit);accepted &= pass_;details.append({'craft':n,'product':typ,'grid_pass':int(pass_.sum()),'tolerance_s':limit})
segments=[];ind=np.flatnonzero(accepted)
for chunk in np.split(ind,np.where(np.diff(ind)>1)[0]+1):
 if len(chunk):segments.append({'start_ns':int(t[chunk[0]]),'end_ns':int(t[chunk[-1]]),'seconds_grid_extent':round((len(chunk)-1)*.1,2)})
mid=int((start+end)//2);coords=[]
for tt,rr in positions:
 assert tt[0]<=mid<=tt[-1]; coords.append([float(np.interp(0,tt.astype(float)-mid,rr[:,k])) for k in range(3)])
coords=np.array(coords);edges=[float(np.linalg.norm(coords[i]-coords[j])) for i in range(4) for j in range(i+1,4)];matrix=(coords[1:]-coords[0]).T;v=abs(float(np.linalg.det(matrix)))/6
out={'records':records,'common_span_utc':[str(cdflib.cdfepoch.to_datetime(start)),str(cdflib.cdfepoch.to_datetime(end))],'grid_0p1s_points':len(t),'joint_good_grid_count':int(accepted.sum()),'streams':details,'joint_segments':segments,'position_mid_utc':str(cdflib.cdfepoch.to_datetime(mid)),'position_gse_km':coords.tolist(),'edge_km':edges,'tetra_volume_km3':v,'condition_number':float(np.linalg.cond(matrix)),'qualification':'nearest sample grid check only, not proof of continuous good support; position interpolation is metadata'}
(OUT/'metadata.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('records','joint_segments','position_gse_km')},indent=2));print('segments',len(segments),'longest_s',max(z['seconds_grid_extent'] for z in segments))

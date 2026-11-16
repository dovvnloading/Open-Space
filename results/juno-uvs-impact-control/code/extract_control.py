"""Bounded HTTP range acquisition of a published control; not a blind search."""
import hashlib,json,pathlib,subprocess,tempfile
import numpy as np
from astropy.time import Time
URL='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-26/UVS_S01_639769969_2020101_P26OBS_V01.FIT'
ROOT=pathlib.Path(__file__).resolve().parents[1]; CACHE=ROOT/'range-cache';CACHE.mkdir(exist_ok=True)
N=14476358; OFFSET=4320000; WIDTH=84
names=['hack','x','y','angle','wave','et','spin','pulse','eclat','eclon','weight','lat','lon','tangent','incidence','emission','digital','analog','local']
dtype=np.dtype(list(zip(names,['>i4']*3+['>f4']*2+['>f8','>f4','>i4']+['>f4']*2+['>f8']+['>f4']*8)))
receipts=[]
def rows(a,b):
 start=OFFSET+a*WIDTH; end=OFFSET+b*WIDTH-1;p=CACHE/f'{start}-{end}.bin';head=CACHE/f'{start}-{end}.http'
 if not p.exists():
  subprocess.run(['curl','--fail','-Ls','--retry','2','--max-time','45','-r',f'{start}-{end}','-D',str(head),URL,'-o',str(p)],check=True)
 data=p.read_bytes();assert len(data)==(b-a)*WIDTH
 assert f'bytes {start}-{end}/1246985280' in head.read_text()
 receipts.append({'start':start,'end':end,'sha256':hashlib.sha256(data).hexdigest()})
 return np.frombuffer(data,dtype=dtype)
target=(Time('2020-04-10T12:57:10',scale='utc').tdb.jd-2451545)*86400
# Interpolation-assisted lower-bound search; exact final bound requires sorted times.
def lower(t):
 lo,hi=0,N-1;tl=float(rows(lo,lo+1)['et'][0]);th=float(rows(hi,hi+1)['et'][0])
 for _ in range(40):
  if hi-lo<2:return hi
  mid=max(lo+1,min(hi-1,int(lo+(hi-lo)*(t-tl)/(th-tl))))
  if hi-lo<2000: mid=(lo+hi)//2
  tm=float(rows(mid,mid+1)['et'][0]);print('probe',mid,tm-target,flush=True)
  if tm<t:lo,tl=mid,tm
  else:hi,th=mid,tm
 raise RuntimeError('search failed')
a=lower(target-2);b=lower(target+2);assert (b-a)*WIDTH<20_000_000
x=rows(a,b);assert np.all(np.diff(x['et'])>=0)
np.savez_compressed(ROOT/'control-photons.npz',photons=x)
json.dump({'url':URL,'target_et':target,'first_row':a,'last_row_exclusive':b,'records':len(x),'receipts':receipts},open(ROOT/'acquisition.json','w'),indent=2)
print('DONE',a,b,len(x),flush=True)

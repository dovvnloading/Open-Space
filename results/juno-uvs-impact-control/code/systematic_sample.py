"""Frozen HK-selected sample; native FITS layout and range-verified extraction."""
import pathlib,subprocess,re,math,io,json,hashlib,concurrent.futures
import threading
locks={}
import numpy as np
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1];out=p/'systematic';out.mkdir(exist_ok=True)
base='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/'
runs={2:[(145,148),(428,435),(729,739)],3:[(3464,3470),(4257,4262),(4646,4655),(5066,5076),(5364,5367)]}
def get(name,a,b):
 with locks.setdefault((name,a,b),threading.Lock()):return _get(name,a,b)
def _get(name,a,b):
 f=out/f'{name[:7]}-{a}-{b}.bin';http=f.with_suffix('.http')
 if not f.exists() or f.stat().st_size!=b-a+1:subprocess.run(['curl','--fail','-Ls','--retry','1','--max-time','90','-r',f'{a}-{b}','-D',str(http),base+name,'-o',str(f)],check=True)
 z=f.read_bytes();assert len(z)==b-a+1;assert f'bytes {a}-{b}/' in http.read_text();return z

def header(z,start=0):
 e=next(k+80 for k in range(start,len(z),80) if z[k:k+8]==b'END     ');hs=math.ceil((e-start)/2880)*2880;h=fits.Header.fromstring(z[start:start+hs].decode(),sep='');n=abs(h['BITPIX'])//8
 for j in range(h['NAXIS']):n*=h[f'NAXIS{j+1}']
 n+=h.get('PCOUNT',0);return h,hs,n

def table(z):return fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+z))[1].data

def run(s):
 name=f'UVS_S{s:02}_805506529_2025192_P74OBS_V01.FIT';head=get(name,0,4319999);o=0
 for j in range(2):
  h,hs,n=header(head,o)
  if j==1:frames=table(head[o:o+hs+math.ceil(n/2880)*2880]);np.savez_compressed(out/f's{s}-frames.npz',**{k:frames[k] for k in frames.names})
  o+=hs+math.ceil(n/2880)*2880
 # photon header may fit at end of initial range; otherwise fetch
 ph=get(name,o,o+287999);h,hs,n=header(ph);offset=o+hs;N=h['NAXIS2'];width=h['NAXIS1'];types={'J':'>i4','E':'>f4','D':'>f8'};dt=np.dtype([(h[f'TTYPE{i}'],types[h[f'TFORM{i}'].strip()]) for i in range(1,h['TFIELDS']+1)]);assert dt.itemsize==width
 (out/f's{s}-photon-header.txt').write_text(h.tostring())
 def rows(a,b):return np.frombuffer(get(name,offset+width*a,offset+width*b-1),dtype=dt)
 first,last=float(rows(0,1)['EPHEMERIS_TIME'][0]),float(rows(N-1,N)['EPHEMERIS_TIME'][0])
 def locate(t):
  lo,hi,tl,th=0,N-1,first,last
  for k in range(24):
   guess=(lo+hi)//2;a=max(lo,guess-4096);b=min(hi+1,guess+4096);x=rows(a,b);q=x['EPHEMERIS_TIME'];assert np.all(np.diff(q)>=0)
   if q[0]<=t<=q[-1]:return a+int(np.searchsorted(q,t))
   if q[-1]<t:lo,tl=b-1,float(q[-1])
   else:hi,th=a,float(q[0])
  raise RuntimeError('timestamp bound not found')
 hk=np.load(p/'clean-screen'/(name+'.hk.npz'));result=[]
 times=[float(hk['EPHEMERIS_TIME'][k]) for r in runs[s] for k in r]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: located=dict(zip(times,ex.map(locate,times)))
 for a,b in runs[s]:
  bounds=[float(hk['EPHEMERIS_TIME'][a]),float(hk['EPHEMERIS_TIME'][b])]
  i,j=[located[v] for v in bounds]
  assert (j-i)*width<100_000_000;x=rows(i,j);assert np.all(np.diff(x['EPHEMERIS_TIME'])>=0);assert np.all((x['EPHEMERIS_TIME']>=bounds[0])&(x['EPHEMERIS_TIME']<bounds[1]));np.savez_compressed(out/f's{s}-{a}-{b}.npz',photons=x)
  r={'file':name,'hk_indices':[a,b],'et_bounds':bounds,'nominal_seconds':bounds[1]-bounds[0],'rows':[i,j],'photon_count':len(x),'data_offset':offset,'row_width':width,'sha256':hashlib.sha256(x.tobytes()).hexdigest()};result.append(r);(out/f's{s}-acquisition.json').write_text(json.dumps(result,indent=2));print(json.dumps(r),flush=True)
 # dynamic ancillary table retrieval
 ao=offset+math.ceil(n/2880)*2880;az=get(name,ao,ao+287999);ah,ahs,an=header(az);assert ah['EXTNAME']=='Ancillary Data';az=get(name,ao,ao+ahs+math.ceil(an/2880)*2880-1);anc=table(az);np.savez_compressed(out/f's{s}-ancillary.npz',**{k:anc[k] for k in anc.names})
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(run,[2,3]))

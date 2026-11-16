import numpy as np,pathlib,subprocess,json,hashlib,concurrent.futures,threading
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1]/'clean-screen';name='UVS_S02_805506529_2025192_P74OBS_V01.FIT';url='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/'+name
h=fits.Header.fromstring((p/'s02-photon-header.txt').read_text(),sep='');types={'J':'>i4','E':'>f4','D':'>f8'};dt=np.dtype([(h[f'TTYPE{i}'],types[h[f'TFORM{i}'].strip()]) for i in range(1,20)]);N=h['NAXIS2'];receipts=[];lock=threading.Lock()
def rows(a,b):
 s=4320000+84*a;e=4320000+84*b-1;f=p/f'range-{s}-{e}.bin';hdr=p/f'range-{s}-{e}.http'
 if not f.exists():subprocess.run(['curl','--fail','-Ls','--max-time','60','--retry','1','-r',f'{s}-{e}','-D',str(hdr),url,'-o',str(f)],check=True)
 z=f.read_bytes();assert len(z)==(b-a)*84;assert f'bytes {s}-{e}/' in hdr.read_text()
 with lock:receipts.append({'start':s,'end':e,'sha256':hashlib.sha256(z).hexdigest()})
 return np.frombuffer(z,dtype=dt)

hk=np.load(p/(name+'.hk.npz'));bounds=[float(hk['EPHEMERIS_TIME'][729]),float(hk['EPHEMERIS_TIME'][739])]
# Frozen timestamp bounds from the metadata-selected interval.
a,b=12431751,12441538
assert (b-a)*84<100_000_000;x=rows(a,b);np.savez_compressed(p/'clean-photons.npz',photons=x);json.dump({'source':url,'first_row':a,'last_row_exclusive':b,'et_bounds':bounds,'rows':len(x),'receipts':receipts},open(p/'clean-acquisition.json','w'),indent=2);print('DONE',a,b,len(x),flush=True)

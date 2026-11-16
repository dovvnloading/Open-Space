import pathlib,json,time,subprocess,io,hashlib
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1]/'newer-pj74';u='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/UVS_S01_805506529_2025192_P74OBS_V01.FIT'
for i in [1,6,8]:
 for _ in range(150):
  infos=json.load(open(p/'hdu-layout.json'))
  if len(infos)>i:break
  time.sleep(2)
 a=infos[i];f=p/f'aux-{i}.bin';header=(p/f'header-{i}.txt').read_text();start=a['data_offset'];end=start+a['data_bytes']-1
 subprocess.run(['curl','--fail','-Ls','--retry','1','--max-time','40','-r',f'{start}-{end}',u,'-o',str(f)],check=True);b=f.read_bytes();assert len(b)==a['data_bytes'];raw=header.encode()+b+b'\0'*((-len(b))%2880)
 hdu=fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+raw))[1]
 print('HDU',i,repr(hdu.columns),flush=True)
 npz={name:hdu.data[name] for name in hdu.columns.names};import numpy as np;np.savez_compressed(p/f'aux-{i}.npz',**npz)
 json.dump({'source':u,'start':start,'end':end,'sha256':hashlib.sha256(b).hexdigest()},open(p/f'aux-{i}-receipt.json','w'),indent=2)

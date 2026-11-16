"""Replay pinned auxiliary ranges; derive and validate FITS table shapes."""
import pathlib,subprocess,io,math,re,hashlib,json
import numpy as np
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1]/'clean-screen';name='UVS_S02_805506529_2025192_P74OBS_V01.FIT';url='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/'+name;receipts=[]
for nm,start,end in [('frames',4210560,4314239),('ancillary',1258608960,1261608959)]:
 f=p/f's02-{nm}.bin';http=p/f's02-{nm}.http';subprocess.run(['curl','--fail','-Ls','--retry','1','--max-time','60','-r',f'{start}-{end}','-D',str(http),url,'-o',str(f)],check=True);b=f.read_bytes();assert len(b)==end-start+1;assert f'bytes {start}-{end}/' in http.read_text();e=next(k+80 for k in range(0,len(b)-80,80) if b[k:k+8]==b'END     ');hs=math.ceil(e/2880)*2880;h=fits.Header.fromstring(b[:hs].decode(),sep='');size=h['NAXIS1']*h['NAXIS2'];raw=b[:hs+math.ceil(size/2880)*2880];d=fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+raw))[1].data;np.savez_compressed(p/f's02-{nm}.npz',**{k:d[k] for k in d.names});receipts.append({'source':url,'start':start,'end':end,'sha256':hashlib.sha256(b).hexdigest()})
b=(p/(name+'.tail')).read_bytes();found=False
for match in re.finditer(b'XTENSION=',b):
 j=match.start()
 try:
  e=next(k+80 for k in range(j,len(b)-80,80) if b[k:k+8]==b'END     ');hs=math.ceil((e-j)/2880)*2880;h=fits.Header.fromstring(b[j:j+hs].decode(),sep='')
  if h.get('EXTNAME')!='Mask Information':continue
  d=fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+b[j:]))[1].data;np.savez_compressed(p/'s02-masks.npz',**{k:d[k] for k in d.names});found=True
 except (UnicodeDecodeError,StopIteration):continue
assert found
(p/'clean-aux-receipts.json').write_text(json.dumps(receipts,indent=2))

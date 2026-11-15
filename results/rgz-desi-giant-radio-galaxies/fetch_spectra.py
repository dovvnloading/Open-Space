"""Read the public DESI DR1 coadd by FITS byte ranges for two target IDs."""
import urllib.request as u
import numpy as np,concurrent.futures,hashlib,io,pathlib,csv,time
def load(args):
 id,p=args
 url=f'https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/main/dark/{p//100}/{p}/coadd-main-dark-{p}.fits'
 def get(a,n):
  for i in range(3):
   try:
    req=u.Request(url,headers={'Range':f'bytes={a}-{a+n-1}'})
    with u.urlopen(req,timeout=50) as r:
     if not r.headers.get('Content-Range','').startswith(f'bytes {a}-'):raise ValueError(str(r.headers.get('Content-Range')))
     d=r.read(n)
    if len(d)!=n:raise ValueError('short')
    return d
   except Exception:
    if i==2:raise
 def header(pos):
  h={};offset=0
  while True:
   d=get(pos+offset,2880);offset+=len(d)
   for i in range(0,2880,80):
    c=d[i:i+80].decode()
    if c[:8].strip()=='END':return h,offset
    if c[8:10]=='= ':h[c[:8].strip()]=c[10:].split('/')[0].strip().strip("'").strip()
 off=0;row=None;out={}
 while 'Z_MASK' not in out:
  h,skip=header(off);off+=skip;name=h.get('EXTNAME','P');axis=int(h['NAXIS']);n=0 if axis==0 else int(np.prod([int(h[f'NAXIS{i+1}']) for i in range(axis)]))*abs(int(h['BITPIX']))//8+int(h.get('PCOUNT',0))
  if name=='FIBERMAP':
   assert h['TTYPE1']=='TARGETID' and h['TFORM1']=='K'
   width=int(h['NAXIS1']);d=get(off,n);ids=np.ndarray(shape=(int(h['NAXIS2']),),dtype='>i8',buffer=d,strides=(width,));idx=np.where(ids==id)[0];assert len(idx)==1;row=int(idx[0])
  if '_' in name and name.split('_',1)[-1] in ['WAVELENGTH','FLUX','IVAR','MASK'] and name[0] in 'BRZ':
   typ={'-64':'>f8','-32':'>f4','32':'>i4'}[h['BITPIX']];width=int(h['NAXIS1'])*np.dtype(typ).itemsize;a=off+(0 if axis==1 else row*width);out[name]=np.frombuffer(get(a,width),dtype=typ)
  off+=((n+2879)//2880)*2880
  print(id,name,flush=True)
 rows=[]
 for b in 'BRZ':
  for v in zip(*[out[b+'_'+k] for k in ['WAVELENGTH','FLUX','IVAR','MASK']]): rows.append((b,*v))
 with open(str(id)+'_spec.csv','w') as st:
  w=csv.writer(st);w.writerow(['band','wavelength','flux','ivar','mask']);w.writerows(rows)
if __name__=='__main__':
 list(concurrent.futures.ThreadPoolExecutor(2).map(load,[(39628055064347506,19649),(39628405439727850,4915)]))

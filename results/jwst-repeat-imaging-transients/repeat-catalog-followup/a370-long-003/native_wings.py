"""Independent aperture measurements in public Level-2 exposures at a370_long-003."""
import argparse, concurrent.futures,hashlib
from io import BytesIO
import time
from urllib.request import urlopen,Request
from pathlib import Path
import numpy as np,pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
RA=39.96591521;DEC=-1.59954268
def fetch(url,lo,hi):
 for i in range(3):
  try:
   with urlopen(Request(url,headers={'Range':f'bytes={lo}-{hi}'}),timeout=60) as r:
    if r.status!=206 or not r.getheader('Content-Range','').startswith(f'bytes {lo}-{hi}/'):raise ValueError('range response')
    v=r.read(); assert len(v)==hi-lo+1;return v
  except Exception:
   if i==2:raise
   time.sleep(i+1)
def head(url,off):
 v=fetch(url,off,off+57599); f=BytesIO(v);h=fits.Header.fromfile(f);return h,off+f.tell()
def next_h(h,off):
 n=np.prod([h['NAXIS'+str(i+1)] for i in range(h['NAXIS'])])*abs(h['BITPIX'])//8
 return off+int(np.ceil(n/2880))*2880
def one(ob):
 p=ob.obs_id;root=p.split('_')[0];url=f'https://stpubdata.s3.amazonaws.com/jwst/public/{root[:7]}/{root}/{p}_cal.fits'
 rc=dict(product=p,url=url,archive_mjd=ob.t_min)
 try:
  pre=fetch(url,0,115199);fh=BytesIO(pre);pr=fits.Header.fromfile(fh);s=fits.Header.fromfile(fh);off=fh.tell()
  w=WCS(s);x,y=w.all_world2pix(RA,DEC,0,quiet=True);x=float(x);y=float(y)
  rc.update(mjd=pr['EXPMID'],filter=pr['FILTER'],header_sha256=hashlib.sha256(pre).hexdigest())
  dim=s['NAXIS1'];ny=s['NAXIS2'];scale=np.sqrt(s['PIXAR_A2']);r=int(np.ceil(.8/scale));x0=round(x)-r;y0=round(y)-r
  rc['edge_arcsec']=min(x,y,dim-1-x,ny-1-y)*scale
  if min(x0,y0)<0 or x0+2*r+1>dim or y0+2*r+1>ny:return rc
  def patch(h,start):
   n0=y0*dim*4;nb=(2*r+1)*dim*4;b=fetch(url,start+n0,start+n0+nb-1)
   ar=np.frombuffer(b,dtype='>f4' if h['BITPIX']==-32 else '>i4').reshape(2*r+1,dim)[:,x0:x0+2*r+1]
   if h.get('BZERO'):ar=ar.astype(np.int64)+int(h['BZERO'])
   rc[h['EXTNAME']+'_sha256']=hashlib.sha256(b).hexdigest()
   return ar
  im=patch(s,off); prevh,prevstart=s,off; get={}
  for j in range(5):
   h,start=head(url,next_h(prevh,prevstart))
   if h['EXTNAME'] in ['ERR','DQ']:get[h['EXTNAME']]=patch(h,start)
   if len(get)==2:break
   prevh,prevstart=h,start
  err=get['ERR'];dq=get['DQ'];unit=s['PIXAR_SR']*1e12
  yy,xx=np.mgrid[y0:y0+2*r+1,x0:x0+2*r+1];radius=scale*np.hypot(xx-x,yy-y)
  ok=np.isfinite(im)&np.isfinite(err)&(err>0)&((dq&3)==0);ann=(radius>=.45)&(radius<.65)&ok
  if ann.sum()<50:return rc
  sky=np.median(im[ann])
  for ri in (12,18,25,35):
   wg=np.zeros_like(im,dtype=float)
   for sx in (-.375,-.125,.125,.375):
    for sy in (-.375,-.125,.125,.375): wg+=scale*np.hypot(xx+sx-x,yy+sy-y)<=ri/100
   wg/=16; core=wg>0;valid=core&ok; ar=wg[valid].sum()/wg.sum()
   rc[f'coverage_{ri}']=ar;rc[f'dq_{ri}']=int(np.bitwise_or.reduce(dq[core]))
   if ar>=.98:
    rc[f'flux_{ri}_ujy']=unit*np.sum(wg[valid]*(im[valid]-sky))/ar
    rc[f'err_{ri}_ujy']=unit*np.sqrt(np.sum((wg[valid]*err[valid])**2))/ar
  for inner,outer in ((0,.07),(.07,.12),(.12,.18),(.18,.30),(.30,.45)):
   inside=(radius>=inner)&(radius<outer);k=f'{round(inner*100):02d}_{round(outer*100):02d}'
   rc['bad_'+k]=int(np.count_nonzero(inside&(~ok)))
   rc['saturated_'+k]=int(np.count_nonzero(inside&((dq&2)>0)))
   rc['npix_'+k]=int(inside.sum())
   # A missing bright core pixel cannot be reconstructed by area scaling.
   vv=inside&ok;rc['annular_ujy_'+k]=unit*(im[vv]-sky).sum() if (inside&(~ok)).sum()==0 else np.nan
   if (inside&(~ok)).sum()==0: rc['annular_err_ujy_'+k]=unit*np.sqrt((err[vv]**2).sum())
  return rc
 except Exception as e:rc['error']=str(e);return rc
if __name__=='__main__':
 a=pd.read_csv('l2.csv')
 # Six primary wide bands in the oldest, January and August visits, plus the December 2024/January 2026 independent programs.
 q=a[a.obs_id.str.startswith(('jw01208','jw05324014','jw05890','jw06882'))]
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
  results=[]
  for r in pool.map(one,q.itertuples()):
   results.append(r);pd.DataFrame(results).to_csv('native_wings.csv',index=False,float_format='%.10g')

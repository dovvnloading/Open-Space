"""Fetch frozen candidate spectra; inspect identity and exposure lists only, not flux."""
import csv,hashlib,json,sys,urllib.request,concurrent.futures
from pathlib import Path
from astropy.io import fits
import numpy as np
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);rows=list(csv.DictReader(open(sys.argv[2],newline='')))[:32]
def get(job):
 r,kind=job;sid=int(r['sdss_id']);m=int(r[kind+'_mjd']);cid=r[kind+'_catalogid'];f=f'spec-112359-{m}-{cid}.fits';u=f'https://dr20.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/spectra/daily/lite/112XXX/112359/{m}/{f}';p=out/f
 try:
  data=urllib.request.urlopen(u,timeout=60).read();p.write_bytes(data)
  with fits.open(p,memmap=False) as h:
   s=h['SPALL'].data[0]; meta=dict(sdss_id=int(s['SDSS_ID']),ra=float(s['RACAT']),dec=float(s['DECCAT']),mjd_list=str(s['MJD_LIST']).split(),tai_list=str(s['TAI_LIST']).split(),z=float(s['Z']),warning=int(s['ZWARNING']),filename=f,kind=kind,url=u,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
  if meta['sdss_id']!=sid:raise ValueError('SDSS_ID mismatch')
  return meta
 except Exception as e:return dict(sdss_id=sid,kind=kind,url=u,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:result=list(ex.map(get,[(r,k) for r in rows for k in ['baseline','dev','hold']]))
Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n');print('retrieved',sum('sha256' in x for x in result),'of',len(result),'errors',sum('error'in x for x in result))

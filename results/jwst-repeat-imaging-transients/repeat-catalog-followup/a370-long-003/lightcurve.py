"""Measure the nearest Level-3 catalog entry at a370_long-003 in every covering NIRCam product."""
import concurrent.futures
from io import BytesIO
from pathlib import Path
import sys,hashlib,time
from urllib.request import urlopen
import numpy as np,pandas as pd
RA=39.96591521;DEC=-1.59954268
def product_url(p):
 prog,obs=p.split('_')[0].split('-');t=p.split('_')[1][0]
 return f'https://stpubdata.s3.amazonaws.com/jwst/public/{prog}/L3/{t}/{obs}/{p}_cat.ecsv'
def one(ob):
 url=product_url(ob.obs_id);rec=dict(obs_id=ob.obs_id,mjd_start=ob.t_min,mjd_end=ob.t_max,filter=ob.obs_id.split('-')[-1],url=url)
 for attempt in range(3):
  try:
   with urlopen(url,timeout=100) as r:raw=r.read()
   a=pd.read_csv(BytesIO(raw),sep=r'\s+',comment='#',usecols=lambda s:s in ['label','sky_centroid.ra','sky_centroid.dec','aper_total_flux','aper_total_flux_err','aper50_flux','aper70_flux','is_extended','nn_dist','xcentroid','ycentroid','semimajor_sigma','semiminor_sigma','fwhm'])
   dx=(a['sky_centroid.ra']-RA)*3600*np.cos(np.radians(DEC));dy=(a['sky_centroid.dec']-DEC)*3600
   ds=np.hypot(dx,dy).dropna().sort_values();j=ds.index[0];r=a.loc[j]
   rec.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),catalog_sources=len(a),separation_arcsec=ds.iloc[0],second_sep=ds.iloc[1],source_id=r.label,ra=r['sky_centroid.ra'],dec=r['sky_centroid.dec'],matched=ds.iloc[0]<.15)
   if rec['matched']:
    rec.update(flux_ujy=r.aper_total_flux*1e6,error_ujy=r.aper_total_flux_err*1e6,aper50_ujy=r.aper50_flux*1e6,aper70_ujy=r.aper70_flux*1e6,extended=r.is_extended,fwhm_px=r.get('fwhm',np.nan))
   return rec
  except Exception as e:
   if attempt==2:rec['error']=str(e);return rec
   time.sleep(1+attempt)
if __name__=='__main__':
 obs=pd.read_csv('l3.csv')
 obs['public_at_2026_11_16']=pd.to_datetime(obs.datarelease)<pd.Timestamp('2026-11-17T00:00')
 with concurrent.futures.ThreadPoolExecutor(max_workers=7) as pool:out=list(pool.map(one,obs[obs.public_at_2026_11_16].itertuples()))
 obs.merge(pd.DataFrame(out),how='left',on='obs_id').to_csv('all_epoch_catalogs.csv',index=False,float_format='%.10g')

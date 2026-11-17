import pathlib,csv,numpy as np,pyarrow.parquet as pq
from astropy.coordinates import SkyCoord,search_around_sky
from astropy import units as u
p=pathlib.Path('D:/Open-Space/campaigns/vast-tde-expanded-20261117');rows=list(csv.DictReader((p/'TDE_catalogue_all.csv').open(encoding='utf-8-sig',newline='')));raw=pq.read_table('D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_source_df.parquet',columns=['source','wavg_ra','wavg_dec','wavg_uncertainty_ew','wavg_uncertainty_ns','n_measurements','n_forced']).to_pydict();radio=SkyCoord(np.array(raw['wavg_ra'])*u.deg,np.array(raw['wavg_dec'])*u.deg);results=[]
for i,r in enumerate(rows):
 name=r['AT name'].strip();date=r['Discovery date (UT)'].strip();status='eligible pre-2023-06-01' if date and date[:10]<'2023-06-01' else 'later or missing discovery date';base=dict(row=i+1,name=name,ztf=r['ZTF name'],alternative=r['Alternative name'],discovery=date,ra=r['RA'],dec=r['DEC'],status=status,source_id='',separation_arcsec='',second_nearest_arcsec='',uncertainty_ew='',uncertainty_ns='',n_measurements='',n_forced='')
 if status.startswith('eligible'):
  try:
   coord=SkyCoord(r['RA'],r['DEC'],unit=(u.hourangle,u.deg));ix,sep,_=coord.match_to_catalog_sky(radio);ix=int(np.asarray(ix).ravel()[0]);sep=float(np.asarray(sep.arcsec).ravel()[0]);_,sep2,_=coord.match_to_catalog_sky(radio,nthneighbor=2);base['second_nearest_arcsec']=float(np.asarray(sep2.arcsec).ravel()[0]);base['separation_arcsec']=sep
   if sep<=30:base.update(source_id=raw['source'][ix],uncertainty_ew=raw['wavg_uncertainty_ew'][ix],uncertainty_ns=raw['wavg_uncertainty_ns'][ix],n_measurements=raw['n_measurements'][ix],n_forced=raw['n_forced'][ix],status='candidate association; unverified')
   else:base['status']='no catalog counterpart within 30 arcsec'
  except Exception as exc:base['status']='coordinate error: '+str(exc)
 results.append(base)
with (p/'expanded-metadata.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
from collections import Counter
print('total',len(results),Counter(r['status'] for r in results));print(*[r for r in results if r['source_id']],sep='\n')

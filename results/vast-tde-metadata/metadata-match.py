import pathlib,csv,io,hashlib,json,pyarrow.parquet as pq,numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
p=pathlib.Path('D:/Open-Space/campaigns/vast-tde-metadata-20261117');lines=[x for x in (p/'ztf-table1.tsv').read_text(encoding='utf-8').splitlines() if x.strip() and not x.startswith('#')]; rows=list(csv.DictReader(io.StringIO('\n'.join([lines[0]]+lines[3:])),delimiter='\t'));rows=[{k:v.strip() for k,v in r.items()} for r in rows if r['recno'].strip().isdigit()];assert len(rows)==30
f=p/'VAST_Extragal_DR1_source_df.parquet';tab=pq.read_table(f,columns=['source','wavg_ra','wavg_dec','wavg_uncertainty_ew','wavg_uncertainty_ns','n_measurements','n_forced']);d=tab.to_pydict();coords=SkyCoord(np.asarray(d['wavg_ra'])*u.deg,np.asarray(d['wavg_dec'])*u.deg);out=[]
for r in rows:
 row=dict(name=r['Name'],ztf=r['ZTF'],spectroscopic_class=r['Spec'],redshift=r['zspec'],coordinate_source='VizieR SIMBAD-added coordinates',ra_deg=r['_RA'],dec_deg=r['_DE'],disposition='no radio match within 30 arcsec')
 if r['_RA'] and r['_DE']:
  c=SkyCoord(float(r['_RA'])*u.deg,float(r['_DE'])*u.deg);ix,sep,_=c.match_to_catalog_sky(coords);ix=int(np.asarray(ix).ravel()[0]);separation=float(np.asarray(sep.arcsec).ravel()[0]);row.update(nearest_arcsec=separation)
  if separation<=30:row.update(vast_source=int(d['source'][ix]),n_measurements=int(d['n_measurements'][ix]),n_forced=int(d['n_forced'][ix]),disposition='candidate match; association unverified')
 else:row['disposition']='missing position'
 out.append(row)
with (p/'cohort-metadata.csv').open('w',newline='') as f:
 fields=['name','ztf','spectroscopic_class','redshift','coordinate_source','ra_deg','dec_deg','nearest_arcsec','vast_source','n_measurements','n_forced','disposition'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
print('30 positions; within30',sum('vast_source' in r for r in out));print('\n'.join(str(r) for r in out if 'vast_source' in r))

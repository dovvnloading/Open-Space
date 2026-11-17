import csv,pathlib,numpy as np,pyarrow.parquet as pq
from astropy.coordinates import SkyCoord,search_around_sky
from astropy import units as u
p=pathlib.Path('D:/Open-Space/campaigns/vast-tde-expanded-20261117');source='D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_source_df.parquet';measurement='D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_measurement_df.parquet';d=pq.read_table(source,columns=['source','wavg_ra','wavg_dec']).to_pydict(); coords=SkyCoord(np.asarray(d['wavg_ra'])*u.deg,np.asarray(d['wavg_dec'])*u.deg); parent=list(csv.DictReader((p/'TDE_catalogue_all.csv').open(encoding='utf-8-sig')));parents=SkyCoord([r['RA'] for r in parent],[r['DEC'] for r in parent],unit=(u.hourangle,u.deg));names=(p/'vast-table2-names.txt').read_text(encoding='utf-8').splitlines();import re
pat=re.compile(r'ASKAP J(\d{6}(?:\.\d)?)([+−-])(\d{6})');parsed=[pat.fullmatch(x).groups() for x in names];published=SkyCoord([f'{h[:2]}:{h[2:4]}:{h[4:]}' for h,s,z in parsed],[f'{"-" if s in "−-" else "+"}{z[:2]}:{z[2:4]}:{z[4:]}' for h,s,z in parsed],unit=(u.hourangle,u.deg));e=list(csv.DictReader((p/'epoch-metadata.csv').open()));out=[]
for sid in [14043503,13894502,14481388]:
 ix=d['source'].index(sid);target=coords[ix];sep=target.separation(coords).deg;inds=np.where((sep<=1)&(sep>60/3600))[0];near=coords[inds];_,psep,_=near.match_to_catalog_sky(parents);_,vsep,_=near.match_to_catalog_sky(published);inds=inds[(psep.arcsec>60)&(vsep.arcsec>30)]; candidate_ids=[int(d['source'][i]) for i in inds];targets=[x for x in e if int(x['source'])==sid];target_images={int(x['image_id']) for x in targets};t=pq.read_table(measurement,columns=['source','image_id'],filters=[('source','in',candidate_ids)]).to_pydict();coverage={}
 for a,b in zip(t['source'],t['image_id']):
  if b in target_images:coverage.setdefault(a,set()).add(b)
 candidates=sorted((float(sep[i]),int(d['source'][i]),len(coverage.get(int(d['source'][i]),set()))) for i in inds if len(coverage.get(int(d['source'][i]),set()))>=4)[:200]
 for rank,(dist,control,n) in enumerate(candidates,1):out.append(dict(target_source=sid,rank=rank,control_source=control,separation_arcsec=round(dist*3600,4),common_images=n))
 print(sid,'eligible controls',len(candidates))
with (p/'control-metadata.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['target_source','rank','control_source','separation_arcsec','common_images']);w.writeheader();w.writerows(out)

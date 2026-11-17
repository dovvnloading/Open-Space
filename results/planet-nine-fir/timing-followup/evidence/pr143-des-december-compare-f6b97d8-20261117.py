import csv,gzip,json,collections,hashlib
from pathlib import Path
base=Path('D:/Open-Space/campaigns'); paths={'published':Path('results/planet-nine-fir/archives/planet_draws.csv.gz'),'historical':base/'pr143-des-historical-581c9c4-ea2d60a-20261117/run/planet_draws.csv.gz','dec1':base/'pr143-des-dec1-f6b97d8-296a210-20261117/run/planet_draws.csv.gz','dec2':base/'pr143-des-dec2-f6b97d8-b000d53-20261117/run/planet_draws.csv.gz'}
rows={name:{(r['experiment'],r['draw']):r for r in csv.DictReader(gzip.open(path,'rt'))} for name,path in paths.items()};flags=['hit_mass5_R2.71_p0.10','hit_mass5_R2.94_p0.30','hit_mass10_R3.70_p0.30'];ans={}
for phase in ('dec1','dec2'):
 p=paths[phase].parent;o=json.loads((p/'draw_outcomes.json').read_text());z=json.loads((p/'difference_result.json').read_text());ans[phase]={'status':json.loads((p/'run_status.json').read_text()),'coverage':{'/'.join(k):v for k,v in collections.Counter((r['experiment'],r['state_status'],r['image_status']) for r in o).items()},'flags':{},'dense_peaks':json.loads((p/'dense_peaks.json').read_text()),'difference':{k:z[k] for k in ('dense_max_snr','dense_area_pixels','dense_eligible_pixels','exposure_ids','mjd_mid','orbit_grid','mask_valid','max_snr','n_above','threshold')},'hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.iterdir() if f.is_file()}}
 for exp in ('published_centroids','centroid_stress'):
  for flag in flags:
   keys=[k for k in rows[phase] if k[0]==exp];ans[phase]['flags'][exp+'/'+flag]={'count':sum(int(rows[phase][k][flag]) for k in keys),'vs_historical':[k[1] for k in keys if rows[phase][k][flag]!=rows['historical'][k][flag]],'vs_published':[k[1] for k in keys if rows[phase][k][flag]!=rows['published'][k][flag]]}
p=base/'pr143-des-december-comparison-f6b97d8-20261117.json';p.write_text(json.dumps(ans,indent=2)+'\n');print(json.dumps(ans,indent=2))

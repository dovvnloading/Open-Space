import csv,gzip,json,hashlib,collections
from pathlib import Path
root=Path('D:/Open-Space/campaigns/pr143-des-historical-581c9c4-ea2d60a-20261117/run'); old=Path('results/planet-nine-fir/archives/planet_draws.csv.gz');oldrows=list(csv.DictReader(gzip.open(old,'rt')));newrows=list(csv.DictReader(gzip.open(root/'planet_draws.csv.gz','rt')));out=json.loads((root/'draw_outcomes.json').read_text());flags=['hit_mass5_R2.71_p0.10','hit_mass5_R2.94_p0.30','hit_mass10_R3.70_p0.30'];a={(r['experiment'],r['draw']):r for r in oldrows};b={(r['experiment'],r['draw']):r for r in newrows};print('identities',len(a),len(b),set(a)==set(b));print('status',collections.Counter((r['experiment'],r['state_status'],r['image_status']) for r in out))
for exp in ['published_centroids','centroid_stress']:
 for flag in flags:
  keys=[k for k in a if k[0]==exp];counts=collections.Counter((a[k][flag],b[k][flag]) for k in keys);changes=[k[1] for k in keys if a[k][flag]!=b[k][flag]];print(exp,flag,'old/new',sum(int(a[k][flag]) for k in keys),sum(int(b[k][flag]) for k in keys),'cross',dict(counts),'changes',changes)
print('status', (root/'run_status.json').read_text().strip())
for f in ['run_manifest.json','run_status.json','mapped_states.json','draw_outcomes.json','planet_draws.csv.gz','difference_result.json','model_recovery.json']:
 p=root/f;print(f,len(p.read_bytes()),hashlib.sha256(p.read_bytes()).hexdigest())

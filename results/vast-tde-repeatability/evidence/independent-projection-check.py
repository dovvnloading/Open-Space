import csv,gzip,json,hashlib,math
from pathlib import Path
import pyarrow.parquet as pq
repo=Path('D:/Open-Space/worktrees/vast-holdout-review-20261117');base=repo/'results/vast-tde-expanded';ev=repo/'results/vast-tde-repeatability/evidence'
raw=Path('D:/Open-Space/campaigns/vast-tde-metadata-20261117/VAST_Extragal_DR1_measurement_df.parquet');img=Path('D:/Open-Space/campaigns/vast-tde-expanded-20261117/image_df.csv')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(raw)=='a3af873784b4b3d5a727d6d50890a48a5fd1bbc41ebe035106ce83d47cfb9ca6'
assert sha(img)=='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195'
images={int(r['image_id']):r for r in csv.DictReader(img.open(newline='',encoding='utf-8'))};packet=json.loads((base/'proposed-packet.json').read_text());pools=list(csv.DictReader((base/'control-metadata.csv').open()));manifest=json.loads((ev/'projection-manifest.json').read_text());cols=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings'];checks=[]
def equal(a,b):
    if isinstance(a,float) and isinstance(b,float) and math.isnan(a) and math.isnan(b):return True
    return a==b
for sid,cfg in packet['sources'].items():
    frozen=json.loads((ev/f'{sid}-development.json').read_text())['result']
    for phase in ['development','holdout']:
        ids={int(i) for g in cfg[phase] for i in g}
        ctr={int(r['control_source']) for r in pools if r['target_source']==sid} if phase=='development' else {int(i) for i in frozen['control_baselines']}
        sources=ctr|{int(sid)}
        table=pq.read_table(raw,columns=cols,filters=[('source','in',list(sources)),('image_id','in',list(ids))])
        assert set(table['source'].to_pylist())<=sources and set(table['image_id'].to_pylist())<=ids
        native=table.to_pylist();expected={}
        for r in native:
            r['beam_bmaj']=images[int(r['image_id'])]['beam_bmaj'];r['beam_bmin']=images[int(r['image_id'])]['beam_bmin'];key=(int(r['source']),int(r['image_id']));assert key not in expected;expected[key]=r
        path=ev/f'{sid}-{phase}.jsonl.gz';assert sha(path)==manifest['projections'][path.name]['sha256']
        with gzip.open(path,'rt',encoding='utf-8') as f:deposited=[json.loads(line) for line in f]
        assert len(deposited)==len(native)==manifest['projections'][path.name]['rows']
        for r in deposited:
            key=(int(r['source']),int(r['image_id']));n=expected.pop(key)
            assert n.keys()==r.keys() and all(equal(v,r[k]) for k,v in n.items()),(path.name,key)
        assert not expected
        checks.append({'projection':path.name,'rows':len(native),'raw_values_and_joined_beam_metadata':'exact','sha256':sha(path)})
result={'measurement_sha256':sha(raw),'image_metadata_sha256':sha(img),'checked':checks,'scope':'Post-unblinding independent row provenance check only; unchanged authorized source/image sets and same observed data.'}
out=Path('D:/Open-Space/campaigns/vast-projection-rootcheck-74ad5a5-20261117');out.mkdir(exist_ok=False)
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
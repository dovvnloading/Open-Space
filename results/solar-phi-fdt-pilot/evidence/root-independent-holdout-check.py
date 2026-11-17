"""Root validation: manual bilinear arithmetic on authorized development only."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
CODE = Path('D:/Open-Space/worktrees/solar-phi-fdt-20261117/results/solar-phi-fdt-pilot')
RAW = Path('D:/Open-Space/campaigns/phi-fdt-4568cd0-metadata-terra61')
sys.path.insert(0, str(CODE))
from geometry import read_header
from footprints import footprints

def sample(image, mapping, shift=(0, 0)):
    x = np.asarray(mapping['x']) + shift[0]
    y = np.asarray(mapping['y']) + shift[1]
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    ax, ay = x - ix, y - iy
    corners = np.stack([image[iy, ix], image[iy, ix+1], image[iy+1, ix], image[iy+1, ix+1]])
    assert np.isfinite(corners).all() and np.max(np.abs(corners)) < 4900
    value = ((1-ay)*(1-ax)*corners[0] + (1-ay)*ax*corners[1]
             + ay*(1-ax)*corners[2] + ay*ax*corners[3])
    return float(np.median(value / np.asarray(mapping['mu'])))


dev=Path('D:/Open-Space/campaigns/phi-fdt-development-24ec48d-rootreview-20261117')
prediction_raw=(dev/'prediction.json').read_bytes()
assert hashlib.sha256(prediction_raw).hexdigest()=='973bd0f0a701399fb6a0540a81efd69750d1765895be1ad6e5b16eb61c47fdee'
p=json.loads(prediction_raw)
auth=json.loads((ROOT/'holdout-authorization.json').read_bytes())
assert auth['authorized'] is True and auth['prediction_sha256']==hashlib.sha256(prediction_raw).hexdigest()
original_raw=(ROOT/'primary-result.json').read_bytes()
assert hashlib.sha256(original_raw).hexdigest()=='04824608556edc05bfc206a99f741afd9a1880f7c2d081ebae1d3818e4381917'
original=json.loads(original_raw)
rows=[]
for i,(name,digest,mapping) in enumerate(zip(p['holdout_obs_ids'],p['holdout_sha256'],p['tracked_footprint'][2:])):
    path=RAW/name;assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    image=fits.getdata(path,0).astype(float)
    value=sample(image,mapping)
    excess=p['predicted_sign']*(value-p['baseline_median_gauss']-p['control_median_gauss'])
    assert abs(value-original['scans'][i]['tracked_proxy_gauss'])<1e-10
    rows.append({'name':name,'manual_median_gauss':value,'signed_excess_gauss':excess,'passes':bool(excess>p['threshold_gauss'])})
result={'prediction_sha256':hashlib.sha256(prediction_raw).hexdigest(),'primary_result_sha256':hashlib.sha256(original_raw).hexdigest(),'threshold_gauss':p['threshold_gauss'],'manual_bilinear_results':rows,'both_pass':all(r['passes'] for r in rows),'meaning':'Independently coded arithmetic reproduces failure of the frozen conditional proxy prediction. Same pixels and physical limitations.'}
assert not result['both_pass']
with (ROOT/'root-independent-holdout-review.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))

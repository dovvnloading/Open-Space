"""Independent offline check of frozen packet, g42 stencils and decisions."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
EVIDENCE=ROOT/'evidence'
def read(name):return json.loads((EVIDENCE/name).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def close(a,b):
    if not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-9):raise AssertionError(f'{a} != {b}')

def verify():
    packet=ROOT/'PACKET.json'
    if sha(packet)!='4827108044a5c777728cd85c455be20eaf7610f866059f384279c40b9404138b':raise AssertionError('packet hash')
    cfg=json.loads(packet.read_text(encoding='utf-8'))
    for name,digest in cfg['file_sha256'].items():
        if sha(ROOT/name)!=digest:raise AssertionError('frozen file '+name)
    for name,expected in [('development.json','f64370e6128241a487baae281c193db26a37fe449692d7f9fab7fd0aecbeb93f'),('prediction.json','973bd0f0a701399fb6a0540a81efd69750d1765895be1ad6e5b16eb61c47fdee'),('primary-result.json','04824608556edc05bfc206a99f741afd9a1880f7c2d081ebae1d3818e4381917'),('secondary.json','6aafc629b0c61c130ebd6e9b4fd317d8084b012bf26cc13b151bafece71afb59')]:
        if sha(EVIDENCE/name)!=expected:raise AssertionError('evidence hash '+name)
    d=read('development.json');p=read('prediction.json');primary=read('primary-result.json');secondary=read('secondary.json')
    rows=d['rows'];valid=[r for r in rows if r['status']=='usable'];assert len(rows)==25 and len(valid)==24 and [r['patch_id'] for r in rows if r['status']!='usable']==['g13']
    top=max(valid,key=lambda x:abs(x['delta']));assert top['patch_id']==p['patch_id']=='g42'
    controls=np.array([r['delta'] for r in valid if r is not top]);median=float(np.median(controls));scale=max(1.,float(1.4826*np.median(np.abs(controls-median))))
    close(median,p['control_median_gauss']);close(scale,p['control_scale_gauss']);close(3*scale,p['threshold_gauss']);assert abs(top['delta']-median)>=5*scale
    for v in top['registration_shifts']:assert np.sign(v-median)==p['predicted_sign']
    with np.load(EVIDENCE/'g42-stencils.npz',allow_pickle=False) as z:
        s=z['stencils'];x=z['fraction_x'];y=z['fraction_y'];mu=z['mu'];assert s.shape==(4,317,4) and x.shape==y.shape==mu.shape==(4,317)
        interp=(1-y)*((1-x)*s[:,:,0]+x*s[:,:,1])+y*((1-x)*s[:,:,2]+x*s[:,:,3])
        med=np.median(interp/mu,axis=1)
    expected=[top['baseline'],top['development']]+[x['tracked_proxy_gauss'] for x in primary['scans']]
    for x,y in zip(med,expected):close(x,y)
    for i,entry in enumerate(primary['scans']):
        score=p['predicted_sign']*(med[i+2]-p['baseline_median_gauss']-median)
        close(score,entry['signed_excess_gauss']);assert entry['passes_prediction']==(score>p['threshold_gauss'])
    assert primary['both_pass']==all(x['passes_prediction'] for x in primary['scans'])==False
    for i in range(2):
        for base in ('early','late'):
            values=np.array([r['scans'][i]['change_from_'+base] for r in secondary['rows'] if r['patch_id']!='g42' and r['scans'][i]['status']=='usable'])
            rec=next(c for c in secondary['controls'] if c['obs_id']==primary['scans'][i]['obs_id'] and c['reference']==base)
            assert len(values)==rec['count'];m=np.median(values);close(m,rec['median_change_gauss']);close(1.4826*np.median(abs(values-m)),rec['robust_mad_gauss'])
    print('24/25 eligible; g42 medians:',', '.join(f'{v:.6f}' for v in med))
    print('Held-out scores:',', '.join(f"{x['signed_excess_gauss']:.6f}" for x in primary['scans']),'threshold',p['threshold_gauss'],'both_pass',primary['both_pass'])
    print('Packet, observational stencil interpolation and reported controls passed.')

if __name__=='__main__':verify()

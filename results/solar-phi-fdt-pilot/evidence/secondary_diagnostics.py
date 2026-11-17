"""Predeclared post-primary diagnostics, separate from frozen test."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
sys.path.insert(0,'D:/Open-Space/worktrees/solar-phi-fdt-20261117/results/solar-phi-fdt-pilot')
from geometry import read_header
from footprints import footprints
from development import sample,verify_packet

ROOT=Path('D:/Open-Space/campaigns/phi-fdt-4568cd0-metadata-terra61')
DEV=Path('D:/Open-Space/campaigns/phi-fdt-development-24ec48d-rootreview-20261117')
OUT=Path('D:/Open-Space/campaigns/phi-fdt-holdout-24ec48d-973bd0f-20261117')
NAMES=['solo_L2_phi-fdt-blos_20260813T102509_V03_0648130504.fits','solo_L2_phi-fdt-blos_20260813T162509_V03_0648130506.fits','solo_L2_phi-fdt-blos_20260815T162509_V03_0648150506.fits','solo_L2_phi-fdt-blos_20260815T192509_V03_0648150507.fits']
HASHES=['9bada5a8e35cbc2fc995afd23c169d9633efb8c50befb10393fa99e58e53eb63','691123a01f7ca4e3fcff911d09e74bf98294774d2e8d154dd5a2f254e29905a7','783784008c4f4f4c0d017e89f94d22d327a6f4a987b1af0cbbf303087d358664','e9a974446e3e77d44bf08d1652d58050ccfd8974923926cd76410ae80e06e1cd']

def filehash(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    resultfile=OUT/'secondary.json'
    if resultfile.exists():raise FileExistsError('secondary result must be new')
    if verify_packet()!='4827108044a5c777728cd85c455be20eaf7610f866059f384279c40b9404138b':raise ValueError('packet mismatch')
    if filehash(OUT/'primary-result.json')!='04824608556edc05bfc206a99f741afd9a1880f7c2d081ebae1d3818e4381917':raise ValueError('primary result absent or changed')
    if filehash(DEV/'development.json')!='f64370e6128241a487baae281c193db26a37fe449692d7f9fab7fd0aecbeb93f':raise ValueError('development mismatch')
    for n,sha in zip(NAMES,HASHES):
        if filehash(ROOT/n)!=sha:raise ValueError('input hash mismatch '+n)
    headers=[read_header(ROOT/n) for n in NAMES];patches=footprints(headers[1],headers)
    old=json.loads((DEV/'development.json').read_text());orig={r['patch_id']:r for r in old['rows']}
    images=[]
    for n in NAMES[2:]:
        with fits.open(ROOT/n,memmap=True) as hdus:images.append(np.array(hdus[0].data,dtype=float,copy=True))
    rows=[]
    for p in patches:
        prior=orig[p['patch_id']];row={'patch_id':p['patch_id'],'geometry_eligible':p['eligible'],'geometry_reasons':p['reasons'],'development_status':prior['status'],'early':prior.get('baseline'),'late':prior.get('development'),'scans':[]}
        for image,mapping,name in zip(images,p['maps'][2:],NAMES[2:]):
            entry={'obs_id':name,'status':'geometry_excluded'}
            if p['eligible'] and prior['status']=='usable':
                try:
                    value=sample(image,mapping);entry.update(status='usable',proxy=value,change_from_early=value-prior['baseline'],change_from_late=value-prior['development'],shifts={})
                    for label,offset in [('x_minus',(-1,0)),('x_plus',(1,0)),('y_minus',(0,-1)),('y_plus',(0,1))]:
                        try:
                            shifted=sample(image,mapping,offset);entry['shifts'][label]={'status':'usable','proxy':shifted,'change_from_early':shifted-prior['baseline'],'change_from_late':shifted-prior['development']}
                        except ValueError as exc:entry['shifts'][label]={'status':'invalid','reason':str(exc)}
                except ValueError as exc:entry.update(status='invalid_stencil',reason=str(exc))
            row['scans'].append(entry)
        rows.append(row)
    controls=[]
    for scan in range(2):
        for baseline in ('early','late'):
            vals=np.array([r['scans'][scan]['change_from_'+baseline] for r in rows if r['patch_id']!='g42' and r['scans'][scan]['status']=='usable'])
            med=float(np.median(vals)) if len(vals) else None;mad=float(1.4826*np.median(np.abs(vals-med))) if len(vals) else None
            controls.append({'obs_id':NAMES[scan+2],'reference':baseline,'count':len(vals),'median_change_gauss':med,'robust_mad_gauss':mad})
    output={'primary_sha256':filehash(OUT/'primary-result.json'),'packet_sha256':verify_packet(),'development_sha256':filehash(DEV/'development.json'),'diagnostic_source_sha256':filehash(Path(__file__)),'rows':rows,'controls':controls,'interpretation':'Secondary diagnostics only; no candidate selection, threshold change, inference rescue or independent false-alarm calibration.'}
    resultfile.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8');print('wrote',resultfile)

if __name__=='__main__':main()

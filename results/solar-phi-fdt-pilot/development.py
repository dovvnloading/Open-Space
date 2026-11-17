"""Development-only candidate selection on two approved August 13 maps.

Do not run on observed images before the exact packet is reviewed. Holdout IDs
are rejected before opening any primary image array. No discovery claim.
"""
import argparse
import hashlib
import json
import subprocess
import re
import sys
from pathlib import Path
import numpy as np
import scipy
import astropy
import sunpy
from scipy.ndimage import map_coordinates
from astropy.io import fits
from geometry import read_header
from footprints import footprints,detector_fixed
from metadata_audit import audit

ALLOWED=('solo_L2_phi-fdt-blos_20260813T102509_V03_0648130504.fits','solo_L2_phi-fdt-blos_20260813T162509_V03_0648130506.fits')
HASHES=('9bada5a8e35cbc2fc995afd23c169d9633efb8c50befb10393fa99e58e53eb63','691123a01f7ca4e3fcff911d09e74bf98294774d2e8d154dd5a2f254e29905a7')

def sample(image,mapping,shift=(0.,0.)):
    x=np.asarray(mapping['x'])+shift[0];y=np.asarray(mapping['y'])+shift[1];mu=np.asarray(mapping['mu'])
    if image.ndim!=2 or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('invalid sample geometry')
    ix=np.floor(x).astype(int);iy=np.floor(y).astype(int)
    if np.any(ix<0)|np.any(iy<0)|np.any(ix+1>=image.shape[1])|np.any(iy+1>=image.shape[0]):raise ValueError('stencil outside image')
    stencil=np.stack((image[iy,ix],image[iy+1,ix],image[iy,ix+1],image[iy+1,ix+1]))
    if not np.isfinite(stencil).all() or np.any(np.abs(stencil)>=4900):raise ValueError('invalid or near-saturated stencil')
    if not np.isfinite(mu).all() or np.any(mu<.65):raise ValueError('invalid projection')
    val=map_coordinates(image,[y,x],order=1,mode='constant',cval=np.nan,prefilter=False)
    return float(np.median(val/mu))

def assess(first,second,patches):
    rows=[]
    for p in patches:
        row={'patch_id':p['patch_id'],'eligible_geometry':p['eligible'],'status':'geometry_excluded','reasons':p['reasons']}
        if p['eligible']:
            try:
                a=sample(first,p['maps'][0]);b=sample(second,p['maps'][1]); delta=b-a
                shifts=[sample(second,p['maps'][1],shift)-a for shift in ((-1,0),(1,0),(0,-1),(0,1))]
                row.update(status='usable',baseline=a,development=b,delta=delta,registration_shifts=shifts)
            except ValueError as exc:row.update(status='invalid_stencil',detail=str(exc))
        rows.append(row)
    usable=[r for r in rows if r['status']=='usable']
    if len(usable)<8:return {'status':'inconclusive_fewer_than_eight_patches','rows':rows,'candidate':None}
    candidate=max(usable,key=lambda r:abs(r['delta'])) # original grid order breaks ties
    controls=np.array([r['delta'] for r in usable if r is not candidate]);median=float(np.median(controls));scale=max(1.,float(1.4826*np.median(np.abs(controls-median))))
    signal=candidate['delta']-median
    robust=all(np.sign(v-median)==np.sign(signal) for v in candidate['registration_shifts'])
    admitted=abs(signal)>=5*scale and robust
    return {'status':'candidate_proposed' if admitted else 'no_candidate','candidate':candidate['patch_id'] if admitted else None,'raw_top_patch':candidate['patch_id'],'control_median':median,'control_scale':scale,'threshold_scale':5.,'registration_sign_stable':robust,'grid_count':len(rows),'geometry_eligible':sum(r['eligible_geometry'] for r in rows),'usable':len(usable),'invalid_stencils':sum(r['status']=='invalid_stencil' for r in rows),'rows':rows}

def verify_packet():
    folder=Path(__file__).resolve().parent
    packet=json.loads((folder/'PACKET.json').read_text(encoding='utf-8'))
    versions={'numpy':np.__version__,'scipy':scipy.__version__,'astropy':astropy.__version__,'sunpy':sunpy.__version__}
    if versions!=packet['package_versions'] or list(sys.version_info[:2])!=packet['python_major_minor']:raise ValueError('software version mismatch')
    for name,digest in packet['file_sha256'].items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=digest:raise ValueError('packet byte mismatch: '+name)
    return hashlib.sha256((folder/'PACKET.json').read_bytes()).hexdigest()

def read_dev(paths):
    if len(paths)!=2:raise ValueError('exactly two development maps required')
    images=[];headers=[]
    for path,expected_name,expected_sha in zip(paths,ALLOWED,HASHES):
        path=Path(path)
        if path.name!=expected_name:raise ValueError('only frozen development filenames are permitted')
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=expected_sha:raise ValueError('development hash mismatch')
        header=read_header(path)
        if header['FILENAME']!=expected_name:raise ValueError('FITS filename mismatch')
        headers.append(header)
    # Only after every development file passed identity and SHA checks.
    for path in paths:
        with fits.open(path,memmap=True) as h:images.append(np.array(h[0].data,dtype=float,copy=True))
    return images,headers

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('early');parser.add_argument('late');parser.add_argument('--holdout-meta-a',required=True);parser.add_argument('--holdout-meta-b',required=True);parser.add_argument('--run-id',required=True);parser.add_argument('--code-sha',required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    if args.out.exists():raise FileExistsError('output must be new')
    packet_sha=verify_packet()
    if not re.fullmatch(r'[a-zA-Z0-9_-]+',args.run_id):raise ValueError('invalid run ID')
    actual=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    if args.code_sha!=actual:raise ValueError('code SHA does not match checkout')
    subprocess.run(['git','diff','--exit-code','HEAD','--','results/solar-phi-fdt-pilot'],check=True,stdout=subprocess.DEVNULL)
    hold_names=('solo_L2_phi-fdt-blos_20260815T162509_V03_0648150506.fits','solo_L2_phi-fdt-blos_20260815T192509_V03_0648150507.fits')
    hold_hashes=('783784008c4f4f4c0d017e89f94d22d327a6f4a987b1af0cbbf303087d358664','e9a974446e3e77d44bf08d1652d58050ccfd8974923926cd76410ae80e06e1cd')
    hold_headers=[]
    for filename,name,digest in zip((args.holdout_meta_a,args.holdout_meta_b),hold_names,hold_hashes):
        path=Path(filename)
        if path.name!=name or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise ValueError('holdout metadata file mismatch')
        h=read_header(path)
        if h['FILENAME']!=name:raise ValueError('holdout FITS identity mismatch')
        hold_headers.append(h) # header only: holdout HDU 0 .data is never accessed
    # Compare acquisition metadata for all four before touching either science image.
    provenance=[audit(Path(name)) for name in (args.early,args.late,args.holdout_meta_a,args.holdout_meta_b)]
    for i in range(len(provenance)):
        for j in range(i):
            if set(provenance[i]['timestamps'])&set(provenance[j]['timestamps']):raise ValueError('shared acquisition times')
            if set(provenance[i]['operand_ids'])&set(provenance[j]['operand_ids']):raise ValueError('shared recorded processing operand')
            if provenance[j]['date_end'] >= provenance[i]['date_beg']:raise ValueError('overlapping or unordered scan windows')
    for item in provenance:item.pop('timestamps')
    images,headers=read_dev((args.early,args.late))
    patches=footprints(headers[1],headers+hold_headers)
    result=assess(*images,patches)
    if result['candidate'] is not None:
        p=next(z for z in patches if z['patch_id']==result['candidate'])
        result['tracked_footprint']=[{'filename':header['FILENAME'],'x':mapping['x'].tolist(),'y':mapping['y'].tolist(),'mu':mapping['mu'].tolist()} for header,mapping in zip(headers+hold_headers,p['maps'])]
        reference=p['maps'][1];controls=[]
        for header in headers+hold_headers:
            z=detector_fixed(headers[1],header,reference)
            controls.append({'filename':header['FILENAME'],'eligible':z['eligible'],'invalid_nodes':z['invalid_nodes'],'x':z['x'].tolist() if z['eligible'] else None,'y':z['y'].tolist() if z['eligible'] else None,'mu':z['mu'].tolist() if z['eligible'] else None})
        try:
            m0=detector_fixed(headers[1],headers[0],reference)
            detector_baseline=sample(images[0],m0) if m0['eligible'] else None
            detector_development=sample(images[1],reference)
        except ValueError:
            detector_baseline=detector_development=None
        result['detector_frame_control']={'reference_pxbeg1':headers[1]['PXBEG1'],'reference_pxbeg2':headers[1]['PXBEG2'],'baseline_median':detector_baseline,'development_median':detector_development,'maps':controls,'interpretation':'fixed detector location traces a different evolving solar surface; diagnostic only'}
    result.update({'run_id':args.run_id,'code_sha':actual,'packet_sha256':packet_sha,'acquisition_metadata':provenance,'development_obs_ids':list(ALLOWED),'development_sha256':list(HASHES),'holdout_obs_ids':list(hold_names),'holdout_sha256':list(hold_hashes),'holdout_pixels_read':False,'status_note':'DEVELOPMENT ONLY: freeze a full numeric prediction and review before any holdout pixel access.'})
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

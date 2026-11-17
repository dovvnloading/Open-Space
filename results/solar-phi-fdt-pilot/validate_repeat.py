"""Future held-out evaluator; not authorized until a prediction is frozen."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from astropy.io import fits
from geometry import read_header
from development import sample,verify_packet


def evaluate(p,images):
    if p['status']!='prediction_frozen' or len(images)!=2:raise ValueError('no usable frozen prediction')
    tracked=p['tracked_footprint'][2:];detector=p['detector_frame_control'];control_maps=detector['maps'][2:];out=[]
    for image,expected,mapping,ctrl in zip(images,p['holdout_obs_ids'],tracked,control_maps):
        if mapping['filename']!=expected or ctrl['filename']!=expected:raise ValueError('held-out footprint identity mismatch')
        value=sample(image,mapping);signal=p['predicted_sign']*(value-p['baseline_median_gauss']-p['control_median_gauss'])
        det=None
        if ctrl['eligible'] and detector['baseline_median'] is not None:
            try:det=sample(image,ctrl)-detector['baseline_median']
            except ValueError:pass
        out.append({'obs_id':expected,'tracked_proxy_gauss':value,'signed_excess_gauss':signal,'passes_prediction':bool(signal>p['threshold_gauss']),'fixed_detector_change_gauss':det})
    return {'status':'evaluated','scans':out,'both_pass':all(x['passes_prediction'] for x in out),'meaning':'conditional signed radial-proxy persistence screen; detector control is diagnostic; no p-value or true flux-emergence claim'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prediction',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True);ap.add_argument('--holdout-a',type=Path,required=True);ap.add_argument('--holdout-b',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError('output must be new')
    packet=verify_packet();raw=a.prediction.read_bytes();p=json.loads(raw);auth=json.loads(a.authorization.read_text(encoding='utf-8'))
    head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    if auth.get('authorized') is not True or auth.get('prediction_sha256')!=hashlib.sha256(raw).hexdigest() or auth.get('code_sha')!=head or p.get('development_code_sha')!=head or p.get('development_packet_sha256')!=packet:raise ValueError('authorization, code or frozen prediction mismatch')
    subprocess.run(['git','diff','--exit-code','HEAD','--','results/solar-phi-fdt-pilot'],check=True,stdout=subprocess.DEVNULL)
    paths=(a.holdout_a,a.holdout_b)
    # Validate both hashes and headers before opening either held-out science array.
    for path,name,digest in zip(paths,p['holdout_obs_ids'],p['holdout_sha256']):
        if path.name!=name or hashlib.sha256(path.read_bytes()).hexdigest()!=digest or read_header(path)['FILENAME']!=name:raise ValueError('holdout identity/hash mismatch')
    images=[]
    for path in paths:
        with fits.open(path,memmap=True) as h:images.append(np.array(h[0].data,dtype=float,copy=True))
    result=evaluate(p,images);result.update({'prediction_sha256':hashlib.sha256(raw).hexdigest(),'authorization':auth,'code_sha':head,'packet_sha256':packet})
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

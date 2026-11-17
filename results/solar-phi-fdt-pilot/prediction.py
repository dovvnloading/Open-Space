"""Freeze a development-selected signed persistence prediction; no FITS access."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def construct(d,development_sha):
    if d.get('holdout_pixels_read') is not False:raise ValueError('holdout isolation not established')
    if d['status']!='candidate_proposed':return {'status':'no_prediction','development_result_sha256':development_sha,'reason':d['status']}
    rows=[r for r in d['rows'] if r['patch_id']==d['candidate']]
    if len(rows)!=1 or rows[0]['status']!='usable':raise ValueError('candidate inconsistent')
    row=rows[0];signal=row['delta']-d['control_median'];scale=d['control_scale']
    if not np.isfinite([signal,scale,row['baseline']]).all() or abs(signal)<5*scale or scale<1:raise ValueError('invalid development candidate')
    return {'status':'prediction_frozen','development_result_sha256':development_sha,'development_run_id':d['run_id'],'development_code_sha':d['code_sha'],'development_packet_sha256':d['packet_sha256'],'patch_id':d['candidate'],'observable':'median(B_LOS / mu), radial-field proxy in Gauss','baseline_obs_id':d['development_obs_ids'][0],'baseline_median_gauss':row['baseline'],'control_median_gauss':d['control_median'],'control_scale_gauss':scale,'predicted_sign':int(np.sign(signal)),'threshold_gauss':3*scale,'holdout_obs_ids':d['holdout_obs_ids'],'holdout_sha256':d['holdout_sha256'],'tracked_footprint':d['tracked_footprint'],'detector_frame_control':d['detector_frame_control'],'rule':'for both scans: predicted_sign * (holdout_median - baseline_median - control_median) > threshold_gauss','interpretation':'conditional persistence screen; not a calibrated p-value or magnetic-flux-emergence claim'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('development',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('output must be new')
    raw=a.development.read_bytes();result=construct(json.loads(raw),hashlib.sha256(raw).hexdigest())
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

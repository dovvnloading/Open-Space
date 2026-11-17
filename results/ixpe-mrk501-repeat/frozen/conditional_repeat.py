"""Frozen conditional repeatability screen; no download functionality.

Do not run with holdout data until the reviewed packet is authorized.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from centroid_rule import baseline as centroid_for_data
from development import ROOT, CAL, CENTER, primary, process

HERE=Path(__file__).resolve().parent
OBS='05004901'
THRESHOLD=5.7e-7

def load(path):return json.loads(Path(path).read_text(encoding='utf-8'))

def verify_receipt(root,receipt,manifest):
    records={x['relative']:x for x in receipt}
    needed={x['relative']:x for x in manifest['files']}
    if len(records)!=len(receipt) or len(needed)!=len(manifest['files']) or set(records)!=set(needed):raise ValueError('receipt file set mismatch')
    for name,spec in needed.items():
        rec=records[name]
        if rec['url']!=spec['url'] or rec['bytes']!=spec['bytes']:raise ValueError('receipt URL or expected size mismatch')
        p=root/name
        if not p.is_file() or p.stat().st_size!=rec['bytes']:raise ValueError('missing/wrong size '+name)
        sha=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):sha.update(block)
        if sha.hexdigest()!=rec['sha256']:raise ValueError('hash mismatch '+name)

def validate_gti(a,b):
    if not len(a) or not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(b<=a) or np.any(np.diff(a)<0) or np.any(a[1:]<b[:-1]):raise ValueError('invalid GTI')

def gti_list(root,obsid):
    out=[]
    for du in (1,3):
        with fits.open(root/f'event_l2/ixpe{obsid}_det{du}_evt2_v01.fits.gz') as h:
            hdr=h['EVENTS'].header;d=h['EVENTS'].data;g=h['GTI'].data
            if hdr['TIMESYS']!='TT' or hdr['MJDREFI']!=57754 or abs(hdr['MJDREFF']-0.00080074074074)>1e-13:raise ValueError('time reference mismatch')
            if hdr.get('OBS_ID')!=obsid:raise ValueError('observation ID mismatch')
            t=np.asarray(d['TIME']);a=np.asarray(g['START']);b=np.asarray(g['STOP'])
            validate_gti(a,b)
            if not np.isfinite(t).all():raise ValueError('invalid event times')
            idx=np.searchsorted(a,t,side='right')-1
            if np.any((idx<0)|(t>b[np.maximum(idx,0)])):raise ValueError('event outside GTI')
            if not all(np.isfinite(d[col]).all() for col in ('Q','U','X','Y')):raise ValueError('non-finite event columns')
            flags=np.asarray(d['STATUS2'],dtype=bool);unexpected=flags.copy();unexpected[:,[1,12,15]]=False
            if unexpected.any():raise ValueError('unexpected STATUS2 flag')
            out.extend(zip(a.tolist(),b.tolist()))
    return out

def ensure_disjoint(a,b):
    for x,y in a:
        for u,v in b:
            if max(x,u)<min(y,v):raise ValueError('baseline and holdout GTIs overlap')

def compare(baseline,holdout):
    out={}
    for name in ('nominal','xy_centroid'):
        a=baseline[name];b=holdout[name]
        d=np.array([b['q']-a['q'],b['u']-a['u']]);c1=np.asarray(a['cov_qu'],dtype=float);c2=np.asarray(b['cov_qu'],dtype=float)
        for cov in (c1,c2):
            if cov.shape!=(2,2) or not np.isfinite(cov).all() or not np.allclose(cov,cov.T) or np.linalg.eigvalsh(cov).min()<=0:raise ValueError('invalid component covariance')
        c=c1+c2
        if not np.isfinite(d).all() or not np.isfinite(c).all() or np.linalg.eigvalsh(c).min()<=0:raise ValueError('invalid comparison covariance')
        chi=float(d@np.linalg.solve(c,d));p=float(np.exp(-chi/2))
        out[name]={'difference':d.tolist(),'covariance':c.tolist(),'chi2':chi,'conditional_p':p,'screening_crossing':p<THRESHOLD}
    return {'threshold':THRESHOLD,'results':out,'interpretation':'Conditional Gaussian screening only; systematics and independent confirmation remain unresolved.'}

def verify_packet():
    packet=load(HERE/'conditional_packet.json')
    if packet['threshold_conditional_p']!=THRESHOLD or version('ixpeobssim')!='33.2.1':raise ValueError('packet threshold or software mismatch')
    for name,expected in packet['repository_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=expected:raise ValueError('packet file mismatch '+name)
    for rel,expected in packet['caldb_sha256'].items():
        if hashlib.sha256((CAL/rel).read_bytes()).hexdigest()!=expected:raise ValueError('calibration mismatch '+rel)
    for rel,expected in packet['baseline_event_sha256'].items():
        if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=expected:raise ValueError('baseline event mismatch '+rel)

def run(root,receipt):
    verify_packet()
    manifest=load(HERE/'holdout_inputs.json');verify_receipt(root,receipt,manifest)
    # These calls inspect holdout event science arrays; execute only after authorization.
    holdout_gtis=gti_list(root,OBS);baseline_gtis=gti_list(ROOT,'05004899');ensure_disjoint(baseline_gtis,holdout_gtis)
    pos=centroid_for_data(root=root,obsid=OBS)
    center=SkyCoord(float(pos['ra_deg'])*u.deg,float(pos['dec_deg'])*u.deg)
    h={'nominal':primary(CENTER,root,OBS,(OBS,)),'xy_centroid':primary(center,root,OBS,(OBS,))}
    for item in h.values():
        for du in item['du_sufficient']:
            if du['source_count']<1000 or du['background_count']<100:raise ValueError('insufficient source or background counts')
    b=load(HERE/'frozen_baseline.json')
    try:diagnostic=process(2,CENTER,root,OBS,(OBS,))
    except ValueError as exc:diagnostic={'failed':str(exc)}
    return {'obsid':OBS,'receipt':receipt,'packet_sha256':hashlib.sha256((HERE/'conditional_packet.json').read_bytes()).hexdigest(),'du2_diagnostic':diagnostic,'centroid':{'ra_deg':float(pos['ra_deg']),'dec_deg':float(pos['dec_deg'])},'summary':h,'comparison':compare(b,h)}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);parser.add_argument('--receipt',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
    output=Path(args.output)
    if output.exists():raise FileExistsError('output must be new')
    result=run(Path(args.root),load(args.receipt))
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

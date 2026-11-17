"""Conditional VAST screen; functions only. No implicit archive or network access.
Inputs must already be phase-filtered by IDs; parquet decoding may read whole
compressed column chunks. This is analysis isolation, not physical byte isolation.
"""
import math
import numpy as np

MIN_CONTROLS=20
MIN_DEV_SESSIONS=4
CONTROL_MIN_SNR=5.0
CONTROL_MAX_FRACTIONAL_MAD=0.10
GAIN_FLOOR=0.03
Z_THRESHOLD=3.0

def valid(r):
    required=('source','image_id','flux_peak','flux_peak_err','local_rms','beam_bmaj','beam_bmin','forced','flag_c4','has_siblings')
    if any(k not in r for k in required):return False,'missing field'
    try:
        x,s,rms,b1,b2=(float(r[k]) for k in ('flux_peak','flux_peak_err','local_rms','beam_bmaj','beam_bmin'))
        if not np.isfinite([x,s,rms,b1,b2]).all() or min(s,rms,b1,b2)<=0:return False,'invalid flux/error/rms/beam'
        if rms>1.82:return False,'local RMS exceeds 1.82 mJy/beam'
        if r['forced'] not in (True,False):return False,'missing forced flag'
        if r['flag_c4'] is True or r['has_siblings'] is True:return False,'fit flag'
        if not r['forced'] and (r['flag_c4'] is not False or r['has_siblings'] is not False):return False,'missing detection flags'
    except (TypeError,ValueError):return False,'malformed value'
    return True,''

def phase_rows(rows,allowed_images,allowed_sources):
    seen=set();kept=[];drops=[]
    for r in rows:
        key=(int(r['source']),int(r['image_id']))
        if key in seen:raise ValueError('duplicate source/image')
        seen.add(key)
        if key[0] not in allowed_sources or key[1] not in allowed_images:raise ValueError('phase or source leak')
        ok,why=valid(r)
        (kept if ok else drops).append(r if ok else {'source':key[0],'image_id':key[1],'reason':why})
    return kept,drops

def stable_controls(rows,control_ids,dev_sessions):
    im_to_session={i:j for j,g in enumerate(dev_sessions) for i in g}; chosen={}
    for sid in control_ids:
        per={}
        for r in rows:
            if int(r['source'])==sid and not r['forced'] and float(r['flux_peak'])/float(r['flux_peak_err'])>=CONTROL_MIN_SNR:
                per.setdefault(im_to_session[int(r['image_id'])],[]).append(float(r['flux_peak']))
        vals=np.array([np.median(v) for v in per.values()])
        if len(vals)<MIN_DEV_SESSIONS:continue
        med=float(np.median(vals));mad=1.4826*float(np.median(np.abs(vals-med)))
        if med>0 and mad/med<=CONTROL_MAX_FRACTIONAL_MAD:chosen[sid]=med
    if len(chosen)<MIN_CONTROLS:raise ValueError('too few stable development controls')
    return chosen

def image_gains(rows,controls,image_ids):
    gains={}
    for image in image_ids:
        ratios=[float(r['flux_peak'])/controls[int(r['source'])] for r in rows if int(r['image_id'])==image and int(r['source']) in controls and not r['forced']]
        if len(ratios)<MIN_CONTROLS:raise ValueError('insufficient frozen controls for image')
        g=float(np.median(ratios));mad=1.4826*float(np.median(np.abs(np.array(ratios)-g)))
        if not np.isfinite(g) or g<=0:raise ValueError('invalid gain')
        gains[image]={'gain':g,'fractional_sigma':max(GAIN_FLOOR,1.2533*mad/(g*np.sqrt(len(ratios)))),'n':len(ratios)}
    return gains

def session_measurements(rows,source,sessions,gains):
    out=[]
    for group in sessions:
        r=[x for x in rows if int(x['source'])==source and int(x['image_id']) in group]
        if not r:raise ValueError('missing source session')
        vals=np.array([float(x['flux_peak'])/gains[int(x['image_id'])]['gain'] for x in r]);errs=np.array([float(x['flux_peak_err'])/gains[int(x['image_id'])]['gain'] for x in r]);value=float(np.mean(vals))
        # No 1/sqrt(n) gain for overlapping images: maximum quoted error and
        # the largest fractional gain uncertainty remain fully correlated.
        variance=float(max(errs)**2+(max(gains[int(x['image_id'])]['fractional_sigma'] for x in r)*abs(value))**2)
        out.append((value,variance))
    return out

def baseline(session_values):
    if len(session_values)<MIN_DEV_SESSIONS:raise ValueError('insufficient development sessions')
    x=np.array([v for v,_ in session_values]);v=np.array([q for _,q in session_values]);w=1/v;mean=float(np.sum(x*w)/np.sum(w));chi=float(np.sum((x-mean)**2/v));var=float(max(1.,chi/(len(x)-1))/np.sum(w));return {'mean':mean,'variance':var,'chi2':chi,'dof':len(x)-1}

def compare(pred,holdout):
    if len(holdout)!=2:raise ValueError('two holdout sessions required')
    z=[(x-pred['mean'])/math.sqrt(v+pred['variance']) for x,v in holdout]
    return {'z':z,'both_at_least_3':all(q>=Z_THRESHOLD for q in z),'conditional_only':True}

def develop(source,packet,rows,control_ids):
    config=packet['sources'][str(source)]
    images={i for group in config['development'] for i in group}
    allowed=set(control_ids)|{source}
    good,drops=phase_rows(rows,images,allowed)
    selected=stable_controls(good,set(control_ids),config['development'])
    gains=image_gains(good,selected,images)
    session=session_measurements(good,source,config['development'],gains)
    pred=baseline(session)
    return {'source':source,'control_baselines':selected,'development_gains':gains,'development_sessions':session,'prediction':pred,'drops':drops,'holdout_images':config['holdout']}

def evaluate_holdout(frozen,rows):
    source=int(frozen['source']);groups=frozen['holdout_images'];images={i for g in groups for i in g}
    controls={int(k):float(v) for k,v in frozen['control_baselines'].items()}
    good,drops=phase_rows(rows,images,set(controls)|{source})
    gains=image_gains(good,controls,images)
    sessions=session_measurements(good,source,groups,gains)
    return {'result':compare(frozen['prediction'],sessions),'holdout_gains':gains,'holdout_sessions':sessions,'drops':drops}

def freeze_development(path,result,packet_sha256,code_sha256):
    """Write an immutable record for later review, after authorized development."""
    import json,hashlib
    from pathlib import Path
    p=Path(path);data={'packet_sha256':packet_sha256,'code_sha256':code_sha256,'result':result}
    raw=(json.dumps(data,sort_keys=True,indent=2)+'\n').encode()
    with p.open('xb') as f:f.write(raw)
    return hashlib.sha256(raw).hexdigest()

def load_reviewed_development(path,expected_sha256,packet_sha256,code_sha256):
    import json,hashlib
    from pathlib import Path
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_sha256:raise ValueError('frozen development digest mismatch')
    data=json.loads(raw)
    if data['packet_sha256']!=packet_sha256 or data['code_sha256']!=code_sha256:raise ValueError('development lineage mismatch')
    return data['result']

"""Fixed-window injection experiment for the extended, locked-v1 DESI search."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import io
import json
from pathlib import Path

import numpy as np

from injection_window import REST_MIN, REST_MAX, read_window
from laser_search import inject_line, scan_target
from recovery import _matches
from summarize_survey import compact_rows
from survey_injections import local_response_fwhm, noise_scale


def star(job):
    path,headers,n = job
    meta, records = read_window(path)
    tid = int(meta['targetid'])
    expids = [int(x['expid']) for x in meta['exposures']]
    base = scan_target(records,expids)
    seed = tid*4  # same b-arm convention, different predeclared wavelength domain
    wave = np.random.default_rng(seed).uniform(REST_MIN,REST_MAX,n)
    usable = any(max(a,REST_MIN) < min(b,REST_MAX) for a,b in base.coverage)
    scales = [noise_scale(e) for e in records] if usable else []
    if not usable or None in scales:
        return {'targetid':tid,'arm':'b','seed':seed,'expids':'|'.join(map(str,expids)),
                'trials':n,'n_joint_coverage':0,'n_joint_reference_eligible':0,
                'states_by_draw':'GGG'*n,'status':'no_joint_coverage'}
    fluxes = [m*max(scales) for m in (8.,24.,72.)]
    positions=[]
    for w in wave:
        ok=any(a<=w<=b for a,b in base.coverage)
        confused=any(_matches(c,w) for c in base.candidates)
        fwhm=[local_response_fwhm(e,w*e.doppler) for e in records] if ok else []
        reasons, fitted = [], []
        for flux in fluxes:
            fit=None
            if not ok:
                reason='outside_joint_coverage'
            elif confused:
                reason='baseline_coincidence'
            else:
                injected=[inject_line(e,w,flux) for e in records]
                result=scan_target(injected,expids)
                near=[c for c in result.candidates if _matches(c,w)]
                if near:
                    c=min(near,key=lambda c:abs(c.rest_A-w))
                    hits={h.expid:h for h in c.hits}
                    fit=[[float(hits[e.expid].observed_A),float(hits[e.expid].rest_fwhm_A*e.doppler)]
                         for e in records]
                    reason='recovered'
                else:
                    reason='missed_detection'
            reasons.append(reason)
            fitted.append(fit)
        positions.append(dict(rest_A=float(w),eligible=ok,fwhm_A=fwhm,reasons=reasons,fitted=fitted))
    epochs=[dict(expid=int(e.expid),doppler=float(e.doppler),rv_error=float(e.velocity_error_kms)) for e in records]
    raw=dict(targetid=tid,arm='b',seed=seed,epoch=epochs,flux17=fluxes,trials=positions,
             range_rest_A=[REST_MIN,REST_MAX])
    row,levels,ledger=compact_rows(raw,meta,headers)
    row.update(status='measured',states_by_draw=''.join(r['states'] for r in ledger))
    for k,level in enumerate(levels):
        row[f'flux{k}_17']=level['flux17']
        row[f'recovered{k}']=level['recovered']
    return row


def run(windows,header_file,output,workers,n):
    windows=Path(windows)
    with open(header_file) as f:
        header={(int(r['expid']),r['camera']):float(r['heliocor']) for r in csv.DictReader(f)}
    files=sorted(windows.glob('*.npz'))
    if not files:
        raise ValueError('no target cutouts')
    jobs=((path,header,n) for path in files)
    if workers==1:
        rows=list(map(star,jobs))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows=list(pool.map(star,jobs,chunksize=4))
    keys=list(dict.fromkeys(k for r in rows for k in r))
    output=Path(output)
    output.parent.mkdir(exist_ok=True,parents=True)
    tmp=output.with_suffix('.partial')
    with tmp.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0) as binary:
            with io.TextIOWrapper(binary,encoding='utf-8',newline='') as text:
                writer=csv.DictWriter(text,fieldnames=keys,lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
    tmp.replace(output)
    return len(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('windows')
    p.add_argument('headers')
    p.add_argument('output')
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--trials',type=int,default=24)
    a=p.parse_args()
    print(run(a.windows,a.headers,a.output,a.workers,a.trials),flush=True)

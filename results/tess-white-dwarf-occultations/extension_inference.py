"""Parallel, restartable injection/recovery and per-trigger image measurements on the frozen extension."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import wd_search
from ffi_images import describe
from recovery import load_visits,occurrence,run_star_cell


def process(job):
    curve_path,cache,out=job
    data=json.loads(curve_path.read_text())
    if data['status']!='included':
        return
    source=data['source_id']
    summary_path=out/(source+'.json')
    if summary_path.is_file():
        return
    row=dict(path=data['lc_file'],cadence_s=data['cadence_s'],light_fraction=data['target_fraction'])
    visits=load_visits([row],curve_path.parent)
    # Check input identity again at the boundary to inference.
    if hashlib.sha256((curve_path.parent/row['path']).read_bytes()).hexdigest()!=data['lc_sha256']:
        raise ValueError('Corrupted photometric projection: '+source)
    sums, trials = [], []
    for p in [.25,1.40794]:
        for r in [1.,11.2]:
            summary, events=run_star_cell(source,visits,[data['search']],p,r,data['mass_solar'],
                                         data['radius_solar'],50,337109)
            sums.append(summary)
            trials.extend(events)
    with (out/(source+'.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,trials[0],lineterminator='\n')
        writer.writeheader()
        writer.writerows(trials)
    if data['search']['periodic'] or data['search']['isolated']:
        vet=describe(data,wd_search,cache,curve_path.parent)
        (out/(source+'-images.json')).write_text(json.dumps(vet,allow_nan=False)+'\n')
    summary_path.write_text(json.dumps(sums,allow_nan=False)+'\n')
    print('INJECTED',source,[(r['recovered'],r['possible_host']) for r in sums],flush=True)


def aggregate(curves,work,record):
    with (record/'extension_visits.csv').open() as f:
        stars=sorted(r['target'] for r in csv.DictReader(f))
    sums,trials,images=[],[],[]
    for source in stars:
        sums+=json.loads((work/(source+'.json')).read_text())
        with (work/(source+'.csv')).open() as f:
            entries=list(csv.DictReader(f))
        for r in sums[-4:]:
            relevant=[t for t in entries if (float(t['period_days']),float(t['radius_earth']))==(r['period_days'],r['radius_earth'])]
            assert sorted(int(t['trial']) for t in relevant)==list(range(50))
            assert sum(int(t['recovered']) for t in relevant)==r['recovered']
        trials+=entries
        path=work/(source+'-images.json')
        if path.is_file():
            images.append(json.loads(path.read_text()))
    if len(sums)!=4*len(set(stars)) or not sums:
        raise ValueError('Missing steady-state target set')
    for name,rows in [('extension_recovery_stars.csv',sums),('extension_recovery_trials.csv.gz',trials)]:
        buf=io.StringIO()
        writer=csv.DictWriter(buf,rows[0],lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
        content=buf.getvalue().encode()
        (record/name).write_bytes(gzip.compress(content,mtime=0) if name.endswith('.gz') else content)
    (record/'extension_limits.json').write_text(json.dumps(occurrence(sums),indent=2)+'\n')
    (record/'extension_image_checks.json.gz').write_bytes(gzip.compress(json.dumps(images,allow_nan=False).encode(),mtime=0))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('curves',type=Path)
    parser.add_argument('cache',type=Path)
    parser.add_argument('work',type=Path)
    parser.add_argument('record',type=Path)
    parser.add_argument('--workers',type=int,default=5)
    parser.add_argument('--partial',action='store_true',help='process currently complete acquisitions and defer global limits')
    args=parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    jobs=[(p,args.cache,args.work) for p in sorted(args.curves.glob('[0-9]*.json'))]
    with ProcessPoolExecutor(args.workers) as pool:
        list(pool.map(process,jobs))
    if not args.partial:
        aggregate(args.curves,args.work,args.record)

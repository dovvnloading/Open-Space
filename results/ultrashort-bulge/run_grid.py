"""Screen downloaded curves; intermediate checkpoints are reusable JSON records."""
import argparse
import csv
import json
import multiprocessing as mp
import tarfile
from pathlib import Path

import numpy as np
import fit

HERE = Path(__file__).resolve().parent


def records(path):
    with open(path) as handle:
        yield from csv.DictReader(handle)


def inputs(cache):
    for row in records(HERE/'ogle/catalog.csv'):
        name=row['event']
        year,_,number=name.split('-')
        path=cache/'.ogle-cache/raw'/year/f'blg-{number}.dat'
        if path.is_file():
            yield 'OGLE-'+name, path, float(row['t0_hjd']), float(row['tau_days'])
    index=cache/'staging/kmt/index.jsonl'
    if index.is_file():
        for line in index.read_text().splitlines():
            r=json.loads(line)
            path=cache/'staging/kmt'/(r['event']+'.tar.gz')
            if path.is_file():
                yield r['event'],path,2400000+float(r['t0'] or 0),float(r['te'] or 0)
    for path in sorted((cache/'.survey-cache/curves').glob('*.npz')):
        meta=cache/'.survey-cache/records'/(path.stem+'.json')
        if meta.is_file():
            r=json.loads(meta.read_text())
            yield path.stem,path,float(r.get('page_t0_jd') or np.nan),float(r.get('page_te_d') or np.nan)


def read_curve(job):
    name,path,t0,te=job
    if name.startswith('OGLE-'):
        t,mag,emag,*_=fit.lens.read_phot(path)
        f,e=fit.lens.to_flux(mag,emag)
        return t+2450000.,f,e,np.zeros(len(t),int),False
    if name.startswith('KMT-'):
        streams=[]
        with tarfile.open(path) as tar:
            for member in tar.getmembers():
                if not member.isfile() or not member.name.endswith('_I.pysis'):
                    continue
                arr=np.loadtxt(tar.extractfile(member),usecols=(0,1,2),ndmin=2)
                if len(arr):
                    streams.append(arr)
        if not streams:
            raise ValueError('no I-band pySIS stream')
        t,df,e=np.concatenate(streams).T
        g=np.repeat(np.arange(len(streams)),[len(a) for a in streams])
        return t,-df,e,g,True
    with np.load(path) as inp:
        t,f,e=inp['data'].T
    return t,f,e,np.zeros(len(t),int),True


def work(job):
    name,path,t0,te=job
    try:
        year=int(next(part for part in name.split('-') if part.isdigit() and len(part)==4))
        t,f,e,g,diff=read_curve(job)
        res=fit.screen(t,f,e,g,t0,te,year,diff)
        return dict(event=name,status='ok',**res)
    except (ValueError,TypeError,IndexError,OSError,KeyError) as exc:
        return dict(event=name,status='failed',error=str(exc))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path.cwd())
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--export',type=Path)
    args=parser.parse_args()
    args.checkpoint.mkdir(parents=True,exist_ok=True)
    if args.export:
        rows=[json.loads(p.read_text()) for p in sorted(args.checkpoint.glob('*.json'))]
        cols=list(dict.fromkeys(k for row in rows for k in row))
        with args.export.open('w') as f:
            w=csv.DictWriter(f,cols);w.writeheader();w.writerows(rows)
        return
    jobs=[j for j in inputs(args.cache) if not (args.checkpoint/(j[0]+'.json')).is_file()]
    print('Queued',len(jobs),flush=True)
    fit.setup()
    with mp.get_context('fork').Pool(args.workers) as pool:
        for res in pool.imap_unordered(work,jobs,chunksize=1):
            path=args.checkpoint/(res['event']+'.json')
            tmp=path.with_suffix('.tmp')
            tmp.write_text(json.dumps(res,default=lambda v:v.item())+'\n')
            tmp.replace(path)


if __name__=='__main__':
    main()

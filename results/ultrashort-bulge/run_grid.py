"""Screen downloaded curves; intermediate checkpoints are reusable JSON records."""
import argparse
import csv
import json
import multiprocessing as mp
import tarfile
from pathlib import Path

import numpy as np
import fit
from kmt.flux import baseline_caps

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
            if not path.is_file():
                path=cache/'staging/kmt/fallback'/(r['event']+'.tar.gz')
            if path.is_file():
                def estimate(v):
                    try:
                        value=float(v)
                        return value if value>0 else np.nan
                    except (TypeError,ValueError):
                        return np.nan
                yield r['event'],path,2400000+estimate(r['t0']),estimate(r['te'])
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
        fallback=path.parent.name=='fallback'
        suffix='_I.diapl' if fallback else '_I.pysis'
        streams=[]
        with tarfile.open(path) as tar:
            for member in tar.getmembers():
                if not member.isfile() or not member.name.endswith(suffix):
                    continue
                arr=np.loadtxt(tar.extractfile(member),usecols=(0,1,2),ndmin=2)
                if len(arr):
                    streams.append(arr)
        if not streams:
            raise ValueError('no I-band stream')
        t,df,e=np.concatenate(streams).T
        g=np.repeat(np.arange(len(streams)),[len(a) for a in streams])
        if fallback:
            return t+2400000.,df,e,g,True
        return t,-df,e,g,True
    with np.load(path) as inp:
        t,f,e=inp['data'].T
    return t,f,e,np.zeros(len(t),int),not name.startswith('Gaia')


def work(item):
    job, bound_kmt = item
    name,path,t0,te=job
    try:
        year=(2000+int(name[4:6]) if name.startswith('Gaia') else int(name.split('-')[1]))
        t,f,e,g,diff=read_curve(job)
        caps=None
        if bound_kmt and name.startswith('KMT') and path.parent.name!='fallback':
            raw=np.asarray(list(baseline_caps(path).values()))
            good=(np.isfinite(t)&(t>2450000)&(t<2470000)&np.isfinite(f)&np.isfinite(e)&(e>0))
            caps=np.where(np.isfinite(raw),1.15*raw,1e30)[np.unique(g[good])]
        res=fit.screen(t,f,e,g,t0,te,year,diff,caps)
        return dict(event=name,status='ok',**res)
    except (ValueError,TypeError,IndexError,OSError,KeyError) as exc:
        return dict(event=name,status='failed',error=str(exc))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path.cwd())
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--kmt-bound',action='store_true',help='KMT-only sensitivity with a cap of 1.15 times the reference baseline')
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
    jobs=[(j,args.kmt_bound) for j in inputs(args.cache) if (not args.kmt_bound or j[0].startswith('KMT')) and not (args.checkpoint/(j[0]+'.json')).is_file()]
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

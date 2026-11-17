"""Screen downloaded curves; intermediate checkpoints are reusable JSON records."""
import argparse
import csv
import json
import gzip
import platform
import sys
import multiprocessing as mp
import tarfile
from pathlib import Path

import numpy as np
import fit
from kmt.flux import baseline_caps
from replay_guard import external, fresh, verified, sha, identity,verify_receipt

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
    parser.add_argument('--grid-library',type=Path,required=True)
    parser.add_argument('--export',type=Path)
    args=parser.parse_args()
    args.cache=external(args.cache)
    if args.export:
        raise ValueError('legacy in-place checkpoint export disabled; use a verified new output')
    for name in ('ogle/catalog.csv','kmt/catalog-2024.tsv','kmt/catalog-2025.tsv'):verify_receipt(name)
    jobs=[(j,args.kmt_bound) for j in inputs(args.cache) if (not args.kmt_bound or j[0].startswith('KMT'))]
    verify_receipt('tables/grid-fits.csv.gz')
    with gzip.open(HERE/'tables/grid-fits.csv.gz','rt',newline='') as f:
        frozen={r['event']:r['input_sha256'] for r in csv.DictReader(f)}
    expected={name for name in frozen if not args.kmt_bound or name.startswith('KMT-')}
    if {j[0] for j,_ in jobs}!=expected:raise ValueError('incomplete or changed frozen event membership')
    target=external(args.checkpoint)
    if target==args.cache or args.cache in target.parents or target in args.cache.parents:raise ValueError('input/output overlap')
    for job,_ in jobs:
        name,path,_,_=job
        if not (name.startswith('OGLE-') or name.startswith('KMT-')):
            raise ValueError('derived survey NPZ cannot be bound to frozen raw bytes: '+name)
        if name not in frozen:raise ValueError('event absent from frozen receipt: '+name)
        verified(path,frozen[name])
    args.grid_library=external(args.grid_library)
    codepaths=[HERE/'run_grid.py',HERE/'fit.py',HERE/'grid.cpp',args.grid_library,HERE/'kmt/flux.py',HERE/'replay_guard.py',HERE/'ogle/catalog.csv',HERE/'kmt/catalog-2024.tsv',HERE/'kmt/catalog-2025.tsv',HERE.parent/'ogle-dark-lenses/lens.py']
    code_hashes={str(p):sha(p) for p in codepaths}
    runtime={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'vbm':getattr(fit.vbm,'__version__','unknown')}
    if 'fork' not in mp.get_all_start_methods():raise RuntimeError('full grid requires Linux/WSL fork')
    args.checkpoint=fresh(args.checkpoint)
    (args.checkpoint/'campaign.json').write_text(json.dumps({'model':'frozen grid','cap':args.kmt_bound,'code_sha256':code_hashes,'runtime':runtime},indent=2)+'\n')
    identities={j[0]:identity(j[1],j[0],bound,code_hashes,runtime,{'t0':float(j[2]) if np.isfinite(j[2]) else None,'te':float(j[3]) if np.isfinite(j[3]) else None}) for j,bound in jobs}
    print('Queued',len(jobs),flush=True)
    fit.setup(args.grid_library)
    with mp.get_context('fork').Pool(args.workers) as pool:
        for res in pool.imap_unordered(work,jobs,chunksize=1):
            path=args.checkpoint/(res['event']+'.json')
            tmp=path.with_suffix('.tmp')
            res['_identity']=identities[res['event']]
            tmp.write_text(json.dumps(res,default=lambda v:v.item())+'\n')
            tmp.replace(path)


if __name__=='__main__':
    main()

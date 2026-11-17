"""Run equal-window continuous alternatives for well-sampled short/large-source leads."""
import argparse
import json
import multiprocessing as mp
from pathlib import Path
from run_grid import inputs
from peak_models import solve


def select(r):
    if 'error' in r or r['peak_n']<8 or r['peak_streams']<2:return False
    if r['fs_te']<.5 or r['ps_te']<.5:return True
    return (r['rho']>1 and r['ps_chi2']-r['fs_chi2']>25 and
            r['flare_chi2']-r['fs_chi2']>0 and r['peak_chi2']/r['peak_n']<5)


def run(item):
    job,seed,width,bin_width,out=item
    try:
        r=solve(job,seed,1.15,5,True,width,bin_width)
    except (ArithmeticError,IndexError,ValueError) as err:
        r={'event':job[0],'error':str(err)}
    p=out/(job[0]+'.json');tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(r,allow_nan=False)+'\n');tmp.replace(p)
    return job[0]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--cache',type=Path,default=Path.cwd())
    a=parser.parse_args();root=Path(__file__).resolve().parent;out=a.cache/'staging/comparisons';out.mkdir(exist_ok=True)
    jobs={j[0]:j for j in inputs(a.cache)};todo=[]
    for p in (a.cache/'staging/continuous').glob('*.json'):
        if (out/p.name).exists() or (root/'peaks'/p.name).exists():continue
        r=json.loads(p.read_text())
        if not select(r) or p.stem not in jobs:continue
        seed=json.loads((a.cache/'staging/bound'/p.name).read_text())
        # Grid seed is common to all competing families. Bins shrink toward
        # short events; a maximum half-window of 100 d allows baseline support.
        effective=min(100.,max(4.,3.5*r['fs_te']*max(1,r['rho'])))
        bin_width=0 if effective<=4 else min(.4,effective/160)
        todo.append((jobs[p.stem],seed,effective,bin_width,out))
    todo.sort(key=lambda x:x[2])
    print('Queued',len(todo),flush=True)
    with mp.Pool(a.workers) as pool:list(pool.imap_unordered(run,todo))

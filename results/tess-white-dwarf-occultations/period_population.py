"""Secondary log-uniform period population on the union of the two disjoint FFI cohorts."""
import argparse
import csv
import gzip
import io
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from recovery import geometry,inject,load_manifest,load_visits,matches,occurrence,trial_seed
from wd_search import search

P_LO,P_HI,RADIUS=.25,3.,11.2
N_TRIALS,SEED=50,909563


def run(job):
    target,rows,root=job
    visits=load_visits(rows,root)
    if len(visits)!=1:
        raise ValueError('Secondary cohort requires exactly one visit per independent star')
    visit=visits[0]
    t,flux,_,dt,_=visit
    mass,rstar=float(rows[0]['mass_solar']),float(rows[0]['radius_solar'])
    def transit_probability(p):
        a,k,_,_=geometry(p,RADIUS,mass,rstar,0.)
        return (1+k)/a
    u0,u1=P_LO**(-2/3),P_HI**(-2/3)
    gbar=transit_probability(1.)*1.5*(u0-u1)/np.log(P_HI/P_LO)
    # The upper bound on total contact duration across the interval is at the
    # upper period endpoint for these detached orbits. Also test the lower endpoint.
    wmax=max(geometry(p,RADIUS,mass,rstar,0.)[-1] for p in [P_LO,P_HI])
    min_cycles=np.floor(np.ptp(t)/P_HI)
    edge=np.inf if min_cycles<=0 else (wmax+2*dt/86400.)/min_cycles
    candidates=[x for x in search(t,flux,dt)['periodic'] if P_LO-edge<=x['period_days']<=P_HI+edge]
    rng=np.random.default_rng(trial_seed(target+'|continuous',P_LO,RADIUS,SEED))
    trials=[]
    for j in range(N_TRIALS):
        p=(u0+rng.random()*(u1-u0))**(-1.5)
        epoch=t[0]+rng.random()*p
        bfrac=rng.random()
        a,k,b,width=geometry(p,RADIUS,mass,rstar,bfrac)
        out=search(t,inject(visit,p,epoch,a,k,b,width),dt)
        found=next((x for x in out['periodic'] if matches(x,visit,p,epoch,width)),None)
        trials.append(dict(target=target,trial=j,period_days=p,epoch_btjd=epoch,impact_fraction=bfrac,
                           duration_min=width*1440,recovered=int(found is not None),
                           recovered_period=found['period_days'] if found else '',
                           recovered_epoch=found['epoch_btjd'] if found else ''))
    summary=dict(target=target,period_days='0.25-3 log-uniform',radius_earth=RADIUS,
                 geometry_probability=float(gbar),trials=N_TRIALS,
                 recovered=sum(x['recovered'] for x in trials),possible_host=int(bool(candidates)),
                 period_edge_buffer_days=float(edge),compatible_periods_days=[x['period_days'] for x in candidates])
    return summary,trials


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('primary_visits',type=Path)
    parser.add_argument('extension_visits',type=Path)
    parser.add_argument('out',type=Path)
    parser.add_argument('--workers',type=int,default=5)
    args=parser.parse_args()
    jobs=[]
    cohorts={}
    for label,src in [('primary',args.primary_visits),('extension',args.extension_visits)]:
        for target,rows in load_manifest(src).items():
            if target in cohorts:
                raise ValueError('Cohorts overlap: '+target)
            cohorts[target]=label
            jobs.append((target,rows,src.parent))
    args.out.mkdir(parents=True,exist_ok=True)
    sums,trials=[],[]
    with ProcessPoolExecutor(args.workers) as pool:
        for s,t in pool.map(run,sorted(jobs,key=lambda x:x[0])):
            s['cohort']=cohorts[s['target']]
            sums.append(s)
            trials.extend(t)
    buf=io.StringIO()
    writer=csv.DictWriter(buf,trials[0],lineterminator='\n')
    writer.writeheader()
    writer.writerows(trials)
    (args.out/'population_trials.csv.gz').write_bytes(gzip.compress(buf.getvalue().encode(),mtime=0))
    (args.out/'population_stars.json').write_text(json.dumps(sums,indent=1)+'\n')
    limit=occurrence(sums)[0]
    limit.pop('period_days')
    limit.update(period_min_days=P_LO,period_max_days=P_HI,host_log_period_distribution='uniform',
                 transit_conditioned_period_density='proportional to P^(-5/3)',
                 possible_host_ids=[s['target'] for s in sums if s['possible_host']])
    (args.out/'population_limit.json').write_text(json.dumps(limit,indent=2)+'\n')
    print(limit)


if __name__=='__main__':
    main()

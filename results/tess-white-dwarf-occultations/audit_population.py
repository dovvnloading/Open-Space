"""Independently reconcile continuous-period draws, geometry integrals and count confidence limits."""
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.stats import beta


def check(record, manifests):
    rows={}
    for manifest in manifests:
        with manifest.open() as stream:
            for row in csv.DictReader(stream):
                assert row['target'] not in rows
                rows[row['target']]=(row,manifest.parent)
    stars=json.loads((record/'population_stars.json').read_text())
    limit=json.loads((record/'population_limit.json').read_text())
    assert len(stars)==len(rows)==limit['targets']
    trials=defaultdict(list)
    with gzip.open(record/'population_trials.csv.gz','rt') as stream:
        for r in csv.DictReader(stream):
            trials[r['target']].append(r)
    estimated=[]; lower=[]; geometries=[]
    for star in stars:
        target=star['target']
        row,root=rows[target]
        with np.load(root/row['path']) as lc:
            t=lc['time']
            start,end=float(t[0]),float(t[-1])
        dt=float(row['cadence_s'])/86400
        mass=float(row['mass_solar']); rs=float(row['radius_solar'])*6.957e8
        def transit_prob(p):
            return (rs+11.2*6.3781e6)/(1.3271244e20*mass*(p*86400/(2*math.pi))**2)**(1/3)
        g=quad(lambda p:transit_prob(p)/(p*math.log(12)),.25,3)[0]
        assert math.isclose(g,star['geometry_probability'],rel_tol=1e-11)
        examples=trials[target]
        assert sorted(int(x['trial']) for x in examples)==list(range(50))
        k=sum(int(x['recovered']) for x in examples)
        assert k==star['recovered'] and len(examples)==star['trials']
        for x in examples:
            p,e,w=float(x['period_days']),float(x['epoch_btjd']),float(x['duration_min'])/1440
            assert .25<=p<=3. and start<=e<start+p and 0<=float(x['impact_fraction'])<=1
            assert int(x['recovered'])==int(bool(x['recovered_period']))
            if x['recovered_period']:
                pp,ee=float(x['recovered_period']),float(x['recovered_epoch'])
                t1=e+round((start-e)/p)*p; t2=e+round((end-e)/p)*p
                k1=round((t1-ee)/pp); nc=round((t2-t1)/p)
                assert abs(ee+k1*pp-t1)<=.5*w+dt+1e-10
                assert abs(ee+(k1+nc)*pp-t2)<=.5*w+dt+1e-10
        estimated.append(k/50);geometries.append(g)
        lower.append(beta.ppf(.025/len(stars),k,51-k) if k else 0.)
    assert sum(s['recovered'] for s in stars)==limit['recovered']
    assert sum(s['trials'] for s in stars)==limit['trials']
    assert set(limit['possible_host_ids'])=={s['target'] for s in stars if s['compatible_periods_days']}
    count=limit['possible_hosts']
    assert count==len(limit['possible_host_ids'])
    assert math.isclose(sum(estimated),limit['effective_transiting_stars'])
    assert math.isclose(np.dot(estimated,geometries),limit['geometry_weighted_stars'])
    for field,factor in [('transiting_host_upper95',1.),('intrinsic_host_upper95',np.array(geometries))]:
        f=limit[field]; probs=np.array(lower)*factor*(1 if f is None else f)
        odds=probs/(1-probs)
        coef=np.array([1.]+[0.]*count)
        for rate in odds:
            coef[1:]+=rate*coef[:-1]
        tail=math.prod(1-probs)*coef.sum()
        assert (tail>=.025-1e-12) if f is None else math.isclose(tail,.025,rel_tol=1e-8)
    print(f"Verified {limit['trials']} secondary draws, all recovered endpoint matches, geometry and both confidence limits.")


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    parser.add_argument('manifests',nargs=2,type=Path)
    args=parser.parse_args()
    check(args.record,args.manifests)

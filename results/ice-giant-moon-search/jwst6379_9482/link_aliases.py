"""Evaluate shared circular orbits across deep visits by enumerating integer revolutions."""
import argparse
import csv
from itertools import product
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'jwst9482'))
from link_visits import make_cube,evaluate,deep_score,trials,geometry,search


def compatible(a,b,dt,phase_step,tol):
    # Both peak radii have finite grid precision and astrometric error.
    lo=max(a['a_km']-tol,b['a_km']-tol,54500)
    hi=min(a['a_km']+tol,b['a_km']+tol,86000)
    if lo>=hi:
        return
    for u,v in product([-.5,0,.5],repeat=2):
        phi1=(a['phase_deg']+phase_step*u)%360
        phi2=(b['phase_deg']+phase_step*v)%360
        dphi=(phi2-phi1)/360
        def spins(r):
            return geometry.frequency('Uranus',r)*dt/(2*np.pi)
        for n in range(int(np.ceil(spins(hi)-dphi)),int(np.floor(spins(lo)-dphi))+1):
            r=brentq(lambda x:spins(x)-(n+dphi),lo,hi)
            yield r,phi1,phi2,n


def read_peaks(path):
    with open(path/'peaks.csv') as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['early','later','early_peaks','later_peaks','ee','out']:
        p.add_argument(name,type=Path)
    p.add_argument('--radius-tolerance',type=float,default=375.)
    a=p.parse_args()
    cube_paths=(a.early,a.later)
    sets=(a.early_peaks,a.later_peaks)
    c=[make_cube(f,a.ee,True) for f in cube_paths]
    ks=[search.validate_cube(x,min_visits=1) for x in c]
    cuts=search.Cuts(psf_support=.8)
    maps=[search.prepare(x,k,cuts) for x,k in zip(c,ks,strict=True)]
    params=[json.loads((f/'parameters.json').read_text()) for f in sets]
    t0,t1=(x['epoch_jd_tdb'] for x in params)
    dt=(t1-t0)*86400
    early,late=(read_peaks(path) for path in sets)
    a.out.mkdir(parents=True,exist_ok=False)
    fields=['early_peak','late_peak','early_sign','late_sign','a_km','phase_first_deg','phase_second_deg','whole_revs',
            'z_first','z_second','first_early','first_late','second_early','second_late','flux_first_Jy','flux_second_Jy']
    summary={}
    with TemporaryDirectory() as tmp,open(a.out/'passing_pairs.csv','x',newline='') as hits:
        tmp=Path(tmp);writer=csv.DictWriter(hits,fieldnames=fields);writer.writeheader()
        def execute(group):
            elem=np.zeros((len(group),6))
            elem[:,0]=[v['a_km'] for v in group]
            elem[:,-1]=[v['phase_first_deg'] for v in group]
            z=[]; e=[];l=[];flux=[];rank=[]
            for j in range(2):
                bank=trials(tmp/'e.csv',cube_paths[j],tmp/'tracks.npz',elem,t0)
                n,i=evaluate(c[j],bank,ks[j],maps[j],cuts)
                signs=np.array([v['early_sign' if j==0 else 'late_sign'] for v in group])
                rr,zz,aa,bb,_=deep_score(n*signs[:,None],i)
                rank.append(rr);z.append(zz);e.append(aa);l.append(bb)
                flux.append(np.divide(n.sum(1),i.sum(1),out=np.zeros(len(n)),where=i.sum(1)>0))
            good=(rank[0]>=1)&(rank[1]>=1)&(z[0]>=7)&(z[1]>=7)
            for k in np.flatnonzero(good):
                row=group[k]
                row.update(z_first=z[0][k],z_second=z[1][k],first_early=e[0][k],first_late=l[0][k],
                           second_early=e[1][k],second_late=l[1][k],flux_first_Jy=flux[0][k],flux_second_Jy=flux[1][k])
                writer.writerow(row)
                key=f"{row['early_sign']:+d}/{row['late_sign']:+d}"
                summary[key]['passing_aliases']+=1
            hits.flush()
        block=[]
        for j,p1 in enumerate(early):
            for k,p2 in enumerate(late):
                sign1,sign2=int(p1['sign']),int(p2['sign'])
                key=f'{sign1:+d}/{sign2:+d}'
                stats=summary.setdefault(key,dict(peak_pairs=0,trial_aliases=0,passing_aliases=0))
                arr=compatible(p1,p2,dt,.75,a.radius_tolerance)
                count=0
                for r,ph1,ph2,n in arr:
                    block.append(dict(early_peak=j,late_peak=k,early_sign=sign1,late_sign=sign2,
                                      a_km=r,phase_first_deg=ph1,phase_second_deg=ph2,whole_revs=n))
                    count+=1
                    if len(block)>=500:
                        execute(block);block=[]
                stats['trial_aliases']+=count
                stats['peak_pairs']+=int(count>0)
        if block:execute(block)
    (a.out/'parameters.json').write_text(json.dumps(dict(epoch1_tdb_jd=t0,epoch2_tdb_jd=t1,days=t1-t0,
        radius_tolerance_km=a.radius_tolerance,signed_counts=summary,
        source_hashes={str(p):search.file_hash(p) for p in [a.early,a.later,a.ee,a.early_peaks/'peaks.csv',a.later_peaks/'peaks.csv',Path(__file__)]}),
        indent=2)+'\n')


if __name__=='__main__':
    main()

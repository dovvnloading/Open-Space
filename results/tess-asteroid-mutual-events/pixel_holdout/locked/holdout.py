"""Score held-out observations once, using an event phase learned on earlier stamps."""
import argparse
import csv
import numpy as np
from pathlib import Path
from score import rot, depth

def quantify(rows, lock, flux_index):
    period=float(lock['period_h']); duration=float(lock['duration_h'])
    pspin={29489:3.1679,3177:3.547245}[int(lock['mpc'])]
    t0=float(lock['training_jd_start'])
    hold=rows[rows[:,0]>float(lock['holdout_jd_after'])]
    # Fixed Gaussian profile is primary; 1.65 pixel aperture is secondary.
    # No search over event or spin period, duration, phase, or extraction width.
    f=hold[:,flux_index]
    use=(np.isfinite(f))&(f>0)&(hold[:,7]<np.quantile(hold[:,7],.97))
    t,y=(hold[use,0]-t0)*24,-2.5*np.log10(f[use])
    trend=np.array([np.median(y[abs(t-x)<=15]) for x in t])
    y=y-trend
    r=rot(t,y,pspin)
    sigma=1.4826*np.median(abs(r-np.median(r)))
    keep=abs(r-np.median(r))<4*sigma
    t,r=t[keep],rot(t[keep],y[keep],pspin)
    ev=depth(t,r,period,duration,float(lock['center_phase']))
    if not len(ev): return [], []
    d=ev[:,1];rng=np.random.default_rng(20261115)
    lo,hi=np.quantile(np.median(rng.choice(d,(10000,len(d)),replace=True),axis=1),[.025,.975])
    summary=[lock['mpc'],lock['sector'],flux_index,t0+t.min()/24,t0+t.max()/24,len(t),
             1.4826*np.median(abs(r-np.median(r))),len(d),sum(d>0),np.median(d),lo,hi]
    table=[[lock['mpc'], flux_index, int(k), val] for k,val in ev]
    return summary,table

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    with (Path(__file__).parent/'phase_lock.csv').open() as f: locks=list(csv.DictReader(f))
    totals=[];events=[]
    for lock in locks:
        a=np.load(args.directory/f"stamps{lock['mpc']}_{lock['sector']}.npz")['rows']
        for j in [8,9]:
            result,table=quantify(a,lock,j);totals.append(result);events.extend(table)
    for filename,labels,table in [('holdout.csv','mpc sector extractor first_jd last_jd points residual_mad_mag events positive median_depth_mag bootstrap_low bootstrap_high',totals),
                                   ('events.csv','mpc extractor cycle depth_mag',events)]:
        with (args.out/filename).open('w',newline='') as f:
            w=csv.writer(f);w.writerow(labels.split());w.writerows(table)

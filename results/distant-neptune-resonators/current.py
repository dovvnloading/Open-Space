"""Validate both candidate resonances with the full, updated MPC Cartesian orbit fits."""
import csv
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.signal import find_peaks
from dynamics import WORK, YEAR, init, integrate, angle, arc

TARGETS=[('2017-DO171',54464562,8,1),('2022-GV6',54276384,29,3)]

def setup(which, n, step):
    name, seed, p, q = which
    fit = json.loads((WORK/(name+'-orbit.json')).read_text())[0]['mpc_orb'][0]
    epoch = fit['epoch_data']['epoch']+2400000.5
    sim,sun = init(epoch,step)
    car=fit['CAR']
    mean=np.array(car['coefficient_values'])
    cov=np.array([[car['covariance'][f'cov{min(i,j)}{max(i,j)}'] for j in range(6)] for i in range(6)])
    sig=np.sqrt(cov.diagonal())
    chol=np.linalg.cholesky(cov / np.outer(sig,sig))
    samples=mean+(np.random.default_rng(seed+1).standard_normal((n,6))@chol.T)*sig
    for state in np.vstack([mean,samples]):
        pos=np.array(sun.xyz)+state[:3]
        vel=np.array(sun.vxyz)+state[3:]*YEAR
        sim.add(x=pos[0],y=pos[1],z=pos[2],vx=vel[0],vy=vel[1],vz=vel[2])
    return sim

def run(task):
    which,n,dt,years,cadence,method=task
    name,_,p,q=which
    sim=setup(which,n,dt)
    sim.integrator=method
    times,data=integrate(sim,years,cadence)
    phi=np.array([np.rad2deg(np.unwrap(angle(data[:,i],p,q))) for i in range(n+1)])
    # Subtract each clone's initial angle; extrema can be combined across directions.
    relative=phi-phi[:,0,None]
    rows=[]
    for i in range(n+1):
        width,center=arc(np.deg2rad(phi[i]))
        first=np.flatnonzero(np.maximum.accumulate(relative[i])-np.minimum.accumulate(relative[i])>=360.)
        peaks,_=find_peaks(phi[i],prominence=25,distance=max(1,int(8000/cadence)))
        periods=np.diff(peaks)*cadence
        rows.append([name,method,n,years,dt,cadence,i,round(width,4),round(center,4),
                     round(min(relative[i]),4),round(max(relative[i]),4),int(abs(times[first[0]])) if len(first) else '',
                     len(peaks),int(np.median(periods)) if len(periods) else ''])
    path=WORK/f'current_{name}_{method}_{n}_{dt}_{years}.csv'
    with path.open('w') as out:
        wr=csv.writer(out);wr.writerow(['object','method','draws','years','step_yr','cadence_yr','clone','arc_deg','center_deg','relative_min_deg','relative_max_deg','first_360_yr','n_peaks','peak_spacing_yr']);wr.writerows(rows)
    if n==0 and dt>0:
        with (WORK/f'trace_{name}_{method}.csv').open('w') as out:
            wr=csv.writer(out);wr.writerow(['years','phi_deg','barycentric_a_au'])
            wr.writerows((int(t),round(v%360,5),round(a,6)) for t,v,a in zip(times,phi[0],data[:,0,0]))
    print(path.name,round(rows[0][7],1),sum(r[10]-r[9]<355 for r in rows[1:]),n,flush=True)

if __name__=='__main__':
    settings=[(128,.5,400000,100,'whfast'),(128,-.5,-400000,100,'whfast'),
              (64,.5,2000000,200,'whfast'),(64,-.5,-2000000,200,'whfast'),
              (8,.25,2000000,200,'whfast'),(8,-.25,-2000000,200,'whfast'),
              (8,.5,10000000,250,'whfast'),(8,-.5,-10000000,250,'whfast'),
              (0,.25,100000,50,'whfast'),(0,.25,100000,50,'ias15')]
    with ProcessPoolExecutor(7) as pool:
        list(pool.map(run, [(o,*s) for s in settings for o in TARGETS]))

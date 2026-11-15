"""Measure transit-phase drift and propagate coplanar Jacobi TTV posterior draws."""
import argparse
import csv
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import rebound

HERE = Path(__file__).resolve().parent
JE = 332946.0487  # Solar mass in Earth masses.
SYSTEMS = {
    2092: ((2092.02, 2092.01, 2092.03), np.array([1, -6, 5])),
    2086: ((2086.01, 2086.02, 2086.03), np.array([1, -2, 1])),
}


def fit_transits():
    times=np.genfromtxt(HERE/'transits.csv',delimiter=',',names=True)
    results=[]
    for host,(ids,coeff) in SYSTEMS.items():
        if not set(ids).issubset(set(times['koi'])):
            continue  # Holczer tabulates only planet c of the Kepler-60 control.
        for lo,hi in [(60,1550),(60,805),(805,1550)]:
            phase=rate=phasevar=ratevar=0.
            counts=[]
            for planet,k in zip(ids,coeff):
                rows=times[(times['koi']==planet)&(times['out']==0)&(times['over']==0)&(times['error_min']>0)&(times['linear_bjd_minus_2454900']>=lo)&(times['linear_bjd_minus_2454900']<=hi)]
                n=rows['n']; count=len(n);counts.append(count)
                t=rows['linear_bjd_minus_2454900']+rows['oc_min']/1440
                n0=np.median(n);matrix=np.column_stack([np.ones(count),n-n0])
                weight=(1440/rows['error_min'])**2
                cov=np.linalg.inv(matrix.T@(matrix*weight[:,None]))
                epoch,period=cov@(matrix.T@(weight*t))
                residual=t-matrix@[epoch,period]
                cov*=max(1., np.dot(weight,residual**2)/(count-2))
                offset=(lo+hi)/2-epoch
                phase+=360*k*(n0+offset/period)
                rate+=360*k/period
                grad=360*k*np.array([-1/period,-offset/period**2])
                phasevar+=grad@cov@grad
                ratevar+=(360*k/period**2)**2*cov[1,1]
            results.append([host,lo,hi,*counts,phase%360,math.sqrt(phasevar),rate*1470,math.sqrt(ratevar)*1470])
    with (HERE/'phase.csv').open('w') as f:
        writer=csv.writer(f,lineterminator='\n')
        writer.writerow(['koi','window_start','window_end','n_b','n_c','n_d','phase_deg','phase_se_deg','drift_deg_per_1470d','drift_se_deg_per_1470d'])
        writer.writerows(results)


def initialize(values, steps):
    sim=rebound.Simulation()
    sim.units=('AU','day','Msun')
    sim.add(m=1)
    pars=values.reshape(-1,5)
    for ratio,period,esinw,ecosw,ttrans in pars:
        ecc=np.hypot(esinw,ecosw)
        omega=np.arctan2(esinw,ecosw)
        ftrans=np.pi/2-omega
        eccanom=2*np.arctan2(np.sqrt(1-ecc)*np.sin(ftrans/2),np.sqrt(1+ecc)*np.cos(ftrans/2))
        mean=eccanom-ecc*np.sin(eccanom)-2*np.pi*(ttrans-780)/period
        sim.add(m=ratio/JE,P=period,e=ecc,omega=omega,M=mean)
    sim.move_to_com()
    sim.integrator='whfast'
    sim.dt=min(pars[:,1])/steps
    # Stop immediately upon any approach inside the smallest initial adjacent mutual Hill radius.
    orbits=sim.orbits()
    radii=[(a.a+b.a)/2*((pars[i,0]+pars[i+1,0])/JE/3)**(1/3) for i,(a,b) in enumerate(zip(orbits,orbits[1:]))]
    sim.exit_min_distance=min(radii)
    sim.exit_max_distance=10*orbits[-1].a
    return sim


def propagate(task):
    host,idx,values,mode,steps,cadence,years=task
    sim=initialize(values,steps); e0=sim.energy()
    max_error=0.; phases=[];state='complete'
    coeff=SYSTEMS[host][1]
    for target in np.arange(0,years*365.25+cadence,cadence):
        try:
            sim.integrate(target,exact_finish_time=0)
        except (rebound.Escape,rebound.Encounter):
            state='close_or_escape'
            break
        max_error=max(max_error,abs(sim.energy()/e0-1))
        ph=sum(k*o.l for k,o in zip(coeff,sim.orbits()))
        phases.append(ph)
    unwrapped=np.rad2deg(np.unwrap(phases))
    span=np.ptp(unwrapped);winding=(unwrapped[-1]-unwrapped[0])/360
    step=np.max(np.abs(np.diff(unwrapped)),initial=0)
    center=((np.min(unwrapped)+np.max(unwrapped))/2)%360
    if state=='complete':
        if max_error>1e-4 or step>90:
            state='unresolved'
        elif span>=720:
            state='circulating' if abs(winding)>=1 else 'multi_turn_excursion'
        elif span<180:
            state='confined'
        else:
            state='intermediate'
    return [host,idx,mode,steps,cadence,round(sim.t/365.25,5),state,round(span,3),round(winding,5),round(step,3),round(center,3),f'{max_error:.6g}']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    fit_transits()
    rows=np.loadtxt(HERE/'posterior.csv',delimiter=',',skiprows=1)
    tasks=[]
    for host in SYSTEMS:
        samples=rows[rows[:,0]==host]
        for mode,steps,cadence,subset in [('nominal',24,45,samples),('refined',48,11.25,samples[:16])]:
            for row in subset:
                tasks.append((host,int(row[1]),row[2:],mode,steps,cadence,300))
    header=['koi','row_0based','run','steps_per_inner_orbit','sampling_days','duration_year','classification','span_deg','net_turns','max_sample_step_deg','arc_midpoint_deg','max_fractional_energy_error']
    lookup={(int(r[0]),int(r[1])):r[2:] for r in rows}
    def write_runs(name, jobs):
        with (HERE/name).open('w') as f, ProcessPoolExecutor(args.workers) as pool:
            writer=csv.writer(f,lineterminator='\n')
            writer.writerow(header)
            results=list(pool.map(propagate,jobs))
            writer.writerows(results)
        return results
    runs=write_runs('orbits.csv', tasks)
    def eligible(result):
        return result[6] not in ('circulating','close_or_escape')
    long_tasks=[(r[0],r[1],lookup[r[0],r[1]],'extended',48,19.5,3000) for r in runs
                if r[2]=='nominal' and r[0]==2092 and eligible(r)]
    extended=write_runs('extended.csv',long_tasks)
    deepest=[(r[0],r[1],lookup[r[0],r[1]],'long',48,45.,30000) for r in extended if eligible(r)]
    write_runs('long.csv',deepest)


if __name__=='__main__': main()

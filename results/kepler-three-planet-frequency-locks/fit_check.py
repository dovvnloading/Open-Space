"""Check the adopted orbital convention by reproducing Kepler-359 transit positions."""
import csv
import numpy as np
from analyze import HERE, initialize

rows=np.loadtxt(HERE/'posterior.csv',delimiter=',',skiprows=1)
measurements=np.genfromtxt(HERE/'transits.csv',delimiter=',',names=True)
data=[]
for event in measurements:
    host=event['koi']
    if host not in (2092.02,2092.01,2092.03) or event['out'] or event['over'] or event['error_min']<=0:
        continue
    planet={2092.02:1,2092.01:2,2092.03:3}[host]
    time=event['linear_bjd_minus_2454900'] + event['oc_min']/1440 - 780
    data.append((time,planet,event['error_min']))
with (HERE/'fit_check.csv').open('w') as f:
    out=csv.writer(f,lineterminator='\n')
    out.writerow(['row_0based','events','median_abs_residual_min','median_abs_normalized_residual'])
    for row in rows[rows[:,0]==2092][:16]:
        sim=initialize(row[2:],24)
        sim.integrator='ias15'
        dt=[]; z=[]
        for time,planet,sigma in sorted(data):
            sim.integrate(time)
            difference=sim.particles[planet]-sim.particles[0]
            if difference.y<0: raise ValueError('Predicted position is not a transit')
            # x/vx is the local first-order timing residual at the observed transit.
            delta=abs(difference.x/difference.vx)*1440
            dt.append(delta);z.append(delta/sigma)
        out.writerow([int(row[1]),len(dt),round(np.median(dt),4),round(np.median(z),4)])

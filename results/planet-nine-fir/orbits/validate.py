"""Numerical cross-check, grid convergence and historical centroid sensitivities."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import propagate as p

HERE = Path(__file__).resolve().parent


def bounds(s,t):
    ra,de = p.sky(s,t)
    return np.array([ra.min(),ra.max(),de.min(),de.max()])


def disturbed(coord, east, north):
    u = p.unit(*coord)
    ra = np.radians(coord[0])
    e = np.array([-np.sin(ra), np.cos(ra), 0])
    q = u + e * east/p.ARCSEC + np.cross(u,e) * north/p.ARCSEC
    q /= p.norm(q)
    return [float(np.degrees(np.arctan2(q[1],q[0])))%360, float(np.degrees(np.arcsin(q[2])))]


def main():
    inputs = json.loads((HERE/'orbit_input.json').read_text())
    t = p.epoch('2026-11-17')
    s = p.family(inputs)
    base = bounds(s,t)
    fine = p.family(inputs,17,41,81)
    dense = bounds(fine,t)
    def acceleration(t,y):
        return np.r_[y[3:], -p.MU*y[:3]/p.norm(y[:3])**3]
    diffs = []
    for j in np.linspace(0,len(s[0])-1,17,dtype=int):
        dt = t-s[2][j]
        y = np.r_[s[0][j],s[1][j]]
        sol = solve_ivp(acceleration,[0,dt],y,method="DOP853",rtol=2e-13,atol=1e-14)
        prec = p.kepler(s[0][j],s[1][j],dt)
        diffs.append(p.norm(sol.y[:3,-1]-prec)/p.norm(prec)*p.ARCSEC)
    pos1 = p.direction(s,p.epoch(inputs['iras_coadd_interval_utc'][0]))
    pos2 = p.direction(s,p.epoch(inputs['iras_coadd_interval_utc'][1]))
    sweep = 2*np.arcsin(p.norm(pos1-pos2)/2)*p.ARCSEC
    models = [p.family(inputs,coadd_endpoint_fraction=f) for f in (.5,1.)]
    cb = np.stack([base,*[bounds(m,t) for m in models]])
    cbounds = np.array([cb[:,0].min(),cb[:,1].max(),cb[:,2].min(),cb[:,3].max()])
    results=[]
    # Native IRAS 18 x 7 arcsec PA 79 deg. Four cardinal ellipse directions.
    # AKARI radius below is an explicit illustrative displacement, no catalog covariance.
    pa=np.radians(79.)
    iras_offsets=[(0,0),(18*np.sin(pa),18*np.cos(pa)),(-18*np.sin(pa),-18*np.cos(pa)),
                  (7*np.cos(pa),-7*np.sin(pa)),(-7*np.cos(pa),7*np.sin(pa))]
    for e1,n1 in iras_offsets:
        for e2,n2 in ((0,0),(40,0),(-40,0),(0,40),(0,-40)):
            changed = dict(inputs, iras_radec_deg=disturbed(inputs['iras_radec_deg'],e1,n1),
                           akari_radec_deg=disturbed(inputs['akari_radec_deg'],e2,n2))
            results.append(bounds(p.family(changed),t))
    results=np.stack(results)
    out = dict(validation_epoch_utc='2026-11-17', independent_integrator_max_position_error_as=max(diffs),
               nominal_bounds_deg=base.tolist(),fine_bounds_deg=dense.tolist(),fine_counts=fine[-1],
               coadd_mixture_union_bounds_deg=cbounds.tolist(),
               centroid_stress_union_bounds_deg=[results[:,0].min(),results[:,1].max(),results[:,2].min(),results[:,3].max()],
               modeled_iras_start_to_end_sweep_arcsec=[sweep.min(),sweep.max()])
    (HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()

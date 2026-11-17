"""Calibrate the native-cadence search by adding physical transits to real aperture flux."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta

from wd_search import flatten, search

# IAU nominal conversions, plus a fixed negligible companion mass and circular orbit.
SOLAR_GM = 1.3271244e20
SOLAR_RADIUS = 6.957e8
EARTH_RADIUS = 6.3781e6
NODES, WEIGHTS = np.polynomial.legendre.leggauss(32)


def geometry(period, radius_earth, mass, star_radius, impact_fraction):
    k = radius_earth * EARTH_RADIUS / (star_radius * SOLAR_RADIUS)
    a = np.cbrt(SOLAR_GM*mass*(period*86400/(2*np.pi))**2)/(star_radius*SOLAR_RADIUS)
    if min(mass, star_radius, period, radius_earth) <= 0 or a <= 1+k:
        raise ValueError('Invalid detached-orbit geometry')
    b = impact_fraction * (1+k)
    chord = max(0., (1+k)**2-b**2)
    width = period/np.pi * np.arcsin(np.sqrt(chord/(a*a-b*b)))
    return float(a), float(k), float(b), float(width)


def overlap(z, k):
    """Occulted fraction of a uniform stellar disk given projected separation/Rstar."""
    depth = np.zeros_like(z)
    depth[z <= abs(1-k)] = min(1., k*k)
    partial = (abs(1-k) < z) & (z < 1+k)
    d = z[partial]
    c1 = np.clip((d*d+1-k*k)/(2*d), -1., 1.)
    c2 = np.clip((d*d+k*k-1)/(2*d*k), -1., 1.)
    root = np.sqrt(np.maximum(0., (-d+1+k)*(d+1-k)*(d-1+k)*(d+1+k)))
    depth[partial] = (np.arccos(c1) + k*k*np.arccos(c2) - .5*root)/np.pi
    return np.clip(depth, 0., 1.)


def inject(visit, period, epoch, a, k, b, width):
    t, flux, baseline, cadence, fraction = visit
    exposure = cadence/86400.
    phase = (t-epoch+.5*period) % period-.5*period
    affected = abs(phase) < .5*(width+exposure)
    removed = np.zeros_like(flux)
    # Integrate just the intersection of this exposure and first-to-fourth contact.
    # This also resolves extremely short grazing transits within a long exposure.
    left = np.maximum(phase[affected]-.5*exposure, -.5*width)
    right = np.minimum(phase[affected]+.5*exposure, .5*width)
    mid = .5*(left+right)[:,None]
    half = .5*(right-left)[:,None]
    angle = (mid+half*NODES)*(2*np.pi/period)
    distance = np.hypot(a*np.sin(angle), b*np.cos(angle))
    removed[affected] = (overlap(distance, k) @ WEIGHTS) * (half[:,0]/exposure)
    return flux - removed*fraction*baseline


def matches(row, visit, period, epoch, width):
    t, _, _, dt, _ = visit
    p, e = row['period_days'], row['epoch_btjd']
    tol = .5*width + dt/86400.
    endpoints = (t[0],t[-1])
    centers = [epoch+np.rint((x-epoch)/period)*period for x in endpoints]
    initial = np.rint((centers[0]-e)/p)
    cycles = np.rint((centers[1]-centers[0])/period)
    return (abs(e+initial*p-centers[0]) <= tol and
            abs(e+(initial+cycles)*p-centers[1]) <= tol)


def load_manifest(path):
    by_star = defaultdict(list)
    with path.open() as stream:
        for row in csv.DictReader(stream):
            by_star[row['target']].append(row)
    if not by_star:
        raise ValueError('Empty visit manifest')
    return by_star


def load_visits(rows, root):
    visits = []
    for row in rows:
        with np.load(root/row['path']) as data:
            t, f = np.asarray(data['time'],float), np.asarray(data['flux'],float)
            ok = np.isfinite(t) & np.isfinite(f)
            if 'quality' in data:
                ok &= data['quality'] == 0
            t, idx = np.unique(t[ok], return_index=True)
            f = f[ok][idx]
        if len(t) < 20 or np.ptp(t) < .1:
            raise ValueError(f"Unusable manifest visit {row['path']}")
        base = flatten(t, f)
        valid = np.isfinite(base)
        if not valid.any():
            raise ValueError(f"No valid baseline in {row['path']}")
        # search() rejects these short segments identically before and after injection.
        base[~valid] = np.median(base[valid])
        fraction, dt = float(row['light_fraction']), float(row['cadence_s'])
        if not 0 < fraction <= 1 or dt <= 0:
            raise ValueError(f"Invalid aperture fraction or cadence in {row['path']}")
        visits.append((t, f, base, dt, fraction))
    return visits


def trial_seed(target, period, radius, seed):
    raw = f'{target}|{period:.12g}|{radius:.12g}|{seed}'.encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], 'big')


def run_star_cell(target, visits, baseline_results, period, radius, mass, star_radius, n, seed):
    rng = np.random.default_rng(trial_seed(target, period, radius, seed))
    ref = min(v[0][0] for v in visits)
    a,k,_,width_max = geometry(period, radius, mass, star_radius, 0)
    # For ANY unknown phase/impact to satisfy the endpoint recovery criterion, the
    # two ephemerides must accumulate <= twice the longest allowed endpoint error.
    maybe = False
    for visit, result in zip(visits, baseline_results):
        t, _, _, dt, _ = visit
        min_cycles = max(0., np.floor(np.ptp(t)/period))
        for row in result['periodic']:
            if min_cycles*abs(row['period_days']-period) <= width_max + 2*dt/86400.:
                maybe = True
    summary = dict(target=target, period_days=period, radius_earth=radius, mass_solar=mass,
                   star_radius_solar=star_radius, geometry_probability=(1+k)/a,
                   trials=n, recovered=0, possible_host=int(maybe))
    trials = []
    for j in range(n):
        epoch = ref + rng.random()*period
        impact_fraction = rng.random()
        a,k,b,width = geometry(period, radius, mass, star_radius, impact_fraction)
        detected = ''
        for visit_index, visit in enumerate(visits):
            t, _, _, dt, _ = visit
            simulated = inject(visit, period, epoch, a,k,b,width)
            result = search(t, simulated, dt)
            if any(matches(row, visit, period, epoch, width) for row in result['periodic']):
                detected = str(visit_index)
                break
        summary['recovered'] += int(bool(detected))
        trials.append(dict(target=target, period_days=period, radius_earth=radius, trial=j,
                           epoch_btjd=epoch, impact_fraction=impact_fraction, duration_min=width*1440,
                           recovered=int(bool(detected)), first_visit_index=detected))
    return summary, trials


def occurrence(summaries):
    """Pointwise joint 95% upper bounds on per-star host fraction for exact size and period."""
    groups = defaultdict(list)
    for row in summaries:
        groups[(row['period_days'], row['radius_earth'])].append(row)
    results = []
    for (p, r), rows in groups.items():
        n = len(rows)
        if len({row['target'] for row in rows}) != n:
            raise ValueError('Repeated star/cell in inference')
        q = np.array([row['recovered']/row['trials'] for row in rows])
        # Simultaneous one-sided lower confidence limits for independent Bernoulli trials.
        lo = np.array([beta.ppf(.025/n, row['recovered'], row['trials']-row['recovered']+1)
                       if row['recovered'] else 0. for row in rows])
        g = np.array([row['geometry_probability'] for row in rows])
        count = sum(row['possible_host'] for row in rows)
        def bound(prob):
            def tail(f):
                # P(K<=count), Poisson-binomial, dropping higher counts at every update.
                pmf = np.array([1.])
                for rate in prob*f:
                    pmf = np.convolve(pmf, [1-rate,rate])[:count+1]
                return pmf.sum()
            return float(brentq(lambda f:tail(f)-.025, 0., 1.)) if tail(1.) < .025 else None
        results.append(dict(period_days=p, radius_earth=r, targets=n, possible_hosts=count,
                       trials=sum(x['trials'] for x in rows), recovered=sum(x['recovered'] for x in rows),
                       effective_transiting_stars=float(sum(q)), geometry_weighted_stars=float(sum(q*g)),
                       transiting_host_upper95=bound(lo), intrinsic_host_upper95=bound(lo*g)))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path, help='CSV: target,path,cadence_s,light_fraction[,mass_solar,radius_solar]')
    parser.add_argument('out',type=Path)
    parser.add_argument('--periods',type=float,nargs='+',default=[.15, 1.4, 4.])
    parser.add_argument('--radii',type=float,nargs='+',default=[1., 4., 11.2])
    parser.add_argument('--trials',type=int,default=32)
    parser.add_argument('--seed',type=int,default=17563)
    parser.add_argument('--target', help='process only this star, for independent parallel shards')
    args = parser.parse_args()
    if args.trials < 1:
        parser.error('Positive --trials required')
    args.out.mkdir(parents=True,exist_ok=True)
    stars = load_manifest(args.manifest)
    select = [args.target] if args.target else list(stars)
    totals = []
    for target in select:
        visits = load_visits(stars[target], args.manifest.parent)
        baseline = [search(t,f,dt) for t,f,_,dt,_ in visits]
        params = {(float(row.get('mass_solar') or .6), float(row.get('radius_solar') or .012))
                  for row in stars[target]}
        if len(params) != 1:
            raise ValueError('Conflicting stellar parameters for '+target)
        mass, rs = params.pop()
        for period in args.periods:
            for radius in args.radii:
                summary, trials = run_star_cell(target,visits,baseline,period,radius,mass,rs,args.trials,args.seed)
                totals.append(summary)
                name = f'{target}-P{period:g}-R{radius:g}'
                (args.out/(name+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
                with (args.out/(name+'.csv')).open('w',newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(trials[0]))
                    writer.writeheader()
                    writer.writerows(trials)
    if not args.target:
        (args.out/'limits.json').write_text(json.dumps(occurrence(totals),indent=2)+'\n')


if __name__ == '__main__':
    main()

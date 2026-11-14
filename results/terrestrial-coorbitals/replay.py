"""Propagate the two measured orbits and reproduce the counts in the research note."""
import argparse
import json
from pathlib import Path

import numpy as np
import rebound

ROOT = Path(__file__).resolve().parent
G = 0.0002959122082855911  # Solar GM in au^3/day^2, solar masses
EARTH_MOON = 1 / 328900.5597
MASS = {
    10: 1., 1: 1 / 6023600., 2: 1 / 408523.718655,
    399: EARTH_MOON / (1 + 1 / 81.30057), 4: 1 / 3098703.59,
    5: 1 / 1047.34863, 6: 1 / 3497.9018, 7: 1 / 22902.98,
    8: 1 / 19412.24, 301: EARTH_MOON / (1 + 81.30057),
}
# With insertion order above, particles[3] and [4] are Earth and Mars.
# (target, solution, signed years, step days, random samples, solver, output years)
CASES = {
    'mars-forward': ('2026 RF64', 'jpl', 4000, 1, 64, 'mercurius', 2),
    'mars-backward': ('2026 RF64', 'jpl', -4000, 1, 64, 'mercurius', 2),
    'mars-esa-forward': ('2026 RF64', 'esa', 2500, 1, 24, 'mercurius', 1),
    'mars-esa-backward': ('2026 RF64', 'esa', -2500, 1, 24, 'mercurius', 1),
    'mars-ias15': ('2026 RF64', 'jpl', 2200, 1, 0, 'ias15', 5),
    'massless-mars': ('2026 RF64', 'jpl', 4000, 1, 0, 'ias15', 2),
    'earth-forward': ('2025 XN6', 'jpl', 200, .7, 128, 'ias15', .5),
    'earth-backward': ('2025 XN6', 'jpl', -120, .5, 24, 'mercurius', .5),
    'earth-esa': ('2025 XN6', 'esa', 220, 1, 64, 'mercurius', 1),
}

def sampled_elements(fit, source, count):
    mean = np.array(fit['elements'])
    if source == 'esa':
        cov = np.zeros((6, 6))
        cov[np.triu_indices(6)] = fit['upper_covariance']
        cov += np.triu(cov, 1).T
        seed = 20261115
    else:
        assert fit['order'] == ['e', 'q', 'tp', 'node', 'peri', 'i']
        cov = np.array(fit['covariance'])
        seed = 20261114
    noise = np.random.default_rng(seed).standard_normal((count, 6))
    return np.vstack([mean, mean + noise @ np.linalg.cholesky(cov).T])


def add_test_particle(sim, elements, source, epoch):
    if source == 'jpl':
        e, q, tp, node, peri, inc = elements
        a = q / (1 - e)
        anomaly = np.sqrt(G / a**3) * (epoch - tp)
    else:
        a, e, inc, node, peri, mean_degrees = elements
        anomaly = np.radians(mean_degrees)
    sim.add(primary=sim.particles[0], a=a, e=e, M=anomaly, inc=np.radians(inc),
            Omega=np.radians(node), omega=np.radians(peri))


def integrate(case):
    name, source, duration, dt, count, method, interval = CASES[case]
    data = json.loads((ROOT / 'orbits.json').read_text())
    obj = data['objects'][name]
    start = obj['jpl']['epoch']
    fit = obj[source]
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = method
    sim.dt = dt
    for body, mass in MASS.items():
        x, y, z, vx, vy, vz = data['planet_vectors'][str(start)][str(body)]
        sim.add(m=mass, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.N_active = sim.N
    if source == 'esa':
        # Transform to the same centre-of-mass frame as the JPL initialisation,
        # then carry the massive bodies to the slightly later ESA epoch.
        sim.move_to_com()
        sim.integrate(fit['epoch'] - start)
        sim.t = 0
    for elements in sampled_elements(fit, source, count):
        add_test_particle(sim, elements, source, fit['epoch'])
    sim.move_to_com()
    if case == 'massless-mars':
        sim.particles[4].m = 0
        sim.move_to_com()
    times = np.arange(0., abs(duration) + interval / 4, interval) * np.sign(duration)
    angles = np.empty((len(times), count + 1))
    planet = 4 if name == '2026 RF64' else 3
    for j, year in enumerate(times):
        sim.integrate(year * 365.25)
        sun = sim.particles[0]
        longitude = sim.particles[planet].orbit(primary=sun).l
        for k, particle in enumerate(list(sim.particles)[sim.N_active:]):
            angle = np.degrees(particle.orbit(primary=sun).l - longitude)
            angles[j, k] = (angle + 180) % 360 - 180
    return 2000 + (fit['epoch'] - 2451545) / 365.25 + times, angles


def trailing_samples(years, angles):
    elapsed = np.abs(years - years[0])
    masks = {}
    for horizon in (1000, 2000, 4000):
        if elapsed[-1] >= horizon:
            block = angles[elapsed <= horizon, 1:]
            masks[horizon] = np.all((block > -180) & (block < 0), axis=0)
    return masks


def summarize(case, years, angles):
    nominal = angles[:, 0]
    result = {'nominal_range_deg': np.round([nominal.min(), nominal.max()], 2).tolist()}
    if case == 'massless-mars' or case.endswith('backward') and case.startswith('earth'):
        result['nominal_unwrapped_span_deg'] = round(float(np.ptp(
            np.degrees(np.unwrap(np.radians(nominal))))), 2)
    elif case.startswith('mars'):
        if angles.shape[1] > 1:
            for horizon, mask in trailing_samples(years, angles).items():
                result[f'samples_trailing_{horizon}yr'] = int(mask.sum())
    else:
        for year in (2065, 2150, 2200):
            block = angles[years <= year, 1:]
            confined = np.all(np.abs(block) < 90, axis=0)
            crossings = np.sum(np.diff(block > 0, axis=0), axis=0)
            result[f'samples_within_90deg_through_{year}'] = int(confined.sum())
            result[f'samples_with_two_crossings_through_{year}'] = int(
                (confined & (crossings >= 2)).sum())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=[*CASES, 'all'])
    parser.add_argument('--save', type=Path, help='Directory for trajectory arrays')
    args = parser.parse_args()
    jpl_mars_masks = {}
    for label in CASES if args.case == 'all' else [args.case]:
        t, phi = integrate(label)
        if args.save:
            args.save.mkdir(exist_ok=True, parents=True)
            np.savez_compressed(args.save / f'{label}.npz', time=t, angle=phi)
        print(label, json.dumps(summarize(label, t, phi)), flush=True)
        if args.case == 'all' and label in ('mars-forward', 'mars-backward'):
            jpl_mars_masks[label] = trailing_samples(t, phi)
    if jpl_mars_masks:
        # The JPL seed and covariance are identical in the two directions.
        forward = jpl_mars_masks['mars-forward']
        backward = jpl_mars_masks['mars-backward']
        joint = {f'samples_trailing_plus_minus_{h}yr': int((mask & backward[h]).sum())
                 for h, mask in forward.items()}
        print('mars-jpl-both-directions', json.dumps(joint))

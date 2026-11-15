"""Fit linked and separate astrometry with Project Pluto's pinned forward model."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
SCALE = np.array([1.e-4] * 3 + [1.e-6] * 3)
SOLAR_GM = 0.0002959122082855911  # AU^3/day^2, DE440 nominal solar GM
PAIRS = {
    'VE138': (('2007 VE138', '2026 VV15'), {'595', 'J95'}),
    'RB15': (('2017 RB15', '2026 RG39'), {'F51', 'I52'}),
    'SA12': (('2026 SA12', '2026 SW11'), {'H36'}),
    'VO': (('2026 VO', '2026 VR1'), {'291'}),
}


def seed_state(elements):
    a, e, m, inc, node, peri, n = (float(elements[key]) for key in
                                   ('a', 'e', 'ma', 'i', 'om', 'w', 'n'))
    m = np.radians(m)
    eccentric_anomaly = m
    for _ in range(15):
        eccentric_anomaly -= ((eccentric_anomaly - e * np.sin(eccentric_anomaly) - m)
                              / (1 - e * np.cos(eccentric_anomaly)))
    c, s = np.cos(eccentric_anomaly), np.sin(eccentric_anomaly)
    b = np.sqrt(1 - e * e)
    rotate = Rotation.from_euler('ZXZ', [node, inc, peri], degrees=True)
    return np.r_[rotate.apply(a * np.array([c-e, b*s, 0])),
                 rotate.apply(np.radians(n) * a / (1-e*c) * np.array([-s, b*c, 0]))]


class Forward:
    """Persistent Find_Orb process. Results are ordered by the observations' UTC times."""

    def __init__(self, build, work, observations, designation, epoch, seed):
        self.obs = observations
        self.work = work
        source = build / 'find_orb'
        for path in source.iterdir():
            if path.name != 'ObsCodes.htm' and path.suffix.lower() in {
                    '.txt', '.htm', '.html', '.sof', '.hdr', '.def', '.json', '.pgm'}:
                (work / path.name).symlink_to(path)
        # Find_Orb's short-form station list takes precedence. S11 postdates the
        # pinned software release. Supply the exact MPC constants used for all sites.
        shutil.copyfile(HERE / 'stations.txt', work / 'ObsCodes.htm')
        shutil.copyfile(source / 'environ.def', work / 'environ.dat')
        obsfile = work / 'astrometry.psv'
        with obsfile.open('w') as handle:
            handle.write('# version=2022\nprovID|mode|stn|obsTime|ra|dec|astCat\n')
            for row in observations:
                handle.write('|'.join((designation, row['mode'], row['site'], row['utc'],
                                       row['ra_deg'], row['dec_deg'], row['catalog'])) + '\n')
        self.log = (work / 'engine.log').open('w+')
        self.proc = subprocess.Popen(
            [str(source / 'fo'), str(obsfile), '-x' + str(work) + '/',
             '-v' + epoch + ',' + ','.join(map(str, seed)),
             'LINUX_JPL_FILENAME=' + str(build / 'de440.bin'), 'PERTURBERS=7fe', '-i'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log,
            text=True, bufsize=1, cwd=work, env={**os.environ, 'FO_FORWARD_MODEL': '1'})
        self.read()

    def read(self):
        for line in self.proc.stdout:
            if 'RES ' in line:
                values = np.fromstring(line.split('RES ')[-1], sep=' ')
                if len(values) != 2 * len(self.obs) or not np.isfinite(values).all():
                    raise ValueError('Invalid Find_Orb residual count or values')
                return values.reshape(-1, 2)
        self.log.seek(0)
        raise RuntimeError('Find_Orb exited: ' + self.log.read())

    def __call__(self, state):
        self.proc.stdin.write(' '.join(map(str, state)) + '\n')
        self.proc.stdin.flush()
        return self.read()

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()
        self.log.close()


def derivative(model, state, delta):
    return np.stack([(model(state + step*SCALE*delta) - model(state - step*SCALE*delta))
                     / (2*delta) for step in np.eye(6)], axis=2)


def fit(model, seed, mask):
    def residual(offset):
        return model(seed + offset*SCALE)[mask].ravel()

    def jacobian(offset):
        return derivative(model, seed + offset*SCALE, 0.002)[mask].reshape(-1, 6)

    solution = least_squares(residual, np.zeros(6), jac=jacobian, x_scale='jac',
                             ftol=1.e-10, xtol=1.e-10, gtol=1.e-6, max_nfev=140)
    if not solution.success:
        raise ValueError(solution.message)
    state = seed + solution.x * SCALE
    return state, model(state)


def summarize(residual):
    distance = np.linalg.norm(residual, axis=1)
    return {'n': len(distance), 'radial_rms_arcsec': float(np.mean(distance**2)**0.5),
            'radial_median_arcsec': float(np.median(distance)),
            'radial_max_arcsec': float(np.max(distance))}


def elements(state):
    r, v = state[:3], state[3:]
    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h)
    ev = np.cross(v, h) / SOLAR_GM - r / np.linalg.norm(r)
    e = np.linalg.norm(ev)
    a = 1 / (2 / np.linalg.norm(r) - v @ v / SOLAR_GM)
    hhat, nhat, ehat = h / np.linalg.norm(h), n / np.linalg.norm(n), ev / e
    peri = np.arctan2(np.cross(nhat, ehat) @ hhat, nhat @ ehat)
    anomaly = np.arctan2(np.cross(ehat, r) @ hhat, ehat @ r)
    eccentric = 2 * np.arctan2(np.sqrt(1-e)*np.sin(anomaly/2),
                               np.sqrt(1+e)*np.cos(anomaly/2))
    return dict(a_au=float(a), e=float(e), i_deg=float(np.degrees(np.arccos(hhat[2]))),
                node_deg=float(np.degrees(np.arctan2(n[1], n[0])) % 360),
                peri_deg=float(np.degrees(peri) % 360),
                mean_anomaly_deg=float(np.degrees(eccentric-e*np.sin(eccentric)) % 360),
                period_days=float(2*np.pi*np.sqrt(a**3 / SOLAR_GM)))


def covariance(model, state, mask, residuals):
    jac = derivative(model, state, 0.001)[mask]
    r = residuals[mask]
    inv = np.linalg.pinv(jac.reshape(-1, 6))
    h_inv = inv @ inv.T
    formal = h_inv * (r*r).sum() / (r.size - 6)
    keys = np.array([(o['site'] + '/' + o['utc'][:10]) for o in model.obs])[mask]
    groups = np.unique(keys)
    meat = np.zeros((6, 6))
    for group in groups:
        selection = keys == group
        score = jac[selection].reshape(-1, 6).T @ r[selection].ravel()
        meat += np.outer(score, score)
    correction = len(groups) / (len(groups) - 1) * (r.size - 1) / (r.size - 6)
    clustered = h_inv @ meat @ h_inv * correction
    outer_scale = np.outer(SCALE, SCALE)
    return {'formal': (formal*outer_scale).tolist(),
            'site_utc_date_clustered': (clustered*outer_scale).tolist(),
            'site_utc_date_groups': len(groups)}


def run_pair(build, key, modes, outdir, rows, seeds):
    names, withheld = PAIRS[key]
    observations = sorted((o for o in rows if o['designation'] in names), key=lambda o:o['utc'])
    labels = np.array([o['designation'] for o in observations])
    sites = np.array([o['site'] for o in observations])
    epoch = seeds[names[1]]['epoch']
    if seeds[names[0]]['epoch'] != epoch:
        raise ValueError('This replay expects catalog seeds at a common epoch')
    answer = {'designations': names, 'epoch_jd_tdb': float(epoch),
              'frame': 'heliocentric J2000 ecliptic and equinox',
              'state_units': ['au']*3 + ['au/day']*3, 'holdout_sites': sorted(withheld),
              'residual_sign': 'computed - observed; RA residual includes cos(dec)'}
    with tempfile.TemporaryDirectory() as tmp:
        engine = Forward(build, Path(tmp), observations, names[1], epoch, seed_state(seeds[names[1]]))
        try:
            for mode in modes:
                mask = np.ones(len(labels), dtype=bool)
                source = names[1]
                if mode == 'hold':
                    mask = ~np.isin(sites, list(withheld))
                if mode in ('old', 'new'):
                    source = names[int(mode == 'new')]
                    mask = labels == source
                state, residuals = fit(engine, seed_state(seeds[source]), mask)
                result = {'state': state.tolist(), 'elements': elements(state),
                          'training': summarize(residuals[mask])}
                if not mask.all():
                    result['withheld'] = summarize(residuals[~mask])
                if mode in ('full', 'new'):
                    result['covariance'] = covariance(engine, state, mask, residuals)
                answer[mode] = result
                with (outdir / (key + '-' + mode + '.csv')).open('w') as handle:
                    writer = csv.writer(handle, lineterminator='\n')
                    writer.writerow(['obsid','designation','utc','site','in_fit','ra_cosdec_arcsec','dec_arcsec'])
                    for o, used, (dx, dy) in zip(observations, mask, residuals):
                        writer.writerow([o['obsid'],o['designation'],o['utc'],o['site'],int(used),f'{dx:.6f}',f'{dy:.6f}'])
        finally:
            engine.close()
    (outdir / (key + '.json')).write_text(json.dumps(answer, indent=2) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('build', type=Path, help='Find_Orb build directory from build.sh')
    ap.add_argument('output', type=Path, help='Output directory')
    ap.add_argument('--pair', choices=list(PAIRS), action='append')
    ap.add_argument('--modes', nargs='+', choices=['full','hold','old','new'],
                    default=['full','hold','old','new'])
    args = ap.parse_args()
    build = args.build.resolve()
    if hashlib.sha256((build / 'de440.bin').read_bytes()).hexdigest() != '29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062':
        raise ValueError('DE440 checksum mismatch')
    with (HERE / 'astrometry.csv').open() as f:
        rows = list(csv.DictReader(f))
    seeds = {r['pdes']: r for r in json.loads((HERE / 'seeds.json').read_text())}
    args.output.mkdir(exist_ok=True, parents=True)
    for key in args.pair or PAIRS:
        run_pair(build, key, args.modes, args.output, rows, seeds)


if __name__ == '__main__':
    main()

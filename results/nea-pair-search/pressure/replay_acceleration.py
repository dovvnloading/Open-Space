"""Test inverse-square radial and transverse acceleration with the pinned orbit engine."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'orbits'))
import replay as orbit

SCALE = np.r_[orbit.SCALE, 1.e-10, 1.e-12]
STEP = np.r_[[0.0005]*6, 0.03, 0.03]


def jacobian(model, state, axes, step_multiplier=1.):
    steps = STEP * step_multiplier
    directions = np.diag(steps * SCALE)
    return np.stack([(model(state + directions[i]) - model(state - directions[i]))
                     / (2 * steps[i]) for i in axes], axis=2)


def fit(model, state, axes, mask, weights):
    x = state.copy()
    for iteration in range(30):
        r = model(x)
        j = jacobian(model, x, axes)
        a = (j * weights[:, None, None])[mask].reshape(-1, len(axes))
        target = (r * weights[:, None])[mask].ravel()
        delta = np.linalg.lstsq(a, -target, rcond=None)[0]
        current = np.linalg.norm(target)
        for power in range(15):
            step = 2.**-power
            trial = x.copy()
            trial[axes] += SCALE[axes] * delta * step
            if np.linalg.norm((model(trial)*weights[:, None])[mask]) <= current + 1.e-7:
                break
        x = trial
        if np.linalg.norm(delta*step) < 1.e-4:
            break
    return x, model(x)


def covariance(jac, residual, groups, scale):
    a = jac.reshape(-1, jac.shape[-1])
    y = residual.ravel()
    inv = np.linalg.pinv(a)
    h = inv @ inv.T
    formal = h * (y @ y) / (len(y) - a.shape[1])
    scores = np.array([jac[groups == k].reshape(-1, len(scale)).T @ residual[groups == k].ravel()
                       for k in np.unique(groups)])
    count = len(scores)
    correction = count / (count - 1) * (len(y)-1)/(len(y)-a.shape[1])
    clustered = h @ scores.T @ scores @ h * correction
    return {'formal': (formal*np.outer(scale, scale)).tolist(),
            'site_date': (clustered*np.outer(scale, scale)).tolist()}


def diagnostic_offsets(jac, residuals, rows, state):
    site = np.array([o['site'] for o in rows])
    date = np.array([o['utc'][:10] for o in rows])
    catalog = np.array([o['catalog'] for o in rows])
    old = date < '2026'
    all_rows = np.ones(len(rows), bool)
    def solve(mask, offsets=None, prior=None):
        j = jac[mask].reshape(-1, 8)
        y = residuals[mask].ravel()
        if offsets is not None:
            nuisance = offsets[np.repeat(mask,2)]
            j = np.column_stack((j, nuisance))
            j = np.vstack((j, np.c_[np.zeros((nuisance.shape[1],8)), np.eye(nuisance.shape[1])/prior]))
            y = np.r_[y, np.zeros(nuisance.shape[1])]
        dx = np.linalg.lstsq(j, -y, rcond=None)[0][:8]
        return (state + dx*SCALE)[6:].tolist()
    out = {label:{v:solve(values != v) for v in np.unique(values)}
           for label, values in [('drop_site', site), ('drop_utc_date', date)]}
    for label, groups in [('old_global', np.array(['old']*len(rows))), ('old_catalog', catalog)]:
        offsets = np.column_stack([np.outer(old & (groups==v), np.eye(2)[i]).ravel()
                                   for v in np.unique(groups[old]) for i in (0,1)])
        out[label] = {str(prior):solve(all_rows, offsets, prior) for prior in (0.2, 0.5, 1.0)}
    return out


def run(build, output, key):
    labels, excluded = orbit.PAIRS[key]
    with (orbit.HERE/'astrometry.csv').open() as f:
        rows = sorted([r for r in csv.DictReader(f) if r['designation'] in labels],
                      key=lambda r:r['utc'])
    nominal = json.loads((orbit.HERE/'fits'/(key+'.json')).read_text())
    start = np.r_[nominal['full']['state'], 0., 0.]
    date = np.array([r['utc'][:10] for r in rows])
    site = np.array([r['site'] for r in rows])
    all_rows = np.ones(len(rows), bool)
    unit = np.ones(len(rows))
    groups = np.char.add(site, date)
    old = np.array([r['designation'] == labels[0] for r in rows])
    axes = {'grav':list(range(6)), 'A1':list(range(7)), 'A2':[*range(6),7], 'A12':list(range(8))}
    cases = [('all', all_rows, unit, m) for m in ['grav','A1','A2','A12']]
    cases += [('sites', ~np.isin(site, list(excluded)), unit, 'A12')]
    if key == 'VE138':
        reported = np.array([max(float(o['rms_ra_arcsec'] or 0.8),
                                 float(o['rms_dec_arcsec'] or 0.8), 0.5) for o in rows])
        balance = 1/np.sqrt([sum(groups==g) for g in groups])
        cases += [('equal_site_date', all_rows, balance, 'A12'),
                  ('reported_floor', all_rows, 1/reported, 'A12'),
                  ('omit_old_703', ~(old & (site=='703')), unit, 'A12')]
        for label, mask in [('first_new2',old|(date=='2026-11-10')),
                            ('first_new8',old|(date<='2026-11-12')),
                            ('omit_last_new',date!='2026-11-14'),('omit_last_old',date!='2007-11-13')]:
            for m in ['grav','A12']:
                cases.append((label,mask,unit,m))
    out = {'designations':labels,'epoch_jd_tdb':nominal['epoch_jd_tdb'],
           'units':'au, au/day, au/day^2 at 1 au (radial outward and prograde transverse)'}
    with tempfile.TemporaryDirectory() as t:
        engine = orbit.Forward(build, Path(t), rows, labels[1], str(nominal['epoch_jd_tdb']), start[:6])
        try:
            for label, mask, weights, mode in cases:
                active = axes[mode]
                seed = start if label=='all' or mode=='grav' else np.array(out['all/A12']['state'])
                state, residuals = fit(engine, seed, active, mask, weights)
                result = {'state':state.tolist(), 'training':orbit.summarize(residuals[mask])}
                if not mask.all():
                    result['withheld'] = orbit.summarize(residuals[~mask])
                if label == 'first_new2' and mode == 'A12':
                    result['fit_status'] = 'iteration limit; poorly constrained; exclude from force evidence'
                j = jacobian(engine, state, active)
                result['covariance'] = covariance((j*weights[:,None,None])[mask],
                                                  (residuals*weights[:,None])[mask], groups[mask], SCALE[active])
                if mode=='A12' and label=='all':
                    result['offset_and_deletion_tests_linearized'] = diagnostic_offsets(j, residuals, rows, state)
                    result['finite_difference_linear_refits'] = {}
                    for step_factor in [0.5,2.]:
                        jac = jacobian(engine,state,active,step_factor).reshape(-1,8)
                        delta = np.linalg.lstsq(jac,-residuals.ravel(),rcond=None)[0]
                        result['finite_difference_linear_refits'][str(step_factor)] = (state+delta*SCALE)[6:].tolist()
                    with (output/(key+'-residuals.csv')).open('w') as h:
                        w=csv.writer(h,lineterminator='\n')
                        w.writerow(['obsid','utc','site','ra_cosdec_arcsec','dec_arcsec'])
                        for row, (east,north) in zip(rows,residuals):
                            w.writerow([row['obsid'],row['utc'],row['site'],f'{east:.6f}',f'{north:.6f}'])
                out[label+'/'+mode] = result
                (output/(key+'.json')).write_text(json.dumps(out,indent=2)+'\n')
        finally:
            engine.close()


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('build',type=Path)
    cli.add_argument('output',type=Path)
    cli.add_argument('--pair',choices=['VE138','RB15'],action='append')
    args = cli.parse_args()
    os.environ['FO_ACCEL']='1'
    args.output.mkdir(parents=True,exist_ok=True)
    for key in args.pair or ['VE138','RB15']:
        run(args.build.resolve(),args.output,key)

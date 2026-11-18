"""Retest the three remaining single-band candidates with near-FS binary starts."""
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from fetch import OUT, write_csv
from flag_refits import SEED, WORK, inputs
from peak_models import solve

NAMES = {'PRIME-2024-BLG-245', 'PRIME-2025-BLG-0239', 'PRIME-2025-BLG-0360'}


def compare(job):
    event = job[0]
    calibrated = not event.startswith('PRIME-2024')
    seed = json.loads((SEED / (event + '.json')).read_text())
    ckpt = Path('staging/flag-calibrated') / (event + '-1.15.json') if calibrated else WORK / (event + '.json')
    r = json.loads(ckpt.read_text())
    seed.update(ps_t0=r['ps_t0'], ps_te=r['ps_te_or_rise'], ps_u0=r['ps_u0_or_fall'],
                fs_t0=r['fs_t0'], fs_te=r['fs_te_or_rise'], fs_u0=r['fs_u0_or_fall'], fs_rho=r['rho'])
    if calibrated:
        lookup = list(csv.DictReader((OUT / 'calibrated_alternatives.csv').open()))
        c = next(v for v in lookup if v['event'] == event)
        hw, dt = float(c['window_halfwidth']), float(c['bin_d'])
        cap = [1.15 * r['baseline_counts']]
    else:
        lookup = list(csv.DictReader((OUT / 'flag_alternatives.csv').open()))
        c = next(v for v in lookup if v['event'] == event)
        hw, dt, cap = float(c['width_d']), float(c['bin_d']), None
    first = solve(job, seed, 1.15, 32, False, hw, dt, source_flux_caps=cap)
    second = solve(job, seed, 1.15, 64, False, hw, dt, previous=first['models'], source_flux_caps=cap)
    models = {**first['models'], **second['models']}
    return dict(event=event, window_halfwidth=hw, bin_d=dt, n_peak_bins=first['n_peak_window'],
                cap_counts=cap[0] if cap else '', fs_chi2=models['fs']['chi2'], fs_te=models['fs']['te_or_rise_d'],
                fs_rho=models['fs']['rho'], binary_chi2=models['bin']['chi2'], binary_q=models['bin']['q'],
                binary_s=models['bin']['s'], binary_rho=models['bin']['rho'], binary_converged=models['bin']['converged'],
                exp_chi2=models['exp']['chi2'], fs_delta_against_binary=round(models['fs']['chi2']-models['bin']['chi2'],3))


if __name__ == '__main__':
    jobs = [j for j in inputs(Path.cwd()) if j[0] in NAMES]
    with ProcessPoolExecutor(3) as pool:
        write_csv(OUT / 'prime_shape.csv', list(pool.map(compare, jobs)))

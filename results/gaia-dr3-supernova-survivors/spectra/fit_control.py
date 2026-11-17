"""Measure the known LP 40-365 star J1109+0001 in DESI DR1.

Run: python fit_control.py [--cache DIRECTORY] [--output DIRECTORY]
Dependencies: Python 3.10+, numpy, scipy. All wavelengths are Angstrom.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy.optimize import least_squares

C = 299792.458
STAR = 'WDJ110907.96+000132.92'
BASE = f'https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/J1109/{STAR}/'
URL = BASE + STAR + '_DESI_DR1_spectra.tar.gz'
EXPECTED_SHA256 = 'bf4bc44dc10b37ef155062be01e530f33fd38d2d615e4fe10d04fcb4bd9c0d70'
# Four isolated, well-spaced Mg I lines with no line-ratio assumption.
MG = [4167.271, 4351.906, 4702.991, 5528.405]
# Ca H can blend with H epsilon and both Ca H and K with interstellar Ca.
CA = [3933.6614, 3968.4673]
LINES = {'Mg I': MG, 'Ca II': CA, 'H alpha': [6562.80], 'H beta': [4861.33],
         'He I': [5875.62], 'Si II': [6347.11, 6371.37], 'Ne I': [6402.25],
         'C II': [4267.18], 'Fe II': [5018.44], 'O I triplet': [7771.94, 7774.17, 7775.39]}


def vacuum(air):
    # Standard dry-air conversion (Morton 2000), optical wavelengths.
    s2 = (1e4 / np.asarray(air)) ** 2
    return air * (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


def factor(v):
    return np.sqrt((1 + v / C) / (1 - v / C))


def data_window(spectra, lo, hi):
    arrays = [a[(a[:, 0] >= lo) & (a[:, 0] <= hi)] for a in spectra]
    d = np.concatenate(arrays)
    d = d[(d[:, 4] == 0) & (d[:, 3] > 0) & np.isfinite(d).all(axis=1)]
    if len(d) < 12:
        raise ValueError(f'Insufficient valid pixels between {lo} and {hi}')
    return d[:, 0], d[:, 1], np.sqrt(d[:, 3])


def velocity_fit(spectra, air, initials=(-200, 0, 200)):
    windows = []
    guess = []
    low = []
    high = []
    for rest in vacuum(np.asarray(air)):
        w, flux, weights = data_window(spectra, rest - 13, rest + 13)
        continuum = np.median(flux)
        windows.append((rest, w, flux, weights))
        guess.extend([continuum, 0., max(0.1, continuum * 0.3), 1.5])
        low.extend([-np.inf, -np.inf, 0., 0.45])
        high.extend([np.inf, np.inf, np.inf, 5.])

    def residual(p):
        res = []
        doppler = factor(p[0])
        for i, (rest, w, flux, weights) in enumerate(windows):
            cont, slope, depth, sigma = p[1 + 4*i:5 + 4*i]
            line = np.exp(-0.5 * ((w - rest * doppler) / sigma) ** 2)
            res.append((cont + slope * (w - rest) - depth * line - flux) * weights)
        return np.concatenate(res)

    fits = [least_squares(residual, [v0] + guess, bounds=([-350] + low, [350] + high),
                          x_scale='jac') for v0 in initials]
    fit = min(fits, key=lambda r: r.cost)
    dof = fit.fun.size - fit.x.size
    cov = np.linalg.pinv(fit.jac.T @ fit.jac)
    err = np.sqrt(cov[0, 0] * max(1, 2 * fit.cost / dof))
    return dict(rv=float(fit.x[0]), formal_sigma=float(err), chi2=float(2*fit.cost), dof=dof)


def ew_window(spectra, wavelengths, v, pad=7):
    # Direct flux-sum pseudo-EW. Exclude a 3-A gap on each side of
    # integration window, then fit a linear continuum to the next 12 A.
    centers = vacuum(np.asarray(wavelengths)) * factor(v)
    lo, hi = min(centers) - pad, max(centers) + pad
    w, flux, weights = data_window(spectra, lo - 15, hi + 15)
    if len(set(w)) != len(w):
        raise ValueError('Multiple arms in EW interval')
    x = np.column_stack((np.ones_like(w), w - np.mean(centers)))
    line = (w >= lo) & (w <= hi)
    side = (w < lo - 3) | (w > hi + 3)
    dx = 0.8
    if np.sum(line) * dx < 0.90 * (hi - lo):
        return None
    design = x[side] * weights[side, None]
    pars = np.linalg.lstsq(design, flux[side] * weights[side], rcond=None)[0]
    cov = np.linalg.inv(design.T @ design)
    redchi = np.sum(((x[side] @ pars - flux[side]) * weights[side])**2)/(side.sum()-2)
    cont = x[line] @ pars
    if np.min(cont) <= 0:
        return None
    ew = dx * np.sum(1.0 - flux[line] / cont)
    grad = dx * np.sum((flux[line]/cont**2)[:, None] * x[line], axis=0)
    var = dx**2 * np.sum((1/(weights[line]*cont))**2) + grad @ cov @ grad
    return float(ew), float(np.sqrt(var * max(1, redchi)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('.'))
    parser.add_argument('--output', type=Path, default=Path('.'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.cache / 'J1109.tar.gz'
    if not path.exists():
        with urlopen(URL, timeout=90) as response:
            data = response.read()
    else:
        data = path.read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    if sha != EXPECTED_SHA256:
        raise ValueError('Archive fingerprint has changed. Inspect the new data before fitting.')
    if not path.exists():
        path.write_bytes(data)
    spectra = {}
    manifest = []
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        for member in tar.getmembers():
            if member.isfile() and member.name.endswith('.dat'):
                raw = tar.extractfile(member).read()
                key = member.name.split('_DESI_DR1_')[1].rsplit('-', 1)[0]
                spectrum = np.loadtxt(io.BytesIO(raw))
                spectra.setdefault(key, []).append(spectrum)
                manifest.append(dict(name=member.name, sha256=hashlib.sha256(raw).hexdigest(),
                                     n_pixels=len(spectrum), bytes=len(raw)))
    receipt = dict(source_id='3804182280735442560', origin=URL, archive_sha256=sha,
                   archive_bytes=len(data), members=manifest)
    (args.output/'spectrum_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    rows = []
    for key, arms in sorted(spectra.items()):
        for label, line_set in [('Mg I', MG), ('Ca II', CA)]:
            rows.append(dict(spectrum=key, line_set=label, **velocity_fit(arms, line_set)))
    for line in MG:
        rows.append(dict(spectrum='coadd', line_set=f'Mg I {line}', **velocity_fit(spectra['coadd'], [line])))
    with (args.output/'rvs.csv').open('w') as fp:
        out = csv.DictWriter(fp, fieldnames=rows[0].keys())
        out.writeheader()
        out.writerows(rows)
    rv = next(r['rv'] for r in rows if r['spectrum']=='coadd' and r['line_set']=='Mg I')
    with (args.output/'features.csv').open('w') as fp:
        out=csv.writer(fp)
        out.writerow(['feature','air_A','pseudo_EW_A','formal_sigma_A','fixed_RV_kms','half_padding_A'])
        for label, waves in LINES.items():
            groups = [waves] if label == 'O I triplet' else [[w] for w in waves]
            for group in groups:
                val = ew_window(spectra['coadd'],group,rv)
                out.writerow([label,'/'.join(map(str,group)), *(val or ('','')),rv,7])


if __name__ == '__main__':
    main()

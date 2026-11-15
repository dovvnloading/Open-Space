"""Check broad emission wings, continuum placement, and unblended hydrogen lines."""
import gzip

import numpy as np

from analyze import C, CENTERS, ROOT, write_csv


def line_fit(s, centers, outer, inner, exclude_pas):
    w = s['wave']
    velocity = (w[:, None] / centers - 1) * C
    vabs = abs(velocity)
    good = (w > 8395) & (w < 8830) & (s['ivar'] > 0) & (s['andmask'] == 0)
    if exclude_pas:
        hydrogen = np.array([8440.3, 8469.6, 8600.75, 8752.89])  # vacuum
        good &= np.all(abs(w[:, None] - hydrogen) > 15, axis=1)
    if inner:
        good &= np.all(vabs >= inner, axis=1)
    t = (w[good] - 8600) / 220
    basis = np.column_stack((np.ones_like(t), t, t**2, (vabs[good] < outer)))
    if np.linalg.matrix_rank(basis) != basis.shape[1]:
        return None
    precision = (basis.T * s['ivar'][good]) @ basis
    cov = np.linalg.inv(precision)
    coeff = cov @ ((basis.T * s['ivar'][good]) @ s['flux'][good])
    count = len(centers)
    amp = coeff[-count:]
    shape_precision = np.linalg.inv(cov[-count:, -count:])
    return {'combined_z': np.sum(shape_precision @ amp) / np.sqrt(shape_precision.sum()),
            **{f'z_{l:.0f}': v for l, v in zip(centers, amp / np.sqrt(cov.diagonal()[-count:]))}}


def local_hydrogen(s, center):
    w, iv = s['wave'], s['ivar']
    d = w - center
    good = (iv > 0) & (s['andmask'] == 0) & (abs(d) < 37)
    lines = (abs(d[good]) < 15)
    X = np.column_stack((np.ones(good.sum()), d[good], lines))
    cov = np.linalg.inv((X.T * iv[good]) @ X)
    coef = cov @ ((X.T * iv[good]) @ s['flux'][good])
    return coef[-1] / np.sqrt(cov[-1, -1])


def main():
    with gzip.open(ROOT / 'spectra.csv.gz', 'rt') as f:
        spectra = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf8')
    rows = []
    for epoch in sorted(set(spectra['epoch'])):
        s = spectra[spectra['epoch'] == epoch]
        for inner in (0, 220):
            for excise in (False, True):
                result = line_fit(s, CENTERS, 650, inner, excise)
                if result is not None:
                    rows.append(dict(epoch=epoch, inner_kms=inner, excise_unblended_paschen=excise,
                                     p14_z=local_hydrogen(s, 8600.75), p12_z=local_hydrogen(s, 8752.89),
                                     **result))
    write_csv(ROOT / 'profile_checks.csv', rows)


if __name__ == '__main__':
    main()

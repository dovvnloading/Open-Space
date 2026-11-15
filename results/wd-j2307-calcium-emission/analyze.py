"""Fit fixed-wavelength Ca II triplets in four public spectra. Run with NumPy."""
import csv
import gzip
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CENTERS = np.array([8500.36, 8544.44, 8664.52])  # vacuum angstroms
C = 299792.458
BAND = (8425, 8735)
HALF_WIDTH = 650.0  # km/s


def fit(wave, flux, ivar, mask, centers=CENTERS, width=HALF_WIDTH, degree=2):
    """Fit polynomial continuum and independent mean line excesses simultaneously."""
    use = ((wave > BAND[0]) & (wave < BAND[1]) & (ivar > 0) & (mask == 0)
           & np.isfinite(flux))
    w, y, iv = wave[use], flux[use], ivar[use]
    coordinate = (w - 8580) / 160
    continuum = np.vander(coordinate, degree + 1, increasing=True)
    lines = np.abs((w[:, None] / centers - 1) * C) < width
    design = np.column_stack((continuum, lines))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        return None
    cov = np.linalg.inv((design.T * iv) @ design)
    coeff = cov @ ((design.T * iv) @ y)
    amplitudes = coeff[-3:]
    line_cov = cov[-3:, -3:]
    precision = np.linalg.inv(line_cov)
    combined = np.sum(precision @ amplitudes) / np.sqrt(precision.sum())
    z = amplitudes / np.sqrt(line_cov.diagonal())
    center_basis = np.vander((centers - 8580) / 160, degree + 1, increasing=True)
    baseline = center_basis @ coeff[:-3]
    widths = 2 * centers * width / C
    ews = widths * amplitudes / baseline
    # Delta-method derivative of total emission equivalent width, including continuum.
    jac = np.r_[-(widths * amplitudes / baseline**2) @ center_basis, widths / baseline]
    ew_error = np.sqrt(jac @ cov @ jac)
    return {
        'combined_z': combined, 'min_line_z': min(z), 'z_8500': z[0],
        'z_8544': z[1], 'z_8665': z[2], 'ew_sum_A': ews.sum(),
        'ew_sum_error_A': ew_error,
        'n8500': int(lines[:, 0].sum()), 'n8544': int(lines[:, 1].sum()),
        'n8665': int(lines[:, 2].sum())
    }


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        out = csv.DictWriter(f, fieldnames=rows[0])
        out.writeheader()
        out.writerows(rows)


def main():
    with gzip.open(ROOT / 'spectra.csv.gz', 'rt') as f:
        spectra = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf8')
    rows = []
    for epoch in sorted(set(spectra['epoch'])):
        s = spectra[spectra['epoch'] == epoch]
        for mask_name, mask in [('and', s['andmask']), ('or', s['ormask']),
                                ('positive_ivar', np.zeros(len(s)))]:
            for degree in (1, 2, 3):
                for width in (500, 650, 800):
                    result = fit(s['wave'], s['flux'], s['ivar'], mask, degree=degree, width=width)
                    if result is not None:
                        rows.append(dict(epoch=epoch, mask=mask_name, degree=degree,
                                         half_width_kms=width, **result))
    write_csv(ROOT / 'epoch_fits.csv', rows)
    for row in rows:
        if (row['mask'], row['degree'], row['half_width_kms']) == ('and', 2, 650):
            print(f"{row['epoch']:12} z={row['combined_z']:.2f} "
                  f"EW={row['ew_sum_A']:.2f} +/- {row['ew_sum_error_A']:.2f} A")


if __name__ == '__main__':
    main()

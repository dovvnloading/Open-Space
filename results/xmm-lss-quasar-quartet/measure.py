"""Replay geometry and line fits from the archived DESI rows and spectral pixels."""
import argparse
import csv
import gzip
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
LINES = [('C IV', 1549.48, 'B'), ('C III]', 1908.73, 'B'), ('Mg II', 2798.75, 'Z')]
REFERENCE_Z = 1.9295
C = 299792.458


def read_csv(name):
    with (HERE / name).open() as stream:
        return list(csv.DictReader(stream))


def extract(coadd):
    from astropy.io import fits
    with fits.open(coadd) as hdus, gzip.open(HERE / 'line-pixels.csv.gz', 'wt') as out:
        writer = csv.writer(out)
        writer.writerow(['member', 'line', 'wavelength', 'flux', 'ivar'])
        for member in read_csv('members.csv'):
            row, = np.flatnonzero(hdus['FIBERMAP'].data['TARGETID'] == int(member['targetid']))
            for line, rest_wave, arm in LINES:
                wave = hdus[arm + '_WAVELENGTH'].data
                flux = hdus[arm + '_FLUX'].data[row]
                ivar = hdus[arm + '_IVAR'].data[row]
                mask = hdus[arm + '_MASK'].data[row]
                keep = (abs(wave / (1 + REFERENCE_Z) - rest_wave) < 110) & (ivar > 0) & (mask == 0)
                writer.writerows((member['member'], line, w, f, v) for w, f, v in zip(wave[keep], flux[keep], ivar[keep]))


def fit_line(wave, flux, ivar, expected):
    def model(p, x):
        baseline, slope, amplitude, center, sigma = p
        return baseline + slope * (x - expected) / 200 + amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)
    amplitude = max(.01, np.quantile(flux, .93) - np.quantile(flux, .2))
    fit = least_squares(lambda p: (model(p, wave) - flux) * np.sqrt(ivar),
                        [np.median(flux), 0, amplitude, expected, 30],
                        bounds=([-np.inf, -np.inf, 0, expected - 70, 5],
                                [np.inf, np.inf, np.inf, expected + 70, 130]), loss='soft_l1')
    return fit.x, model


def main():
    args = argparse.ArgumentParser()
    args.add_argument('--extract', metavar='COADD_FITS', help='regenerate pixels from DR1 main/dark healpix 17689')
    args.add_argument('--plot', action='store_true')
    options = args.parse_args()
    if options.extract:
        extract(options.extract)

    members = read_csv('members.csv')
    ra, dec, z = [np.array([float(m[k]) for m in members]) for k in ['ra_deg', 'dec_deg', 'z']]
    r, d = np.deg2rad([ra, dec])
    xyz = np.stack((np.cos(r) * np.cos(d), np.sin(r) * np.cos(d), np.sin(d)), axis=1)
    zmean = z.mean()
    da = C / 70 / (1 + zmean) * quad(lambda x: (.3 * (1 + x) ** 3 + .7) ** -.5, 0, zmean)[0]
    for i, j in itertools.combinations(range(4), 2):
        theta = 2 * np.arcsin(np.linalg.norm(xyz[i] - xyz[j]) / 2)
        print(members[i]['member'] + members[j]['member'], f'{np.rad2deg(theta)*3600:.2f} arcsec, {theta*da*1e3:.1f} proper kpc')
    print(f'Catalog redshift span: {np.ptp(z)*C/(1+zmean):.1f} km/s')
    with gzip.open(HERE / 'line-pixels.csv.gz', 'rt') as stream:
        pixels = list(csv.DictReader(stream))
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.6), sharex=True, layout='constrained')
    for pos, member in enumerate(members):
        label = member['member']
        ax = axes.flat[pos]
        for line, lam, arm in LINES:
            cols = [[float(row[key]) for row in pixels if row['member'] == label and row['line'] == line]
                    for key in ('wavelength', 'flux', 'ivar')]
            w, f, iv = np.array(cols)
            expected = lam * (1 + REFERENCE_Z)
            p, model = fit_line(w, f, iv, expected)
            print(label, line, f'z={p[3]/lam-1:.5f}, FWHM={2.35482*C*p[4]/p[3]:.0f} km/s')
            if line != 'Mg II':
                continue
            # Display inverse-variance weighted bins; fit uses the original pixels.
            edges = np.arange(w.min(), w.max() + 9.6, 9.6)
            weights = np.histogram(w, edges, weights=iv)[0]
            flux = np.histogram(w, edges, weights=iv*f)[0] / np.maximum(weights, 1e-12)
            ok = weights > 0
            centers = .5 * (edges[1:] + edges[:-1])
            ax.errorbar(centers[ok], flux[ok], yerr=weights[ok]**-.5, color='#486889',
                        fmt='o', ms=2.6, lw=.8, alpha=.8)
            ax.plot(w, model(p, w), color='#a33537', lw=2)
            ax.axvline(expected, color='#555555', ls=':', lw=1)
            ax.set_title(f'{label}    DESI z = {float(member["z"]):.5f}', loc='left')
    fig.supxlabel('Observed vacuum wavelength (Å)')
    fig.supylabel('Flux density (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)')
    fig.suptitle('Broad Mg II emission from four quasar sightlines in XMM-LSS', fontsize=15)
    if options.plot:
        fig.savefig(HERE / 'mgii.png', dpi=135)
    plt.close(fig)


if __name__ == '__main__':
    main()

"""Run the Gaia candidate cross-match, spectral tests and independent light-curve checks."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import least_squares
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent
TARGET = '684208204'
GAIA = '363580870076290816'
C = 299792.458


def rows(data, name):
    with (data / name).open() as source:
        return list(csv.DictReader(source))


def amplitude(time, values, frequency):
    phase = 2 * np.pi * frequency * (time - 1700)
    design = np.column_stack((np.ones(len(phase)), np.cos(phase), np.sin(phase)))
    fit = np.linalg.lstsq(design, values, rcond=None)[0]
    return fit, design @ fit


def nightly_permutation(t, mag, frequency):
    """Shuffle night medians within observing years, preserving seasonal offsets and cadences."""
    night = np.floor(t).astype(int)
    nights = np.unique(night)
    times = np.array([np.mean(t[night == n]) for n in nights])
    values = np.array([np.median(mag[night == n]) for n in nights])
    years = np.floor((times + 2457000 - 2451544.5) / 365.25).astype(int)
    groups = [np.where(years == year)[0] for year in np.unique(years)]
    values = values - np.array([np.mean(values[years == year]) for year in years])
    phase = 2 * np.pi * frequency * (times-1700)
    design = np.column_stack((np.cos(phase), np.sin(phase)))
    for indices in groups:
        design[indices] -= design[indices].mean(axis=0)
    q, _ = np.linalg.qr(design)
    observed = sum((q.T @ values)**2)
    rng = np.random.default_rng(363580)
    exceed = 0
    for _ in range(4999):
        shuffled = values.copy()
        for ix in groups:
            shuffled[ix] = rng.permutation(values[ix])
        exceed += sum((q.T @ shuffled)**2) >= observed
    return {'nights': len(nights), 'exceed': int(exceed), 'draws': 4999}


def ccf(spectra, visit, template, template_rv, window):
    wave, signal = spectra[visit]
    wt, profile = spectra[template + '_COADD_B']
    use = (wave > window[0]) & (wave < window[1])
    wave, signal = wave[use], signal[use]
    speeds = np.arange(-320, 322, 2)
    corr = [np.corrcoef(signal, np.interp(wave / (1 + (v - template_rv) / C), wt, profile))[0, 1]
            for v in speeds]
    peaks = find_peaks(corr, distance=12, prominence=.02)[0]
    return [(int(speeds[i]), float(corr[i])) for i in sorted(peaks, key=lambda i: -corr[i])[:2]]


def broaden(spectra, visit, window):
    # A high-S/N visit of the same star supplies its own line list, widths and sampled instrumental response.
    wt, profile = spectra['1189315023_COADD']
    wave, signal = spectra[visit + '_COADD']
    use = (wave > window[0]) & (wave < window[1])
    wave, signal = wave[use], signal[use]

    def residual(parameters):
        gamma, sep, weight, scale, offset = parameters
        lower = np.interp(wave / (1 + (gamma - sep / 2) / C), wt, profile)
        upper = np.interp(wave / (1 + (gamma + sep / 2) / C), wt, profile)
        return signal - offset - scale * (weight * lower + (1 - weight) * upper)

    def fit(limit, starts):
        return min([least_squares(residual, [0, start, .5, 1, 0],
                                  bounds=([-150, 0, .2, .05, -.5], [150, limit, .8, 2, .5]))
                    for start in starts], key=lambda item: item.cost)

    single = fit(.001, [.0005])
    double = fit(400, [10, 80, 140, 200, 250])
    return {'separation_kms': round(float(double.x[1]), 2), 'rss_reduction': round(1 - double.cost / single.cost, 4)}


def analyze(data):
    screen = rows(data, 'screen.csv')
    target = next(r for r in screen if r['source_id'] == GAIA)
    findings = {'sample_counts': {s: sum(r['sample'] == s for r in screen) for s in ('MS', 'RGB')},
                'sb2_matches': [r['source_id'] for r in screen if r['kovalev_sb2'] == 'True'],
                'known_orbit_matches': [r['source_id'] for r in screen if r['guo_orbit'] == 'True'],
                'target_ms_mass_rank': 1 + sum(r['sample'] == 'MS' and float(r['fit_companion_mass']) > float(target['fit_companion_mass']) for r in screen)}
    spectra = np.load(data / 'spectra.npz')
    controls = rows(data, 'controls.csv')
    measurements = []
    for ext in ['COADD_B', 'B-84126121', 'B-84126145']:
        for ctrl in controls:
            for lo, hi in [(4900, 5155), (5160, 5300)]:
                peaks = ccf(spectra, TARGET + '_' + ext, ctrl['obsid'], float(ctrl['rv_b1']), (lo, hi))
                measurements.append({'target': TARGET, 'exposure': ext, 'template': ctrl['obsid'], 'window': [lo, hi],
                                     'velocity': [x[0] for x in peaks], 'correlation': [round(x[1], 3) for x in peaks]})
    for ctrl in controls:
        for other in controls:
            if ctrl == other:
                continue
            for lo, hi in [(4900, 5155), (5160, 5300)]:
                peaks = ccf(spectra, ctrl['obsid'] + '_COADD_B', other['obsid'], float(other['rv_b1']), (lo, hi))
                measurements.append({'target': ctrl['obsid'], 'exposure': 'COADD_B', 'template': other['obsid'], 'window': [lo, hi],
                                     'velocity': [x[0] for x in peaks], 'correlation': [round(x[1], 3) for x in peaks]})
    findings['ccf'] = measurements
    findings['low_resolution'] = {visit: {f'{lo}-{hi}': broaden(spectra, visit, (lo, hi))
                                         for lo, hi in [(3900, 4800), (4800, 5400), (5100, 5300)]}
                                  for visit in ['392902126', '491706223']}

    tess = np.genfromtxt(data / 'tess.csv', delimiter=',', names=True)
    sectors = []
    for sector in (17, 57, 84):
        d = tess[tess['sector'] == sector]
        t, v = d['btjd'], d['flux'] - np.mean(d['flux'])
        frequency, power = LombScargle(t, v).autopower(minimum_frequency=.08, maximum_frequency=6, samples_per_peak=10)
        sectors.append({'sector': sector, 'n': int(sum(d['n'])), 'bins': len(v), 'local_period': round(float(1 / frequency[np.argmax(power)]), 7)})
    t_all = tess['btjd']
    v_all = np.zeros_like(t_all)
    for rid in (17, 57, 84):
        sel = tess['sector'] == rid
        vals = tess['flux'][sel]
        v_all[sel] = (vals - vals.mean()) / vals.std()
    test_freq = np.linspace(1/.538, 1/.532, 12000)
    power = LombScargle(t_all, v_all).power(test_freq)
    peaks = sorted(find_peaks(power)[0], key=lambda i: -power[i])[:5]
    findings['tess_aliases'] = [[round(1/test_freq[i], 9), round(power[i], 4)] for i in peaks]
    tess_freq = test_freq[peaks[0]]
    ztf = rows(data, 'ztf.csv')
    ztf_result = []
    for oid in sorted({r['oid'] for r in ztf}):
        group = [r for r in ztf if r['oid'] == oid]
        mag = np.array([float(r['mag']) for r in group])
        t = np.array([float(r['hjd']) - 2457000 for r in group])
        ok = np.array([int(r['catflags']) == 0 and float(r['magerr']) < .06 for r in group]) & (abs(mag - np.median(mag)) < .15)
        t, mag = t[ok], mag[ok]
        if len(t) < 15:
            continue
        zfit, pred = amplitude(t, mag, tess_freq)
        ztf_result.append({'oid': oid, 'source': group[0]['position'], 'band': group[0]['filtercode'],
                           'n': len(t), 'semiamp_mag': round(float(np.hypot(*zfit[1:])), 6),
                           'explained_fraction': round(1 - float(np.var(mag-pred)/np.var(mag)), 4)})
        if len(t) >= 100:
            ztf_result[-1]['nightly_fixed_period_null'] = nightly_permutation(t, mag, tess_freq)
        if oid == '650115100004176':
            # Independent ground fit with g-band only, within the period window identified by individual TESS sectors.
            fgrid = np.linspace(1/.55, 1/.52, 20000)
            ls = LombScargle(t, mag)
            power = ls.power(fgrid)
            fground = fgrid[np.argmax(power)]
            findings['ztf_g_period'] = 1 / fground
            for result in sectors:
                d = tess[tess['sector'] == result['sector']]
                result['ground_prediction_flux_correlation'] = round(float(np.corrcoef(-ls.model(d['btjd'], fground), d['flux'])[0, 1]), 4)
                fit, pred = amplitude(d['btjd'], d['flux'], tess_freq)
                result['flux_semiamp'] = round(float(np.hypot(*fit[1:])), 6)
    findings['sectors'] = sectors
    findings['ztf'] = ztf_result
    return findings


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT/'data')
    parser.add_argument('--write', type=Path)
    args = parser.parse_args()
    results = analyze(args.data)
    output = json.dumps(results, indent=2) + '\n'
    if args.write:
        args.write.write_text(output)
    else:
        print(output, end='')

"""Reproduce the Gaia23bop recurrence measurements from archived survey rows."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parent


def load_measurements():
    gaia = pd.read_csv(ROOT / 'gaia.csv')
    gaia['mjd'] = gaia['JD(TCB)'] - 2400000.5
    gaia['value'] = pd.to_numeric(gaia.averagemag, errors='coerce')
    ztf = pd.read_csv(ROOT / 'ztf.csv')
    ztf = ztf.loc[(ztf.catflags == 0) & (ztf.magerr < 0.15)].copy()
    ztf['value'] = ztf.mag
    asas = pd.read_csv(ROOT / 'asassn.csv')
    asas = asas.loc[(asas.quality == 'G') & (asas.flux_err < 0.3) & (asas.flux > 0)].copy()
    asas['mjd'] = asas.jd - 2400000.5
    # Each camera has a different flux baseline in this crowded field.
    asas['value'] = asas.groupby('camera').flux.transform('median') - asas.flux
    ps1 = pd.read_csv(ROOT / 'ps1.csv')
    ps1 = ps1.loc[(ps1.psfFlux > 0) & (ps1.psfFlux / ps1.psfFluxErr > 10)].copy()
    ps1['mjd'] = ps1.obsTime
    ps1['mag'] = -2.5 * np.log10(ps1.psfFlux / 3631)
    ps1['value'] = ps1.mag - ps1.groupby('filterID').mag.transform('median')
    return {'Gaia G': gaia.dropna(subset=['value']), 'ZTF g': ztf.loc[ztf.filtercode == 'zg'],
            'ZTF r': ztf.loc[ztf.filtercode == 'zr'], 'ASAS-SN': asas, 'PS1': ps1}


def aggregate(frame, width):
    return frame.groupby(np.floor(frame.mjd / width))[['mjd', 'value']].median()


def scan(frame):
    t = frame.mjd.to_numpy()
    y = frame.value.to_numpy()
    freq = np.linspace(1 / 2000, 1 / 50, 6000)
    power = lombscargle(t - 58000, y - y.mean(), 2 * np.pi * freq, normalize=True)
    ix = np.argmax(power)
    return float(1 / freq[ix]), float(power[ix]), float(power[(1 / freq >= 210) & (1 / freq <= 240)].max())


def sine_phase(frame, period):
    angle = (frame.mjd.to_numpy() - 58000) * (2 * np.pi / period)
    coeff = np.linalg.lstsq(np.column_stack([np.ones(len(angle)), np.cos(angle), np.sin(angle)]),
                            frame.value.to_numpy(), rcond=None)[0]
    return float(58000 + np.arctan2(coeff[2], coeff[1]) * period / (2 * np.pi))


def fixed_score(frame, period, zero):
    t = frame.mjd.to_numpy()
    y = frame.value.to_numpy()
    predictor = np.cos(2 * np.pi * (t - zero) / period)
    return float(np.corrcoef(y, predictor)[0, 1])


def red_noise_controls(frame):
    """Compare the Gaia peak to equal-cadence, rank-matched OU surrogates."""
    t = frame.mjd.to_numpy()
    y = frame.value.to_numpy()
    # Frequency spacing samples the narrowest peaks several times across.
    freq = np.linspace(1 / 2000, 1 / 50, 700)
    ang = 2 * np.pi * (t - t.mean())[:, None] * freq
    # Orthogonal bases give the same zero-centered LS power as scan().
    shift = 0.5 * np.arctan2(np.sin(2 * ang).sum(0), np.cos(2 * ang).sum(0))
    c, s = np.cos(ang - shift), np.sin(ang - shift)
    c /= np.sqrt((c * c).sum(0))
    s /= np.sqrt((s * s).sum(0))

    def peak(v):
        v = v - v.mean(axis=1, keepdims=True)
        return (((v @ c) ** 2 + (v @ s) ** 2) / (v * v).sum(1)[:, None]).max(1)

    actual = peak(y[None, :])[0]
    rng = np.random.default_rng(60247)
    results = {}
    for tau in [30, 90, 180, 360, 600]:
        trials = np.empty((2000, len(t)))
        trials[:, 0] = rng.standard_normal(len(trials))
        for i, rho in enumerate(np.exp(-np.diff(t) / tau), start=1):
            trials[:, i] = rho * trials[:, i-1] + np.sqrt(1 - rho ** 2) * rng.standard_normal(len(trials))
        trials = np.sort(y)[np.argsort(np.argsort(trials, axis=1), axis=1)]
        results[str(tau)] = (1 + int((peak(trials) >= actual).sum())) / (len(trials) + 1)
    return results


def main():
    frames = load_measurements()
    grouped = {k: aggregate(v, 3 if k == 'ASAS-SN' else 1) for k, v in frames.items()}
    period = scan(grouped['Gaia G'])[0]
    zero = sine_phase(grouped['Gaia G'], period)
    alias = 1 / (1 / period + 1 / 365.25)
    alias_zero = sine_phase(grouped['Gaia G'], alias)
    result = {'Gaia_period_days': period, 'Gaia_faint_phase_MJD': zero, 'annual_alias_days': alias}
    for name, g in grouped.items():
        p, peak, alias_peak = scan(g)
        result[name] = {'input_rows': len(frames[name]), 'bins': len(g), 'measured_period_days': p,
                        'peak_power': peak, 'strongest_210_240d_power': alias_peak,
                        'fixed_Gaia_correlation': fixed_score(g, period, zero),
                        'fixed_alias_correlation': fixed_score(g, alias, alias_zero),
                        'start_mjd': g.mjd.min(), 'end_mjd': g.mjd.max()}
    # ASAS-SN camera systematics can correlate with time. Repeat using only
    # each well-sampled camera with independent temporal coverage.
    for cam, data in frames['ASAS-SN'].groupby('camera'):
        if len(data) < 100:
            continue
        g = aggregate(data, 3)
        result['camera_' + cam] = {'bins': len(g), 'fixed_Gaia_correlation': fixed_score(g, period, zero),
                                  'fixed_alias_correlation': fixed_score(g, alias, alias_zero)}
    result['Gaia_red_noise_peak_exceedance'] = red_noise_controls(grouped['Gaia G'])
    for name in ['Gaia G', 'ZTF g', 'ZTF r', 'ASAS-SN']:
        g = grouped[name]
        phase = (g.mjd.to_numpy() - zero) / period
        dist = np.abs(phase - np.rint(phase))
        y = g.value.to_numpy()
        result[name]['faint_minus_bright_median'] = float(np.median(y[dist < 0.17]) - np.median(y[dist > 0.33]))
        # Leave out each full recurrence interval, retaining only cycles with ≥3 bins.
        cycle = np.rint(phase)
        valid = [k for k in np.unique(cycle) if (cycle == k).sum() >= 3]
        result[name]['leave_one_cycle_fixed_correlation_range'] = [
            float(fn([fixed_score(g[cycle != k], period, zero) for k in valid])) for fn in [min, max]]
    # Equal survey weights keep dense blended-aperture photometry from
    # dominating the consensus. The four band-specific estimates bracket it.
    periods = np.arange(450, 751, 0.25)
    names = ['Gaia G', 'ZTF g', 'ZTF r', 'ASAS-SN']
    band_scores = []
    for name in names:
        g = grouped[name]
        band_scores.append(lombscargle(g.mjd.to_numpy() - 58000, g.value.to_numpy() - g.value.mean(),
                                       2 * np.pi / periods, normalize=True))
    consensus = float(periods[np.argmax(np.mean(band_scores, axis=0))])
    result['equal_weight_consensus_days'] = consensus
    ps = grouped['PS1']
    result['PS1']['fixed_consensus_correlation'] = fixed_score(ps, consensus, sine_phase(grouped['Gaia G'], consensus))
    # Train on Gaia up to its first three well-sampled recurrence intervals.
    # Score all subsequent Gaia nights at their actual times, with no refit.
    g = grouped['Gaia G']
    early, late = g[g.mjd < 58700], g[g.mjd >= 58700]
    p = scan(early)[0]
    result['chronological_Gaia_test'] = {'train_bins': len(early), 'test_bins': len(late),
        'training_period_days': p, 'test_correlation': fixed_score(late, p, sine_phase(early, p))}
    result['ASASSN_chronological_cross_survey'] = {}
    # A discovery curve from a separate ground survey, evaluated on later
    # Gaia and ZTF nights. Adjacent cut dates expose sensitivity to the split.
    for cutoff in [58300, 58500, 58700, 58900]:
        a = grouped['ASAS-SN']
        train = a[a.mjd < cutoff]
        p = scan(train)[0]
        epoch = sine_phase(train, p)
        scores = {}
        for name in ['ASAS-SN', 'Gaia G', 'ZTF g', 'ZTF r']:
            g = grouped[name]
            test = g[g.mjd >= cutoff]
            scores[name] = {'bins': len(test), 'correlation': fixed_score(test, p, epoch)}
        result['ASASSN_chronological_cross_survey'][str(cutoff)] = {
            'train_bins': len(train), 'training_period_days': p, 'future_fixed_scores': scores}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

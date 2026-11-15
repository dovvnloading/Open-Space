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
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

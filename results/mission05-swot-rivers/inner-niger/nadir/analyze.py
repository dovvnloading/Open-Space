"""Reproduce satellite-only flood-peak offsets and stage-excursion matchups."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
MANIFEST = json.loads((ROOT / 'provenance.json').read_text())
EVENT_DATES = [('2024-07-22', '2025-01-26'), ('2025-07-12', '2026-01-16')]

def load(product):
    path = ROOT / product['file']
    expected_hash = product.get('sha256_csv', product.get('sha256'))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_hash
    if product['id'] == 'hydro7695':
        native = pd.read_csv(path, sep=r'\s+', comment='#', header=None)
        data = pd.DataFrame({'time_utc': native[0] + 'T' + native[1],
                             'orthometric_height_m': native[2],
                             'reported_uncertainty_m': native[3]})
    else:
        data = pd.read_csv(path)
    data['time'] = pd.to_datetime(data.time_utc, utc=True)
    return data

def main():
    seasons, matches = [], []
    for product in MANIFEST:
        station = str(product['id'])
        observations = load(product)
        quality = observations[observations.reported_uncertainty_m.between(0, 0.30)].copy()
        for window, months in [('Oct-Feb', [10, 11, 12, 1, 2]), ('Nov-Jan', [11, 12, 1])]:
            sample = quality[quality.time.dt.month.isin(months)]
            year = sample.time.dt.year - (sample.time.dt.month <= 2)
            for year_id, group in sample.groupby(year):
                if year_id > 2025 or len(group) < 2:
                    continue
                peak = group.nlargest(2, 'orthometric_height_m')
                seasons.append({'station': station, 'window': window, 'year': year_id, 'n': len(group),
                                'peak_time': peak.iloc[0].time_utc, 'peak_m': peak.iloc[0].orthometric_height_m,
                                'peak_u_m': peak.iloc[0].reported_uncertainty_m,
                                'second_highest_m': peak.iloc[1].orthometric_height_m})
        if station not in ['13763', '14345', '14657', '15109', '17276']:
            continue
        for low, high in EVENT_DATES:
            dated = []
            for target in (low, high):
                day = pd.Timestamp(target, tz='UTC')
                days = (quality.time.dt.normalize() - day).dt.days
                nearest = quality.loc[days.abs().idxmin()]
                dated.append((nearest, int(days.loc[nearest.name])))
            (a, low_lag), (b, high_lag) = dated
            eligible = max(abs(low_lag), abs(high_lag)) <= 15
            matches.append(dict(station=station, low_target=low, high_target=high,
                                low_time=a.time_utc, high_time=b.time_utc,
                                low_delta_days=low_lag, high_delta_days=high_lag,
                                low_m=a.orthometric_height_m, high_m=b.orthometric_height_m,
                                low_u_m=a.reported_uncertainty_m, high_u_m=b.reported_uncertainty_m,
                                eligible=eligible, rise_m=(b.orthometric_height_m - a.orthometric_height_m)
                                if eligible else np.nan))
    seasons = pd.DataFrame(seasons)
    seasons.to_csv(ROOT / 'campaign_peaks.csv', index=False)
    pd.DataFrame(matches).to_csv(ROOT / 'stage_matchups.csv', index=False)
    results = []
    for name in ['13763', '14657', 'hydro7695']:
        for window in ['Oct-Feb', 'Nov-Jan']:
            for statistic in ['peak_m', 'second_highest_m']:
                choose = seasons[(seasons.window == window) & seasons.year.between(2016, 2024)]
                peaks = choose.pivot(index='year', columns='station', values=statistic)[[name, '14345']].dropna()
                post = (peaks.index >= 2018).astype(int)
                design = np.column_stack([np.ones(len(peaks)), peaks['14345'] - peaks['14345'].mean(), post])
                y = peaks[name].to_numpy()
                beta = np.linalg.lstsq(design, y, rcond=None)[0]
                residual = y - design @ beta
                df = len(peaks) - design.shape[1]
                se = np.sqrt((residual @ residual / df) * np.linalg.inv(design.T @ design)[2, 2])
                margin = t.ppf(0.975, df) * se
                results.append(dict(downstream=name, upstream='14345', window=window, statistic=statistic,
                                    start_year=peaks.index.min(), end_year=peaks.index.max(),
                                    pre_seasons=sum(1-post), event_seasons=sum(post), slope=beta[1],
                                    offset_cm=100*beta[2], classical_95lo_cm=100*(beta[2]-margin),
                                    classical_95hi_cm=100*(beta[2]+margin)))
    pd.DataFrame(results).to_csv(ROOT / 'conditional_offset.csv', index=False)

if __name__ == '__main__':
    main()

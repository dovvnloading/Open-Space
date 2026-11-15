"""Measure the three L1 crossings and test timing against a ground-based index."""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = pd.read_csv(ROOT / 'observations.csv.gz', index_col=0, parse_dates=True)
EARTH_RADIUS_KM = 6371.2
PROTON_MASS_KG = 1.67262192369e-27


def crossing(product, fraction=0.5, smooth_minutes=5):
    field = DATA[product + '.B'].rolling(smooth_minutes, center=True).median()
    left = field.loc['2023-10-06 20:00':'2023-10-06 20:39']
    right = field.loc['2023-10-06 22:21':'2023-10-06 23:00']
    ambient = pd.concat([left, right]).median()
    interval = field.loc['2023-10-06 21:00':'2023-10-06 22:00']
    minimum = interval.idxmin()
    below = interval < fraction * ambient
    groups = below.ne(below.shift()).cumsum()
    selected = interval[(groups == groups.loc[minimum]) & below]
    return selected.index[0], selected.index[-1], minimum, float(interval.min()), float(ambient)


def plasma(product, start, end):
    name = 'ace_swe' if product == 'ace_mag' else 'wind_pp'
    plasma_data = DATA.filter(like=name + '.').copy()
    plasma_data.columns = plasma_data.columns.str.split('.').str[-1]
    if name == 'wind_pp':
        plasma_data['Vp'] = np.sqrt(plasma_data[['VX', 'VY', 'VZ']].pow(2).sum(axis=1, min_count=3))
    plasma_data['Pdyn'] = PROTON_MASS_KG * 1e21 * plasma_data.Np * plasma_data.Vp ** 2
    minutes = pd.Timedelta(minutes=1)
    flanks = pd.concat([plasma_data.loc[start-65*minutes:start-15*minutes],
                        plasma_data.loc[end+15*minutes:end+65*minutes]])
    core = plasma_data.loc[start:end]
    return {k: {'core': float(core[k].median()), 'core_n': int(core[k].count()),
                'ambient': float(flanks[k].median()), 'ambient_n': int(flanks[k].count())}
            for k in ['Np', 'Vp', 'Pdyn']}


def calculate():
    results = {'crossings': {}, 'pair_baselines_km': {}, 'sensitivity': {}}
    positions = {}
    for name, pos in [('wind_mag', 'wind_mag'), ('ace_mag', 'ace_mag'), ('dscovr', 'dscovr_pos')]:
        start, end, minimum, bmin, ambient = crossing(name)
        mins = int((end-start).total_seconds() / 60) + 1
        xyz = DATA[[pos+'.'+c for c in ('PX','PY','PZ')]].loc[start:end].median().to_numpy()
        positions[name] = xyz * (EARTH_RADIUS_KM if name == 'wind_mag' else 1)
        row = {'start_utc': str(start), 'last_minute_utc': str(end), 'min_utc': str(minimum),
               'min_nT': bmin, 'ambient_nT': ambient, 'half_field_duration_min': mins,
               'depth_fraction': 1 - bmin/ambient}
        if name != 'dscovr':
            row['plasma'] = plasma(name, start, end)
            number = row['plasma']['Vp']['core']
            row['radial_convection_scale_km'] = number * mins * 60
            row['ballistic_arrival_at_15Re_utc'] = str(start + pd.Timedelta(seconds=(positions[name][0]-15*EARTH_RADIUS_KM)/number))
        results['crossings'][name] = row
        sens = []
        for f, window in itertools.product([.4, .5, .6], [3,5,7]):
            lo, hi, _, val, bg = crossing(name, f, window)
            sens.append({'threshold': f, 'window': window, 'duration_min': int((hi-lo).total_seconds()/60)+1,
                         'depth_fraction': 1-val/bg})
        results['sensitivity'][name] = sens
    for a,b in itertools.combinations(positions,2):
        delta = positions[a]-positions[b]
        results['pair_baselines_km'][a+'_'+b] = {'total': float(np.linalg.norm(delta)),
                                                  'transverse_to_gse_x': float(np.linalg.norm(delta[1:]))}
    h = DATA['omni.SYM_H'] if 'omni.SYM_H' in DATA else DATA['omni.H']
    # The observed core leading edge at Wind sets the time-of-flight window.
    expected = pd.Timestamp(results['crossings']['wind_mag']['ballistic_arrival_at_15Re_utc'])
    h0 = h.loc[expected-pd.Timedelta('15min'):expected].median()
    peak = h.loc[expected:expected+pd.Timedelta('30min')]
    results['ground'] = {'ambient_median_nT': float(h0), 'peak_nT': float(peak.max()),
                         'first_peak_utc': str(peak.idxmax()), 'rise_nT': float(peak.max()-h0)}
    # Ground coupling uses bow-shock-shifted OMNI field as an additional timing check.
    b = DATA['omni.B'].loc['2023-10-06 22:00':'2023-10-06 22:35']
    results['ground']['first_omni_field_below_3_5_nT_utc'] = str(b[b<3.5].index[0])
    geo = pd.read_csv(ROOT / 'goes.csv', index_col=0, parse_dates=True)
    results['geostationary'] = {}
    for instrument in ['g16', 'g18']:
        # NOAA summary flag 0 = good. The separate valid-arcjet-correction bit is informational.
        field = geo[instrument+'.B'].where(geo[instrument+'.quality'] == 0)
        baseline = field.loc[expected-pd.Timedelta('15min'):expected].median()
        window = field.loc[expected:expected+pd.Timedelta('30min')]
        results['geostationary'][instrument] = {
            'ambient_nT': float(baseline), 'peak_nT': float(window.max()),
            'peak_utc': str(window.idxmax()), 'rise_nT': float(window.max()-baseline),
            'good_minutes_in_window': int(window.count()),
            'longitude_deg': float(geo.loc[window.index,instrument+'.lon_deg'].median())}
    # A weaker feature is visible ahead of the primary depression. These windows
    # were chosen after inspection and are retained separately from the main test.
    precursor = {'upstream': {}, 'at_earth': {}}
    for instrument, row in results['crossings'].items():
        main = pd.Timestamp(row['start_utc'])
        field = DATA[instrument+'.B'].rolling(5, center=True).median()
        before = field.loc[main-pd.Timedelta('32min'):main-pd.Timedelta('8min')]
        moment = before.idxmin()
        measure = {'field_min_utc': str(moment), 'field_min_nT': float(before.min())}
        if instrument != 'dscovr':
            plasma_name = 'ace_swe' if instrument == 'ace_mag' else 'wind_pp'
            density = DATA[plasma_name+'.Np'].loc[moment-pd.Timedelta('4min'):moment+pd.Timedelta('6min')]
            measure['nearby_density_peak_cm3'] = float(density.max())
            delay = pd.Timestamp(row['ballistic_arrival_at_15Re_utc'])-main
            measure['ballistic_minimum_arrival_utc'] = str(moment+delay)
        precursor['upstream'][instrument] = measure
    arrival = pd.Timestamp(precursor['upstream']['wind_mag']['ballistic_minimum_arrival_utc'])
    series = {'SYM-H': h}
    for instrument in ['g16','g18']:
        series[instrument] = geo[instrument+'.B'].where(geo[instrument+'.quality'] == 0)
    for instrument, field in series.items():
        before = field.loc[arrival-pd.Timedelta('13min'):arrival-pd.Timedelta('3min')]
        after = field.loc[arrival:arrival+pd.Timedelta('12min')]
        precursor['at_earth'][instrument] = {'peak_utc': str(after.idxmax()), 'peak_nT': float(after.max()),
                                              'rise_nT': float(after.max()-before.median())}
    results['precursor'] = precursor
    return results


if __name__ == '__main__':
    print(json.dumps(calculate(), indent=2, allow_nan=False))

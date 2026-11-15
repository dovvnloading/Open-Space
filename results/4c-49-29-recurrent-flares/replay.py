"""Recompute the two gamma-ray localizations, field controls and optical contrasts."""
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'inputs'
RADIO = (262.68571533333335, 49.10718957777778)
REFERENCE_MJD = 54682  # FAVA week 1 began 2008-08-04, Monday 00:00 UTC
WEEK_MET = 604800


def rows(name):
    with gzip.open(DATA / (name + '.csv.gz'), 'rt') as f:
        return list(csv.DictReader(f))


def distance(ra, dec):
    delta_ra = math.radians(ra - RADIO[0])
    delta_dec = math.radians(dec - RADIO[1])
    a = math.sin(delta_dec / 2)**2 + math.cos(math.radians(dec)) * math.cos(math.radians(RADIO[1])) * math.sin(delta_ra / 2)**2
    return math.degrees(2 * math.asin(math.sqrt(a)))


def weekly_medians(points):
    bins = defaultdict(list)
    for x, mag in points:
        bins[1 + math.floor((x - REFERENCE_MJD) / 7)].append(mag)
    return {w: statistics.median(mags) for w, mags in bins.items()}


def gamma():
    apa = defaultdict(dict)
    for r in rows('fixed-apertures'):
        week = 1 + round((float(r['time']) - 239859818) / WEEK_MET)
        if week in apa[r['field']]:
            raise ValueError('Duplicate field/week')
        apa[r['field']][week] = r
    phases = {'before': (1, 863), 'between': (865, 913), 'after': (915, 951)}
    counts = {}
    fields = sorted(apa)
    for name, (low, high) in phases.items():
        counts[name] = {}
        for field in fields:
            sample = [v for k, v in apa[field].items() if low <= k <= high]
            n = sum(int(float(v['he_nev'])) for v in sample)
            mu = sum(float(v['he_avnev']) for v in sample)
            counts[name][field] = {'covered_weeks': len(sample), 'observed': n,
                                   'expected': round(mu, 4), 'ratio': round(n / mu, 4)}
    target = apa['target']
    rank = sorted(target, key=lambda k: -float(target[k]['he_sigma']))
    flares = []
    for f in json.loads((DATA / 'flares.json').read_text()):
        w = int(f['week'])
        d = distance(float(f['best_ra']), float(f['best_dec']))
        flares.append({'id': f['flareID'], 'offset_deg': round(d, 5), 'r95_deg': float(f['best_r95']),
                       'TS_low': float(f['le_ts']), 'TS_high': float(f['he_ts']),
                       'fixed_aperture_high_rank': rank.index(w) + 1,
                       'fixed_aperture_high_counts': int(float(target[w]['he_nev'])),
                       'fixed_aperture_high_expected': float(target[w]['he_avnev'])})
    return {'weeks_in_target_lightcurve': len(target), 'flares': flares, 'field_counts_above_800_MeV': counts}


def optical():
    direct = rows('ztf-direct')
    alerts = rows('ztf-alerts')
    summary = {}
    for code, fid in [('zg', '1'), ('zr', '2')]:
        good = [r for r in direct if r['filtercode'] == code and int(r['catflags']) == 0 and float(r['magerr']) < 0.15]
        # One measurement per exposure for this target. OIDs on overlapping fields can repeat a night.
        if len({r['expid'] for r in good}) != len(good):
            raise ValueError('Duplicate exposure IDs')
        weeks = weekly_medians((float(r['mjd']), float(r['mag'])) for r in good)
        pre = [value for w, value in weeks.items() if w < 864]
        candidates = [r for r in alerts if r['fid'] == fid and r['corrected'] == 'True' and r['dubious'] == 'False'
                      and r['isdiffpos'] == '1' and float(r['drb'] or 0) >= 0.9 and float(r['distnr']) < 1
                      and float(r['sigmapsf_corr_ext']) < .15 and distance(float(r['ra']), float(r['dec'])) * 3600 < 1]
        in_week = [r for r in candidates if 61073 <= float(r['mjd']) < 61080]
        calibration_offsets = [abs(float(direct_row['mag']) - float(alert['magpsf_corr']))
                               for alert in candidates if 60723 <= float(alert['mjd']) < 60730
                               for direct_row in good if abs(float(direct_row['mjd']) - float(alert['mjd'])) < 1e-4]
        if not in_week:
            raise ValueError('Missing second-flare optical observations')
        median26 = statistics.median(float(r['magpsf_corr']) for r in in_week)
        base = statistics.median(pre)
        summary[code] = {'pre_2025_weeks': len(pre), 'baseline_weekly_median_mag': round(base, 4),
                         '2025_same_week_direct_median_mag': round(weeks[864], 4),
                         'baseline_weeks_as_bright_as_2025_median': sum(x <= weeks[864] for x in pre),
                         '2025_matched_direct_alert_exposures': len(calibration_offsets),
                         '2025_max_direct_alert_absolute_difference_mag': round(max(calibration_offsets), 4),
                         '2026_same_week_alerts': len(in_week), '2026_same_week_corrected_median_mag': round(median26, 4),
                         'baseline_weeks_as_bright_as_2026_median': sum(x <= median26 for x in pre),
                         '2026_optical_flux_factor_vs_baseline': round(10**(0.4 * (base - median26)), 3),
                         '2026_observations': [{'mjd': float(r['mjd']), 'corrected_mag': float(r['magpsf_corr']),
                                                'mag_err': float(r['sigmapsf_corr_ext']), 'candid': r['candid'],
                                                'radio_separation_arcsec': round(3600 * distance(float(r['ra']), float(r['dec'])), 3)} for r in in_week]}
    return summary


if __name__ == '__main__':
    print(json.dumps({'gamma_ray': gamma(), 'optical': optical()}, indent=2))

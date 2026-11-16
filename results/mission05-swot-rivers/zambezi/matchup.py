"""Pair ZRA daily flow with instantaneous SWOT estimates at the nearest reach.

Usage: python matchup.py DATA_DIRECTORY GAUGE_CSV OUT_DIRECTORY
"""
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from datetime import date, datetime, timezone, timedelta
from pathlib import Path
from coverage import ALGORITHMS, geometry_ok, value

STATIONS = {'chavuma': 'Chavuma', 'nanas_farm': 'Nanas Farm'}
CAT = timezone(timedelta(hours=2))


def write(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main(data, gauges, out):
    out.mkdir(exist_ok=True, parents=True)
    authority = defaultdict(dict)
    with gauges.open() as f:
        for r in csv.DictReader(f):
            station, day = r['station'], date.fromisoformat(r['date'])
            if day in authority[station]:
                raise ValueError(f'Multiple daily flows at {station} {day}')
            authority[station][day] = r
    manifest = json.loads((data / 'manifest.json').read_text())
    matches, summary = [], []
    for source in manifest['hydrocron']:
        if source['rank'] != 1:
            continue
        station = STATIONS[source['site']]
        daily = authority[station]
        with (data / source['file']).open() as f:
            rows = list(csv.DictReader(f))
        for r in rows:
            if r['time_str'] == 'no_data':
                continue
            day = datetime.fromisoformat(r['time_str'].replace('Z','+00:00')).astimezone(CAT).date()
            target = min(daily, key=lambda d: (abs((d - day).days), d))
            lag = (target - day).days
            if abs(lag) > 3:
                continue
            g = daily[target]
            result = dict(site=source['site'], reach_id=r['reach_id'], time_utc=r['time_str'], swot_local_day=day,
                          gauge_date=target, gauge_minus_swot_days=lag, gauge_q_m3s=float(g['q_m3s']),
                          gauge_source_ids=g['source_ids'], cycle_id=r['cycle_id'], pass_id=r['pass_id'],
                          crid=r['crid'], reach_q=r['reach_q'], quality_primary=int(geometry_ok(r)),
                          quality_expanded=int(geometry_ok(r, 2)), granuleUR=r['granuleUR'])
            for algo in ALGORITHMS:
                q = value(r, 'sos_' + algo + '_q')
                result[algo + '_q_m3s'] = q if math.isfinite(q) else ''
            matches.append(result)
        for quality in ('primary', 'expanded'):
            for tolerance in (0, 1, 3):
                group = [r for r in matches if r['site']==source['site'] and r['quality_'+quality]
                         and abs(r['gauge_minus_swot_days'])<=tolerance]
                for algo in ALGORITHMS:
                    # One nearest date per pass, then one pass per gauge date to avoid duplicate reference days.
                    # Ties are broken chronologically. Satellite matches can only enter one day per tolerance.
                    selected = {}
                    for r in sorted(group, key=lambda x:(abs(x['gauge_minus_swot_days']), x['time_utc'])):
                        if r[algo+'_q_m3s'] != '':
                            selected.setdefault(r['gauge_date'],r)
                    ss = list(selected.values())
                    if not ss:
                        continue
                    satellite = [r[algo+'_q_m3s'] for r in ss]
                    gauge = [r['gauge_q_m3s'] for r in ss]
                    difference = [s-g for s,g in zip(satellite,gauge)]
                    ratios = [s/g for s,g in zip(satellite,gauge) if g>0]
                    corr = statistics.correlation(satellite,gauge) if len(ss)>=3 and min(gauge)<max(gauge) and min(satellite)<max(satellite) else math.nan
                    summary.append(dict(site=source['site'], algorithm=algo, quality=quality,
                                        tolerance_days=tolerance, n_pairs=len(ss), n_cycles=len({r['cycle_id'] for r in ss}),
                                        n_positive_differences=sum(d>0 for d in difference),
                                        first_date=min(selected), last_date=max(selected),
                                        median_gauge_m3s=statistics.median(gauge), median_swot_m3s=statistics.median(satellite),
                                        median_signed_error_m3s=statistics.median(difference),
                                        median_absolute_error_m3s=statistics.median(abs(d) for d in difference),
                                        median_ratio=statistics.median(ratios) if ratios else '',
                                        rmse_m3s=math.sqrt(statistics.mean(d*d for d in difference)),
                                        descriptive_r=corr if math.isfinite(corr) else ''))
    write(out/'gauge_match_candidates.csv',matches)
    write(out/'gauge_match_summary.csv',summary)


if __name__ == '__main__':
    main(*map(Path,sys.argv[1:]))

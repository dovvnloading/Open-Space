"""Tabulate SWOT observation and release-specific retrieval coverage for the 12 source files.

Usage: python coverage.py DATA_DIRECTORY OUTPUT_CSV
"""
import csv
import json
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ALGORITHMS = ('consensus', 'hivdi', 'metroman', 'momma', 'sad', 'sic4dvar', 'lakeflow')


def value(row, field):
    try:
        n = float(row[field])
    except (KeyError, ValueError):
        return math.nan
    return n if math.isfinite(n) and n > -1e10 else math.nan


def geometry_ok(row, max_quality=1):
    wse, width = value(row, 'wse'), value(row, 'width')
    return (row['time_str'] != 'no_data' and 0 <= value(row, 'reach_q') <= max_quality
            and row['partial_f'] == '0' and value(row, 'obs_frac_n') >= 0.8
            and 0 <= value(row, 'xovr_cal_q') <= 1 and math.isfinite(wse)
            and width > 0 and 0 < value(row, 'wse_u') <= 0.25
            and 0 < value(row, 'width_u') <= 0.25 * width)


def count_discharge(rows, algorithm):
    return sum(math.isfinite(value(row, f'sos_{algorithm}_q')) for row in rows)


def main(data, output):
    manifest = json.loads((data / 'manifest.json').read_text())
    features = json.loads((data / 'sword_v16_selected.geojson').read_text())['features']
    geo = {f['properties']['reach_id']: f['properties'] for f in features}
    results = []
    for source in manifest['hydrocron']:
        with (data / source['file']).open() as file:
            rows = list(csv.DictReader(file))
        reach = source['reach_id']
        count = Counter(row['reach_q'] for row in rows)
        strict = [row for row in rows if geometry_ok(row)]
        expanded = [row for row in rows if geometry_ok(row, 2)]
        times = sorted(date.fromisoformat(row['time_str'][:10]) for row in strict)
        gaps = [(b - a).days for a, b in zip(times, times[1:])]
        result = dict(site=source['site'], proximity_rank=source['rank'], reach_id=reach,
                      distance_to_vertex_km=round(geo[reach]['distance_km'], 3), n=len(rows),
                      **{f'quality_{q}': count[str(q)] for q in range(4)},
                      geom_q0=sum(geometry_ok(row, 0) for row in rows), geom_primary=len(strict),
                      geom_expanded=len(expanded), distinct_primary_cycles=len({r['cycle_id'] for r in strict}),
                      distinct_primary_pass_ids=len({r['pass_id'] for r in strict}),
                      primary_first=min(times, default=''), primary_last=max(times, default=''),
                      primary_longest_gap_days=max(gaps, default=''),
                      primary_resolved_slope=sum(0 < value(r, 'slope_u') < abs(value(r, 'slope')) / 2
                                                 for r in strict),
                      primary_neg_resolved_slope=sum(0 < value(r, 'slope_u') < -value(r, 'slope') / 2
                                                     for r in strict))
        for algo in ALGORITHMS:
            result[algo + '_all'] = count_discharge(rows, algo)
            result[algo + '_primary'] = count_discharge(strict, algo)
        results.append(result)
    with output.open('w') as file:
        writer = csv.DictWriter(file, fieldnames=results[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(results)


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))

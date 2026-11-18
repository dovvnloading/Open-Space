"""Independent NDBC wave-height comparison for the frozen tandem protocol.

Inputs are compact arrays with buoy timestamps in seconds since 1985-01-01 UTC.
Satellite inputs are already restricted to valid paired wave heights.
"""
from datetime import datetime, timezone

import numpy as np


_EARTH_KM = 6371.0
_EPOCH = datetime(1985, 1, 1, tzinfo=timezone.utc)
_HOLDOUT = (datetime(2026, 11, 1, tzinfo=timezone.utc) - _EPOCH).total_seconds()
_TRANSITION = (datetime(2026, 10, 30, tzinfo=timezone.utc) - _EPOCH).total_seconds()


def _distance_km(latlon, station):
    """Great-circle distance from each pair midpoint to one station."""
    radians = np.pi / 180.0
    lat = latlon[:, 0] * radians
    lon = latlon[:, 1] * radians
    slat, slon = station * radians
    h = np.sin((lat - slat) / 2) ** 2 + np.cos(lat) * np.cos(slat) * np.sin((lon - slon) / 2) ** 2
    return 2 * _EARTH_KM * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


def _period(t):
    return 'reference' if t < _TRANSITION else ('transition' if t < _HOLDOUT else 'holdout')


def _summary(records):
    result = {'station_visits': len(records), 'stations': len({r['station_id'] for r in records})}
    for key in ('a_signed_error_m', 'b_signed_error_m', 'a_absolute_error_m',
                'b_absolute_error_m', 'b_minus_a_absolute_error_m'):
        values = [r[key] for r in records]
        result[key + '_median'] = float(np.median(values)) if values else None
        result[key + '_quartiles'] = [float(x) for x in np.percentile(values, [25, 75])] if values else None
    result['minimum_support_met'] = result['stations'] >= 3 and len(records) >= 20
    return result


def analyze_buoys(inputs, pair_time, pair_latlon, visit_id, a_swh, b_swh):
    """Match satellite pairs to the nearest valid buoy time at each fixed station.

    Return JSON-serializable coverage, station-visit records and paired summaries
    for 50 km and, independently, 25 km. Absolute errors are calculated for
    each matched sample, then the A and B medians are subtracted within a visit.
    Ties in buoy time use the earlier timestamp; duplicate timestamps use the
    first input row. Samples may match more than one geographically close site.
    """
    ids = np.asarray(inputs['station_ids']).astype(str)
    sites = np.asarray(inputs['station_latlon'], dtype=float)
    buoy = np.asarray(inputs['buoy'], dtype=float).reshape(-1, 3)
    times = np.asarray(pair_time, dtype=float)
    coords = np.asarray(pair_latlon, dtype=float)
    visits = np.asarray(visit_id)
    a, b = np.asarray(a_swh, dtype=float), np.asarray(b_swh, dtype=float)
    n = len(times)
    if sites.shape != (len(ids), 2) or coords.shape != (n, 2) or any(len(x) != n for x in (visits, a, b)):
        raise ValueError('inconsistent station or satellite array shapes')
    if not (np.isfinite(times).all() and np.isfinite(coords).all() and np.isfinite(a).all() and np.isfinite(b).all() and (a >= 0).all() and (b >= 0).all()):
        raise ValueError('satellite inputs must have valid coordinates, times, and nonnegative SWH')

    rows = {50: [], 25: []}
    coverage = {radius: {'satellite_pairs': n, 'stations_with_coordinates': 0,
                         'stations_with_valid_buoy_rows': 0, 'station_pair_spatial_candidates': 0,
                         'station_pair_time_matches': 0, 'candidate_station_visits': 0,
                         'supported_station_visits': 0, 'stations': []} for radius in (50, 25)}
    for i, station_id in enumerate(ids):
        valid_site = bool(np.isfinite(sites[i]).all() and abs(sites[i, 0]) <= 90 and abs(sites[i, 1]) <= 180)
        selected = buoy[(buoy[:, 0] == i) & np.isfinite(buoy[:, 1]) & np.isfinite(buoy[:, 2]) & (buoy[:, 2] >= 0)]
        if len(selected):
            selected = selected[np.argsort(selected[:, 1], kind='stable')]
            _, first = np.unique(selected[:, 1], return_index=True)
            selected = selected[np.sort(first)]
        dist = _distance_km(coords, sites[i]) if valid_site else np.full(n, np.inf)
        near_idx = np.full(n, -1, dtype=int)
        dt = np.full(n, np.inf)
        if len(selected) and n:
            idx = np.searchsorted(selected[:, 1], times, side='left')
            left = np.clip(idx - 1, 0, len(selected) - 1)
            right = np.clip(idx, 0, len(selected) - 1)
            choose_right = np.abs(selected[right, 1] - times) < np.abs(selected[left, 1] - times)
            near_idx = np.where(choose_right, right, left)
            dt = np.abs(selected[near_idx, 1] - times)
        for radius in (50, 25):
            spatial = dist <= radius
            matched = spatial & (dt <= 1800)
            c = coverage[radius]
            c['stations_with_coordinates'] += int(valid_site)
            c['stations_with_valid_buoy_rows'] += int(len(selected) > 0)
            c['station_pair_spatial_candidates'] += int(spatial.sum())
            c['station_pair_time_matches'] += int(matched.sum())
            station_records = []
            for visit in np.unique(visits[matched]):
                mask = matched & (visits == visit)
                c['candidate_station_visits'] += 1
                if mask.sum() < 3:
                    continue
                ai, bi = a[mask] - selected[near_idx[mask], 2], b[mask] - selected[near_idx[mask], 2]
                a_abs, b_abs = float(np.median(np.abs(ai))), float(np.median(np.abs(bi)))
                midpoint = float(np.median(times[mask]))
                visit_value = visit.item() if isinstance(visit, np.generic) else visit
                record = {'station_id': str(station_id), 'visit_id': visit_value,
                          'period': _period(midpoint), 'pair_midpoint_time_seconds_since_1985': midpoint,
                          'paired_samples': int(mask.sum()), 'median_distance_km': float(np.median(dist[mask])),
                          'max_distance_km': float(np.max(dist[mask])), 'median_time_offset_seconds': float(np.median(dt[mask])),
                          'a_signed_error_m': float(np.median(ai)), 'b_signed_error_m': float(np.median(bi)),
                          'a_absolute_error_m': a_abs, 'b_absolute_error_m': b_abs,
                          'b_minus_a_absolute_error_m': b_abs - a_abs}
                rows[radius].append(record)
                station_records.append(record)
            c['supported_station_visits'] += len(station_records)
            c['stations'].append({'station_id': str(station_id), 'coordinates_available': valid_site,
                                  'valid_buoy_rows': int(len(selected)), 'spatial_candidates': int(spatial.sum()),
                                  'time_matches': int(matched.sum()), 'supported_station_visits': len(station_records)})
    return {f'{radius}km': {'coverage': coverage[radius], 'station_visits': rows[radius],
                            'summary': {period: _summary([r for r in rows[radius] if r['period'] == period])
                                        for period in ('reference', 'transition', 'holdout')},
                            'overall_summary': _summary(rows[radius])} for radius in (50, 25)}

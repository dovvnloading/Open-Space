"""Reconstruct the light-emitting centroids of GLM groups from two satellites."""
import csv
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
EQUATOR = 6378.137
POLE = 6356.75231414
FIGURES = (EQUATOR + 14., 6362.755)
WASHINGTON = '2026-08-14 07:48:36'


def unit(lat, lon):
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def cloud_point(lat, lon):
    # GLM L2 latitudes encode geodetic latitude on GRS80. Intersect that
    # geocentric direction with the operational lightning ellipsoid.
    phi = np.arctan((POLE / EQUATOR)**2 * np.tan(np.deg2rad(lat)))
    radius = 1. / np.hypot(np.cos(phi) / FIGURES[0], np.sin(phi) / FIGURES[1])
    return (unit(np.rad2deg(phi), lon) * radius).T


def geodetic(point):
    x, y, z = point
    rho = np.hypot(x, y)
    e2 = 1. - (POLE / EQUATOR)**2
    phi = np.arctan2(z, rho * (1. - e2))
    for _ in range(6):
        n = EQUATOR / np.sqrt(1. - e2 * np.sin(phi)**2)
        phi = np.arctan2(z + e2 * n * np.sin(phi), rho)
    return np.rad2deg(phi), np.rad2deg(np.arctan2(y, x)), rho / np.cos(phi) - n


def cross_rays(sats, intercepts):
    rays = [(pt - sat) / np.linalg.norm(pt - sat) for pt, sat in zip(intercepts, sats)]
    distances = np.linalg.lstsq(np.column_stack((rays[0], -rays[1])), sats[1] - sats[0], rcond=None)[0]
    first, second = [sat + dist * ray for sat, dist, ray in zip(sats, distances, rays)]
    return (*geodetic((first + second) / 2.), np.linalg.norm(first - second))


def extract_bin(satdata, start, width, shift=0., good_only=False):
    results = []
    for index, rows in enumerate(satdata):
        if not rows:
            return None
        offset = shift if index else 0.
        take = [r for r in rows if start <= r['seconds'] + offset < start + width
                and (not good_only or r['quality'] == 0.)]
        if not take:
            return None
        energy = np.array([r['joules'] for r in take])
        lat = np.array([r['lat'] for r in take])
        lon = np.array([r['lon'] for r in take])
        results.append((np.average(cloud_point(lat, lon), axis=0, weights=energy), sum(energy), len(energy)))
    return results


def scan(events, sats, start_grid, width, shift=0., good_only=False):
    rows = []
    for start in start_grid:
        found = extract_bin(events, start, width, shift, good_only)
        if found:
            pts, energies, counts = zip(*found)
            rows.append(dict(start=float(start), width=width, e18=float(energies[0]), e19=float(energies[1]),
                             n18=counts[0], n19=counts[1], solution=cross_rays(sats, pts)))
    return rows


def brightest(rows):
    return max(rows, key=lambda row: row['e18'] * row['e19'])


def perturb(point, east, north):
    lat, lon, _ = geodetic(point)
    up = unit(lat, lon)
    east_axis = unit(0, lon + 90.)
    return point + east * east_axis + north * np.cross(up, east_axis)


def main():
    rows = []
    with (ROOT / 'groups.csv').open() as handle:
        for r in csv.DictReader(handle):
            rows.append(dict(event=r['event'], sat=int(r['sat']),
                             **{k: float(r[k]) for k in ('seconds', 'lon', 'lat', 'joules', 'quality')}))
    granules = json.loads((ROOT / 'sources.json').read_text())['granules']
    summaries, trace = {}, []
    for date in sorted({g['event'] for g in granules}):
        meta = sorted((g for g in granules if g['event'] == date), key=lambda g: g['sat'])
        sats = [unit(g['sub_lat'], g['sub_lon']) * (EQUATOR + g['sat_height']) for g in meta]
        tracks = [[r for r in rows if r['event'] == date and r['sat'] == g['sat']] for g in meta]
        output = dict(groups=[len(track) for track in tracks], cneos_km=meta[0]['cneos_height'])
        summaries[date] = output
        # All bins have the same origin at the catalog UTC second.
        binned = scan(tracks, sats, np.arange(-75, 75) * .04, .04)
        if not binned:
            output['status'] = 'no stereo bin'
            continue
        peak = brightest(binned)
        for row in binned:
            lat, lon, h, gap = row['solution']
            trace.append([date, f"{row['start']:.3f}", row['n18'], row['n19'], f"{row['e18']:.7e}",
                          f"{row['e19']:.7e}", f'{lat:.6f}', f'{lon:.6f}', f'{h:.4f}', f'{gap:.4f}'])
        output['peak'] = peak
        if date != WASHINGTON:
            continue
        checks = []
        for w, origin, shift, good in itertools.product((.02, .04, .08, .12), (0., .5), (-.02, 0., .02), (False, True)):
            bins = scan(tracks, sats, (np.arange(-int(3/w), int(3/w)) + origin)*w, w, shift, good)
            if bins:
                checks.append(brightest(bins)['solution'][2])
        output['timing_filter_stress_km'] = [min(checks), max(checks)]
        output['timing_filter_configurations'] = len(checks)
        found = extract_bin(tracks, peak['start'], peak['width'])
        pts = [f[0] for f in found]
        for delta in (5., 10.):
            heights = []
            for e1, n1, e2, n2 in itertools.product((-delta, 0., delta), repeat=4):
                shifted = [perturb(pts[0], e1, n1), perturb(pts[1], e2, n2)]
                heights.append(cross_rays(sats, shifted)[2])
            output[f'centroid_stress_{int(delta)}km'] = [min(heights), max(heights)]
    with (ROOT / 'stereo-bins.csv').open('w', newline='') as handle:
        w = csv.writer(handle, lineterminator='\n')
        w.writerow(['utc_event', 'bin_start_offset_s', 'n18', 'n19', 'energy18_J', 'energy19_J',
                    'latitude_deg', 'longitude_deg', 'ellipsoid_height_km', 'ray_gap_km'])
        w.writerows(trace)
    (ROOT / 'reconstruction.json').write_text(json.dumps(summaries, indent=2) + '\n')


if __name__ == '__main__':
    main()

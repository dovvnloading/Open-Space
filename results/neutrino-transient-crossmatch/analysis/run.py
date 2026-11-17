"""Replay the frozen IceCube native-map search and legacy exclusion audit offline."""
from matplotlib.path import Path as PlotPath
from contours import contains, prepare
from collections import defaultdict
from datetime import datetime, timezone
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'ice'))

OUT = ROOT / 'analysis'
WINDOWS = {'optical': (0, 365), 'xray': (-1, 7), 'gamma': (-7, 7), 'mir': (-180, 180)}
TRIALS = 9999
START = '2024-01-01T00:00:00+00:00'
STOP = '2026-11-17T12:09:02+00:00'
MASKED = {'TXS0506+056', '4FGLJ0509.4+0542', 'AT2019DSG', 'AT2019FDR', 'AT2019AALC',
          'ZTF19AAPREIS', 'ZTF19AATUBSJ', 'ZTF19AAEJTOY'}


def rows(path):
    with (ROOT/path).open(newline='') as f:
        return list(csv.DictReader(f))


def days(t): return datetime.fromisoformat(t.replace('Z', '+00:00')).timestamp()/86400


def write_table(path, data, fields):
    with (OUT/path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(data)


def xyz(ra, dec):
    ra, dec = np.broadcast_arrays(np.deg2rad(ra), np.deg2rad(dec))
    return np.stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)), axis=-1)


def prepared_region(m):
    if 'uniq90' in m:
        mask = prepare(m['uniq90'])
        return lambda ra, dec: contains(mask, ra, dec)
    # Gnomonic projection maps minor great-circle edges exactly onto line segments.
    # Each released IceCat-2 contour fits well within the front hemisphere.
    loops = []
    for loop in m['polygons_deg']:
        vertices = xyz(*np.array(loop).T)
        center = vertices.mean(axis=0)
        center /= np.linalg.norm(center)
        east = np.cross([0., 0., 1.], center)
        east /= np.linalg.norm(east)
        axes = np.column_stack((east, np.cross(center, east), center))
        project = vertices@axes
        assert np.all(project[:, 2] > .5)
        path = PlotPath((project[:, :2]/project[:, 2:]).copy())
        loops.append((axes, path))

    def polygon_contains(ra, dec):
        points = xyz(ra, dec)
        hit = np.zeros(len(points), bool)
        for axes, path in loops:
            v = points@axes
            front = v[:, 2] > 0
            hit[front] |= path.contains_points(v[front, :2]/v[front, 2:])
        return hit
    return polygon_contains


CONTROL_POSITIONS = np.array([[77.358185, 5.6931481], [314.2624159, 14.2045147],
                              [257.2785752, 26.8557578], [231.0694352, 4.855293]])


def sources(archival=False):
    result = defaultdict(list)
    seen = set()
    count = defaultdict(lambda: defaultdict(int))
    paths = {'optical': 'optx/optical_events.csv', 'xray': 'optx/xray_events.csv',
             'gamma': 'gammair/lcr_flare_peaks.csv', 'mir': 'gammair/flaires_peaks.csv',
             'mir_supplement': 'gammair/supplement_2mrs_peaks.csv'}
    if archival:
        paths = {k: v for k, v in paths.items() if k in {"gamma", "mir"}}
    for family, path in paths.items():
        for r in rows(path):
            c = count[family]
            c['inventory_rows'] += 1
            epoch = r.get('utc') or r.get('peak_utc')
            day = days(epoch)
            if day > days(STOP) or ((day < days(START)) != archival):
                continue
            c['in_date_range'] += 1
            tokens = (r['source_id'], r.get('association', ''), r.get(
                'aliases', ''), r.get('published_alias', ''))
            source_vector = xyz(float(r['ra_deg']), float(r['dec_deg']))
            near_control = np.any(np.linalg.norm(xyz(*CONTROL_POSITIONS.T) -
                                  source_vector, axis=1) <= 2*np.sin(np.deg2rad(5/3600)/2))
            if near_control or any(x in z.upper().replace(' ', '') for z in tokens for x in MASKED):
                c['named_control_excluded'] += 1
                continue
            key = (family, r['source_id'], epoch)
            if key in seen:
                raise ValueError(f'Duplicate source episode {key}')
            seen.add(key)
            ra, dec = float(r['ra_deg']), float(r['dec_deg'])
            assert np.isfinite([ra, dec, day]).all() and 0 <= ra < 360 and -90 <= dec <= 90
            result[family].append(dict(source=r['source_id'], ra=ra, dec=dec, epoch=epoch, day=day, alias=r.get(
                'association') or r.get('aliases', ''), catalog=r.get('catalog') or r.get('population') or family))
            c['eligible_epochs'] += 1
    return result, dict(count)


def statistic(hits, areas):
    return np.column_stack([hits.sum(axis=1), np.where(hits, 1/areas, 0.).max(axis=1)])


def summarize(null, obs, factor):
    n = len(null)
    ge = int((null >= obs).sum())
    p = (1+ge)/(1+n)
    return dict(trials=n, exceedances=ge, median=float(np.median(null)), q95=float(np.quantile(null, .95, method='higher')),
                maximum=float(max(null)), p=p, corrected_p=min(1., factor*p), unique_statistics=len(np.unique(null)))


def main(archival=False):
    global OUT
    OUT = Path(__file__).resolve().parent/('archival' if archival else 'modern')
    OUT.mkdir(exist_ok=True)
    registry = rows('ice/gcn_alerts.csv')
    by_id = {r['run_event']: r for r in registry if r['status'] == 'active'}
    maps = [json.loads(line) for line in (ROOT/'ice/gcn_hpd90.jsonl').read_text().splitlines()]
    maps = [m for m in maps if m['run_event'] in by_id]
    assert len(set(m['run_event'] for m in maps)) == len(maps)
    for m in maps:
        m['stream'] = by_id[m['run_event']]['stream']
        m['localization'] = 'GCN HPD'
    # Include the whole first potentially eligible time-shuffle quarter.
    with gzip.open(ROOT/'ice/archival_maps.jsonl.gz', 'rt') as f:
        for line in f:
            m = json.loads(line)
            if archival or days(m['utc']) >= days('2023-07-01T00:00:00Z'):
                maps.append(m)
    # Preserve full quarter/tier strata for archival time permutations, including the late GCN era.
    maps.sort(key=lambda m: m['utc'])
    assert len({m['run_event'] for m in maps}) == len(maps)
    areas = np.array([m.get('published_area_deg2') or m['area_deg2'] for m in maps])
    assert (areas > 0).all()
    masks = [prepared_region(m) for m in maps]
    dates = np.array([days(m['utc']) for m in maps])
    n = len(maps)
    groups = defaultdict(list)
    for i, m in enumerate(maps):
        d = datetime.fromisoformat(m['utc'].replace('Z', '+00:00'))
        groups[(d.year, (d.month-1)//3+1, m['stream'].upper().split('-')[-1])].append(i)
    # Rotate the region east by delta by querying the source at RA-delta in the native cells.
    delta = np.random.default_rng((20261119 if archival else 20261117)).uniform(0, 360, size=(TRIALS, n))
    rng = np.random.default_rng((20261120 if archival else 20261118))
    perms = np.tile(np.arange(n), (TRIALS, 1))
    for ids in groups.values():
        perms[:, ids] = rng.permuted(perms[:, ids], axis=1)
    sx, counts = sources(archival)
    pairs = []
    bounds = []
    summary = {}
    hist = []
    coverage = []
    registry_exclusions = []
    families = ({k: WINDOWS[k] for k in ('gamma', 'mir')} if archival else {
                **WINDOWS, 'mir_supplement': WINDOWS['mir']})
    for family, window in families.items():
        cat = sx[family]
        source_ra = np.array([s['ra'] for s in cat])
        source_dec = np.array([s['dec'] for s in cat])
        t = np.array([s['day'] for s in cat])
        lo, hi = window
        eligible = (dates[:, None]-t >= lo) & (dates[:, None]-t <= hi)
        matched = np.zeros(n, bool)
        ra_hit = np.zeros((TRIALS, n), bool)
        time_hit = np.zeros_like(ra_hit)
        for i, (m, mask) in enumerate(zip(maps, masks)):
            sky = mask(source_ra, source_dec) if len(cat) else np.zeros(0, bool)
            hit = np.flatnonzero(sky & eligible[i])
            matched[i] = len(hit) > 0
            for j in hit:
                s = cat[j]
                pairs.append(dict(family=family, alert=m['alert_id'], source=s['source'], alias=s['alias'], catalog=s['catalog'], ra=s['ra'], dec=s['dec'], epoch=s['epoch'], alert_utc=m['utc'], lag_days=round(
                    dates[i]-t[j], 6), run_event=m['run_event'], localization=m.get('localization', 'archive'), area_deg2=areas[i]))
            lookup = np.any(eligible[:, sky], axis=1)
            time_hit[:, i] = lookup[perms[:, i]]
            # Every eligible position is queried at native HEALPix resolution in every trial.
            for j in np.flatnonzero(eligible[i]):
                ra_hit[:, i] |= mask((source_ra[j]-delta[:, i]) % 360, np.full(TRIALS, source_dec[j]))
        for a in registry:
            if a['status'] != 'active' or a['region'] == 'native_hpd90':
                continue
            needed = ('ra_deg', 'dec_deg', 'ra_plus90', 'ra_minus90', 'dec_plus90', 'dec_minus90')
            if not all(a.get(k) for k in needed):
                if family == 'optical':
                    registry_exclusions.append(a['run_event'])
                continue
            ar, ad, rp, rm, dp, dm = [float(a[k]) for k in needed]
            assert min(rp, rm, dp, dm) >= 0
            dt = days(a['utc'])-t
            dr = (source_ra-ar+180) % 360-180
            possible = (dt >= lo) & (dt <= hi) & (dr >= -rm) & (dr <=
                                                                rp) & (source_dec >= ad-dm) & (source_dec <= ad+dp)
            for j in np.flatnonzero(possible):
                s = cat[j]
                bounds.append(dict(family=family, alert=a['alert_id'], source=s['source'], alias=s['alias'], ra=s['ra'],
                              dec=s['dec'], epoch=s['epoch'], lag_days=round(dt[j], 6), circular=a['refined_circular']))
        coverage.append(dict(family=family, epochs=len(cat), native_alerts=n, alerts_with_temporal_support=int(eligible.any(axis=1).sum()), epochs_with_temporal_support=int(eligible.any(
            axis=0).sum()), time_eligible_pairs=int(eligible.sum()), in_contour_pairs=sum(p['family'] == family for p in pairs), legacy_box_unresolved=sum(p['family'] == family for p in bounds)))
        if not eligible.any() or family not in WINDOWS:
            continue
        obs = statistic(matched[None, :], areas)[0]
        rn = statistic(ra_hit, areas)
        tn = statistic(time_hit, areas)
        for k, label in enumerate(['distinct_alerts', 'maximum_inverse_area_deg-2']):
            ra = summarize(rn[:, k], obs[k], 12 if archival else 8)
            tm = summarize(tn[:, k], obs[k], 12 if archival else 8)
            summary[family+' '+label] = dict(observed=float(obs[k]), ra=ra, time=tm,
                                             larger_corrected_p=max(ra['corrected_p'], tm['corrected_p']))
            for role, v in [('ra', rn[:, k]), ('time', tn[:, k])]:
                vals, frequencies = np.unique(v, return_counts=True)
                for val, num in zip(vals, frequencies):
                    hist.append(dict(family=family, statistic=label, null=role,
                                value=float(val), frequency=int(num)))
    write_table('native_pairs.csv', pairs, ['family', 'alert', 'source', 'alias', 'catalog',
                'ra', 'dec', 'epoch', 'alert_utc', 'lag_days', 'run_event', 'localization', 'area_deg2'])
    write_table('legacy_unresolved.csv', bounds, [
                'family', 'alert', 'source', 'alias', 'ra', 'dec', 'epoch', 'lag_days', 'circular'])
    write_table('coverage.csv', coverage, list(coverage[0]))
    write_table('null_histograms.csv', hist, ['family', 'statistic', 'null', 'value', 'frequency'])
    meta = dict(cutoff=STOP, native_maps=n, scanned_registry=len(registry), unusable_legacy_ids=registry_exclusions, movable_alerts=sum(len(
        g) for g in groups.values() if len(g) > 1), strata={'/'.join(map(str, key)): len(g) for key, g in groups.items()}, source_selection=counts, tests=summary)
    (OUT/'statistics.json').write_text(json.dumps(meta, indent=2, sort_keys=True)+'\n')
    print(f'{n} maps, {len(pairs)} native pairs, {len(bounds)} unresolved bounding-box pairs')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archival', action='store_true')
    main(parser.parse_args().archival)

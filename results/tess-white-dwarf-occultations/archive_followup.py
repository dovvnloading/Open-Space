"""Pack and reconcile the short-cadence screen, pixel remeasurements and external checks."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        return json.load(stream)


def write(path, data):
    body = (json.dumps(data, separators=(',', ':'), allow_nan=False) + '\n').encode()
    path.write_bytes(gzip.compress(body, mtime=0))


def read_products(directory):
    return [read(p) for p in sorted(directory.glob('*.json'))]


def key(tic, event):
    return (str(tic), event['product'], event['kind'], event['epoch_btjd'],
            event.get('period_days'), event['duration_min'])


def event_keys(rows, product):
    return Counter(key(row['tic'], dict(ev, kind=kind, product=row[product]))
                   for row in rows for kind in ('periodic', 'isolated')
                   for ev in row.get('search', row).get(kind, []))


def audit(root):
    sc = root/'short_cadence'
    scans = read(sc/'roster_screen.json.gz')
    pixels = read(sc/'roster_pixels.json.gz')
    optical = read(root/'vetting'/'catalogue_checks.json.gz')
    with (root/'target_pixel_products.csv').open() as stream:
        inventory = list(csv.DictReader(stream))
    expected = Counter((r['tic'], int(r['sector']), int(r['cadence_s'])) for r in inventory)
    assert all(n == 1 for n in expected.values())
    assert Counter((r['tic'], r['sector'], r['cadence_s']) for r in scans) == expected
    for r in scans:
        assert 'error' not in r and 'status' in r
        if 'sha256' in r:
            assert len(r['sha256']) == 64 and r['bytes'] > 0 and r['url'].endswith(r['filename'])
    expected_events = event_keys(scans, 'filename')
    actual_events = Counter()
    for r in pixels:
        assert 'error' not in r and len(r['sha256']) == len(r['gaia_query']['sha256']) == 64
        for e in r['events']:
            actual_events[key(r['tic'], dict(e['trigger'], product=e['trigger']['lc']))] += 1
    assert actual_events == expected_events
    expected_events.update(event_keys(read(sc/'pixel_search.json'), 'filename'))
    expected_events.update(event_keys(read(root/'extension'/'extension_measurements.json'), 'lc_file'))
    for row in read(root/'record'/'historical_images.json'):
        for event in row['image_checks']:
            if event['diagnostic'] in ('individual trigger','periodic trigger'):
                item=dict(event,product=f"tic{row['tic']}-s{row['sector']}.npz",
                          kind='periodic' if event.get('period_days') else 'isolated')
                expected_events[key(row['tic'],item)] += 1
    latest = read(root/'sector110'/'tic55997020'/'search.json')
    expected_events.update(event_keys([latest], 'frames_sha256'))
    assert Counter(key(r['tic'], e) for r in optical for e in r['events']) == expected_events
    assert len({r['tic'] for r in optical}) == len(optical)
    for r in optical:
        assert set(r['queries']) == {'atlas','eb','vsx','simbad','ps','toi','gaiaepoch','ztf'}
        for v in r['queries'].values():
            assert 'error' in v or (len(v['sha256']) == 64 and v['bytes'] > 0)
    curves = {(r['tic'],r['band']):r for r in read(root/'vetting'/'optical_epochs.json.gz')}
    for r in optical:
        for band, count in r['epochs_available'].items():
            data = curves[r['tic'], band]
            assert all(len(data[field])==count for field in ('btjd','flux','flux_error','group'))
            for event in r['events']:
                p, e = event['period_days'], event['epoch_btjd']
                delta = [((t-e+.5*p)%p-.5*p if p else t-e) for t in data['btjd']]
                inside = sum(abs(d)<=event['duration_min']/2880 for d in delta)
                assert inside == event['optical'][band]['matched_points']
    print('Reconciled', len(scans), 'screened or excluded product entries,',
          sum(actual_events.values()), 'deep-screen windows,', len(optical),
          'catalogue-check hosts and', sum(expected_events.values()), 'trigger records.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('.'))
    parser.add_argument('--pack', type=Path, help='Working directory with fast, rosterpixels, catalog, and catalog_cache/')
    args = parser.parse_args()
    root = args.root
    if args.pack:
        work = args.pack
        sc = root/'short_cadence'
        vet = root/'vetting'
        vet.mkdir(exist_ok=True)
        scans = []
        for p in sorted((work/'fast').glob('*.json')):
            row = read(p)
            row['recorded_utc'] = datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat()
            scans.append(row)
        write(sc/'roster_screen.json.gz', scans)
        write(sc/'roster_pixels.json.gz', read_products(work/'rosterpixels'))
        write(vet/'catalogue_checks.json.gz', [read(p) for p in sorted((work/'catalog').glob('[0-9]*.json'))])
        (vet/'known_debris.json').write_bytes((work/'catalog'/'debris_list.json').read_bytes())
        # Preserve the matched optical epochs separately from the much wider raw sky-cone replies.
        # Reuse the service parsing and photometry cuts which produced the folds.
        from catalog_check import parse_epochs
        matches = {r['tic']:r for r in csv.DictReader((root/'tic_match.csv').open())}
        data = []
        for rec in read(vet/'catalogue_checks.json.gz'):
            tic = rec['tic']
            cache = work/'catalog_cache'
            replies = {}
            for k, receipt in rec['queries'].items():
                raw = b''
                if 'error' not in receipt:
                    raw = gzip.decompress((cache / receipt.get('cache_name', tic+'-'+k+'.gz')).read_bytes())
                    assert hashlib.sha256(raw).hexdigest() == receipt['sha256']
                replies[k] = (raw,receipt)
            _, _, epochs = parse_epochs(matches[tic], replies)
            for band, (t,fl,err,group) in epochs.items():
                assert len(t) == rec['epochs_available'][band]
                data.append(dict(tic=tic,band=band,btjd=t.tolist(),flux=fl.tolist(),
                                 flux_error=err.tolist(),group=group.tolist()))
        write(vet/'optical_epochs.json.gz', data)
    audit(root)


if __name__ == '__main__':
    main()

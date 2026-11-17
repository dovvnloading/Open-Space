"""Reacquire the first, middle and last accepted extension stars by numeric Gaia ID.

Re-extract each aperture from MAST and replay its recorded fixed- and variable-period
trials. The catalogue crossmatches are the byte-preserved survey inputs.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from ffi_acquire import download
from ffi_photometry import load_pixel_tools, run
from recovery import geometry, inject, load_visits, matches
from wd_search import search


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        return list(csv.DictReader(stream))


def replay(root, cache, pixel_tools):
    record = root / 'extension'
    measurements = json.loads((record / 'extension_measurements.json').read_text())
    accepted = sorted((r for r in measurements if r['status'] == 'included'),
                      key=lambda r: int(r['source_id']))
    roster = {r['GaiaEDR3']: r for r in table(record / 'extension_parent.csv')}
    manifests = {r['target']: r for r in table(record / 'extension_visits.csv')}
    receipts = {r['file']: r for r in json.loads((record / 'extension_receipts.json').read_text())}
    with gzip.open(record / 'extension_catalogue_responses.json.gz', 'rt') as stream:
        frozen = json.load(stream)
    fixed = table(record / 'extension_recovery_trials.csv.gz')
    continuous = table(record / 'population_trials.csv.gz')
    px = load_pixel_tools(pixel_tools)
    cache.mkdir(exist_ok=True, parents=True)
    results = []
    for index in (0, len(accepted) // 2, len(accepted) - 1):
        original = accepted[index]
        source = original['source_id']
        for suffix in ('gaia.tsv', 'tic82.tsv'):
            name = f'{source}-{suffix}'
            content = frozen[name].encode()
            if hashlib.sha256(content).hexdigest() != receipts[name + '.json']['sha256']:
                raise ValueError('Catalogue response hash mismatch: ' + name)
            (cache / name).write_bytes(content)
        name = f"{source}-s{original['sector']}.fits"
        download(cache / name, receipts[name + '.json']['url'])
        receipt = json.loads((cache / (name + '.json')).read_text())
        star = dict(roster[source], **{k: original[k] for k in ('sector','camera','ccd','coordinate_epoch')})
        fresh = run(star, cache, cache, px)
        for field in ('status', 'tic', 'arrays_sha256', 'target_fraction', 'min_target_fraction',
                      'quality0_points', 'lc_sha256'):
            if fresh[field] != original[field]:
                raise ValueError(f'{source} regenerated {field} differs from archived value')
        row = manifests[source]
        visit, = load_visits([row], cache)
        t, flux, _, dt, _ = visit
        if search(t, flux, dt) != original['search']:
            raise ValueError('Baseline blind search changed for ' + source)
        cells = {}
        for family, trials in [('fixed', fixed), ('continuous', continuous)]:
            for trial in (r for r in trials if r['target'] == source):
                p = float(trial['period_days'])
                r = float(trial['radius_earth']) if family == 'fixed' else 11.2
                epoch = float(trial['epoch_btjd'])
                a, k, b, width = geometry(p, r, float(row['mass_solar']), float(row['radius_solar']),
                                          float(trial['impact_fraction']))
                signal = inject(visit, p, epoch, a, k, b, width)
                periods = search(t, signal, dt)['periodic']
                found = next((event for event in periods if matches(event, visit, p, epoch, width)), None)
                if int(found is not None) != int(trial['recovered']):
                    raise ValueError(f"{source} {family} P={p} R={r} trial {trial['trial']} changed")
                if family == 'continuous' and found:
                    if (found['period_days'], found['epoch_btjd']) != (
                            float(trial['recovered_period']), float(trial['recovered_epoch'])):
                        raise ValueError('Blind recovered ephemeris changed for ' + source)
                cell = f'{p:g} d, {r:g} Earth radii' if family == 'fixed' else '0.25-3 d, 11.2 Earth radii'
                tally = cells.setdefault(cell, {'trials': 0, 'recovered': 0})
                tally['trials'] += 1
                tally['recovered'] += int(found is not None)
        if len(cells) != 5 or any(c['trials'] != 50 for c in cells.values()):
            raise ValueError('Incomplete trial ledger for ' + source)
        results.append(dict(source_id=source, tic=fresh['tic'], sector=fresh['sector'],
                            numeric_id_rank=index+1, original_fits_sha256=original['fits_sha256'],
                            replay_fits_receipt=receipt, verified_arrays_sha256=fresh['arrays_sha256'],
                            verified_lc_sha256=fresh['lc_sha256'], quality_zero_cadences=len(t),
                            baseline_periodic=len(original['search']['periodic']),
                            baseline_isolated=len(original['search']['isolated']), trials=cells))
        print(source, 'verified', flush=True)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('output', type=Path)
    root = Path(__file__).resolve().parent
    parser.add_argument('--pixel-tools', type=Path, default=root.parent/'tess-asymmetric-dips'/'pixels.py')
    args = parser.parse_args()
    results = replay(root, args.cache, args.pixel_tools)
    args.output.write_text(json.dumps(results, indent=2) + '\n')

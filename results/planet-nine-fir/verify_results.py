"""Audit the cross-survey counts and fingerprints from the committed measurements."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads((ROOT / path).read_bytes())


def table(path):
    with gzip.open(ROOT / path, 'rt') as stream:
        return list(csv.DictReader(stream))


def main():
    paired = read('orbits/euclid/paired_summary.json')
    draws = ROOT/'archives/planet_draws.csv.gz'
    assert hashlib.sha256(draws.read_bytes()).hexdigest() == paired['des_input_sha256']
    trial_rows = table('orbits/euclid/paired_trials.csv.gz')
    des = {(r['experiment'], int(r['draw'])): r for r in table('archives/planet_draws.csv.gz')}
    keys = {(r['experiment'], int(r['draw'])) for r in trial_rows}
    assert keys == set(des) and len(keys) == len(trial_rows) == 800
    for experiment, cases in paired['cases'].items():
        rows = [r for r in trial_rows if r['experiment'] == experiment]
        for name, reported in cases.items():
            assert all(r['DES_found_' + name] == des[(experiment, int(r['draw']))]['hit_' + name]
                       for r in rows)
            detections = [(int(r['DES_found_' + name]), int(r['Euclid_found_' + name]),
                           int(r['either_found_' + name])) for r in rows]
            assert all(e == int(bool(a or b)) for a, b, e in detections)
            assert reported['n'] == len(rows) == 400
            for k, v in [('des', sum(a for a,b,e in detections)),
                         ('euclid', sum(b for a,b,e in detections)),
                         ('both', sum(a*b for a,b,e in detections)),
                         ('either', sum(e for a,b,e in detections)),
                         ('neither', sum(1-e for a,b,e in detections))]:
                assert reported[k] == v, (experiment, name, k)
    model = 'mass5_R2.71_p0.10'
    missed_keys = {(r['experiment'], int(r['draw'])) for r in trial_rows
                   if not int(r['either_found_' + model])}
    with (ROOT/'n395_cal/paired_misses.csv').open() as stream:
        blue = list(csv.DictReader(stream))
    assert {(r['experiment'], int(r['draw'])) for r in blue} == missed_keys
    assert len(blue) == 15 and sum(int(r['eligible']) for r in blue) == 7
    assert not any(int(r['recovered'] or 0) for r in blue)
    controls = read('historical/scanpi/all_scan_tracks.json')
    assert controls['grid_size'] == 3617
    assert controls['grid_meeting_flux_threshold_in_both_groups']['0.05'] == 0
    bands = read('spherex/combined.json')
    assert bands['exposures'] == sum(row['exposures'] for row in bands['detectors'])
    assert len(bands['detectors']) == 6 and bands['grid_passes'] == bands['control_passes'] == 0
    stress = read('spherex/centroid_stress/combined.json')
    assert stress['exposures'] == bands['exposures'] == 572
    assert stress['grid_size'] == 15833 and stress['grid_valid'] == 10464
    assert stress['grid_passes'] == stress['control_passes'] == 0
    for target, eligible in [('Chen1', 679), ('Chen2', 622)]:
        result = read(f'other_optical/{target}_result.json')
        assert result['circular_tracks'] == 792 and result['at_least_two_eligible'] == eligible
    bounds = read('extra_candidates/result.json')['necessary_bound_tests']
    assert all(v['required_displacement_lower_bound_au'] > v['bound_displacement_upper_bound_au']
               for v in bounds.values())
    with gzip.open(ROOT/'wise/compact/field_results.jsonl.gz','rt') as f:
        fields = {r['target']:r for r in map(json.loads, f)}
    assert len(fields) == 46
    assert {k:r['threshold_crossings'] for k,r in fields.items() if r['threshold_crossings']} == {'RR-h270':2}
    assert sum(r['eligible'] for k,r in fields.items() if k.startswith('ACT')) == 9411
    for manifest in [ROOT/'orbits/SHA256SUMS',ROOT/'wise_vet/SHA256SUMS', ROOT/'act_vis/SHA256SUMS']:
        for line in manifest.read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            data = (manifest.parent / name.strip().lstrip('*')).read_bytes()
            assert hashlib.sha256(data).hexdigest() == digest, name
    print("Matched DES/Euclid trials, IRAS measurements, six SPHEREx detectors, "
          "46 NEOWISE field results and input manifests verified.")


if __name__ == '__main__':
    main()

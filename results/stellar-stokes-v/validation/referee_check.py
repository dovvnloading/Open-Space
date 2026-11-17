"""Independent offline cross-check of the frozen stellar Stokes V evidence.

Usage (repository checkout):
    python results/stellar-stokes-v/validation/referee_check.py
No network, numpy or astropy needed. It does not recompute archived sky-image pixels;
run ../vlass/measure.py separately for the native-pixel replay.
"""
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def table(path):
    with gzip.open(path, 'rt', newline='') if path.suffix == '.gz' else path.open('r', newline='') as f:
        return list(csv.DictReader(f))


def receipt(path):
    return json.loads(path.read_text())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha_check(path, expected, compressed=False):
    data = path.read_bytes()
    if compressed:
        data = gzip.decompress(data)
    actual = digest(data)
    assert actual == expected, (str(path), actual, expected)


def count_by_sample(rows):
    return [sum(r['sample'] == k for r in rows) for k in
            ('hosts', 'gaia_nearby', 'cns5', 'dwarfs')]


def main():
    v = ROOT / 'lofar/data'
    a = ROOT / 'askap'
    vl = ROOT / 'vlass'
    s = ROOT / 'samples'
    input_receipt = receipt(s / 'receipt.json')
    for name in ('hosts', 'dwarfs'):
        sha_check(s / f'{name}.csv.gz',
                  input_receipt['outputs'][name]['sha256'], compressed=True)
        assert len(table(s / f'{name}.csv.gz')) == input_receipt['outputs'][name]['rows']
    assert (input_receipt['outputs']['hosts']['rows'],
            input_receipt['outputs']['gaia_nearby']['rows'],
            input_receipt['outputs']['cns5']['rows'],
            input_receipt['outputs']['dwarfs']['rows']) == (4862, 48809, 5908, 372)
    assert input_receipt['hosts_with_motion'] == 4480

    lot = receipt(v / 'vlotss_receipt.json')
    sha_check(v / 'vlotss_dr2.fits', lot['input_sha256'])
    assert lot['rows'] == 68
    dr3 = receipt(v / 'dr3_pointings_receipt.json')
    sha_check(v / 'dr3_pointings.csv.gz', dr3['raw_csv_sha256'], compressed=True)
    assert len(table(v / 'dr3_pointings.csv.gz')) == dr3['rows'] == 2551
    stage = receipt(v / 'dr3_archive_status.json')
    assert len(stage) == 3 and all(r['code'] == 202 for r in stage)

    casda = receipt(a / 'casda_receipts.json')
    for item in casda:
        sha_check(a / (item['name'] + '.csv'), item['sha256'])
    assert len(table(a / 'racs_low2_v.csv')) == 221
    fp = receipt(a / 'footprint_receipt.json')
    for key, filename in [('fields', 'v_image_fields.csv.gz'),
                          ('visits', 'host_geometric_visits.csv.gz'),
                          ('hosts', 'host_coverage.csv.gz')]:
        sha_check(a / filename, fp['outputs'][key]['csv_sha256'], compressed=True)
    assert (fp['outputs']['fields']['rows'],
            fp['outputs']['visits']['rows'], fp['outputs']['hosts']['rows']) == (16503, 136617, 4862)
    assert (fp['inputs']['racs']['rows'], fp['inputs']['vast']['rows']) == (7616, 8887)
    assert fp['host_csv_sha256'] == input_receipt['outputs']['hosts']['sha256']

    cov = table(a / 'host_coverage.csv.gz')
    assert len(cov) == 4862
    assert sum(int(x['racs_distinct_sbids']) > 0 for x in cov) == 4478
    assert sum(int(x['racs_good_sbids']) > 0 for x in cov) == 4413
    assert sum(int(x['vast_distinct_sbids']) > 0 for x in cov) == 1320
    assert sum(int(x['vast_good_sbids']) > 0 for x in cov) == 1293
    assert Counter(x['status'] for x in cov) == dict(
        footprint_only=4477, outside_selected_polygons=384, published_detection=1)
    assert sum(bool(x['pixel_limit_mjy_beam']) for x in cov) == 0
    assert sum(bool(x['racs_low2_catalogue_local_5rms_mjy_beam']) for x in cov) == 1
    assert all(int(x['vast_good_sbids']) == 0 or int(x['racs_good_sbids']) > 0 for x in cov)

    lot_rows = [x for x in table(v / 'vlotss_matches.csv') if x['target_position'] == 'on']
    racs_rows = table(a / 'racs_matches.csv')
    assert len(lot_rows) == 35 and len(racs_rows) == 91
    assert count_by_sample(lot_rows) == [1, 20, 14, 0]
    assert count_by_sample(racs_rows) == [1, 51, 39, 0]
    assert len({r['vlotss_name'] for r in lot_rows}) == 21
    assert len({r['radio_id'] for r in racs_rows}) == 46
    assert table(a / 'racs_offset_matches.csv') == []
    assert count_by_sample([r for r in table(v / 'vlotss_matches.csv')
                             if r['target_position'] == 'offset']) == [0, 0, 0, 0]
    lr = receipt(v / 'vlotss_receipt.json')
    rr = receipt(a / 'racs_control_counts.json')
    assert [x['on_pairs'] for x in lr['samples']] == [1, 20, 14, 0]
    assert [x['offset_pairs'] for x in lr['samples']] == [0, 0, 0, 0]
    assert [x['target_pairs'] for x in rr['samples']] == [1, 51, 39, 0]
    assert [x['control_pairs'] for x in rr['samples']] == [0, 0, 0, 0]
    assert rr['control_positions_per_radio_record'] == 32
    known_l = {r['vlotss_name'] for r in table(HERE / 'lofar_known.csv')}
    known_r = {r['radio_id'] for r in table(HERE / 'racs_known.csv')}
    assert len(known_l) == 21 and len(known_r) == 46
    assert {r['vlotss_name'] for r in lot_rows} == known_l
    assert {r['radio_id'] for r in racs_rows} == known_r

    allmatch = table(HERE / 'match_dispositions.csv')
    hostmatch = table(HERE / 'host_v_measurements.csv')
    assert len(allmatch) == 126 and len(hostmatch) == 2
    assert Counter(x['disposition'] for x in allmatch) == {'published': 126}
    assert Counter(x['survey'] for x in allmatch) == {'V-LoTSS DR2': 35, 'RACS-low2': 91}
    assert {x['target'] for x in hostmatch} == {'GJ 1151', 'COCONUTS-2 A'}
    assert all(math.isclose(float(x['five_rms_mjy_beam']),
                            5 * float(x['local_v_rms_mjy_beam']), abs_tol=0.00000051)
               for x in allmatch)
    assert next(x for x in hostmatch if x['target'] == 'GJ 1151')['v_catalogue_value'] == '0.45935628'
    assert next(x for x in hostmatch if x['target'] == 'COCONUTS-2 A')['v_catalogue_value'] == '1.5404146'
    dr = receipt(HERE / 'disposition_receipt.json')
    assert dr['matched_catalogue_rows'] == 126 and dr['unique_racs_sources'] == 46
    assert dr['unique_vlotss_sources'] == 21
    sha_check(HERE / 'match_dispositions.csv', dr['files_sha256']['match_dispositions_csv'])
    sha_check(HERE / 'host_v_measurements.csv', dr['files_sha256']['host_v_measurements_csv'])

    mm = table(vl / 'measurements.csv')
    ctrl = table(vl / 'offsets.csv')
    independent = [r for r in mm if r['independently_counted'] == 'True']
    offset = [r for r in ctrl if r['independently_counted'] == 'True']
    vs = receipt(vl / 'summary.json')
    assert len(mm) == 13 and len(independent) == vs['independent_on'] == 11
    assert len(ctrl) == 104 and len(offset) == vs['independent_controls'] == 88
    assert all(r['sample'] != 'hosts' for r in mm)
    assert sum(abs(float(r['v_snr'])) >= 5 for r in independent) == 0
    assert sum(abs(float(r['v_snr'])) >= 5 for r in offset) == 0
    assert vs['n_5sigma_on'] == vs['n_5sigma_controls'] == 0
    assert math.isclose(max(abs(float(r['v_snr'])) for r in independent),
                        vs['max_abs_snr_on'], abs_tol=1e-10)
    ex = receipt(vl / 'extraction.json')
    assert len(ex['fields']) == vs['files'] == 12
    assert all((vl / packet['path']).exists() for f in ex['fields']
               for g in f['groups'] for packet in (g['I'], g['V']))

    srsc = table(HERE / 'srsc_hosts.csv')
    assert len(srsc) == 28
    assert len({x['hostname'] for x in srsc}) == 9
    assert len({x['hostname'] for x in srsc if x['has_v'] == 'True'}) == 5
    act = table(HERE / 'tess_activity.csv')
    assert len(act) == 6 and all(int(r['overlap_sectors']) == 0 for r in act)
    assert sum(int(r['candidate_positive_excursions']) for r in act) == len(table(HERE / 'tess_excursions.csv')) == 116

    print('PASS: published 21 V-LoTSS and 46 RACS source IDs; 126 on-pairs; 0 offsets;')
    print('PASS: 4,413 GOOD RACS and 1,293 GOOD VAST host polygons; no forced host V limits;')
    print('PASS: 0/11 VLASS target and 0/88 independent control detections;')
    print('PASS: local host V thresholds, raw catalogue/manifest SHA-256, output SHA-256, TESS timing.')


if __name__ == '__main__':
    main()

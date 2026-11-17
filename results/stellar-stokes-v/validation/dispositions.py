"""Replay the matched-source dispositions and host V measurements from frozen tables.

The status is assigned to the *radio catalogue source*, not automatically to
every Gaia component in its search circle. Undetected/uncalibrated positions
do not become limits. Does not download new archive data.
"""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIELDS = ('survey sample target catalogue_id published_counterpart date_or_mjd '
          'separation_arcsec v_catalogue_value v_unit v_catalogue_error_mjy '
          'local_v_rms_mjy_beam five_rms_mjy_beam disposition caveat').split()


def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, records, fields):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    return sha(path.read_bytes())


def main():
    sample_receipt = json.loads((ROOT / 'samples/receipt.json').read_text())
    for label in ('hosts', 'dwarfs'):
        payload = gzip.decompress((ROOT / f'samples/{label}.csv.gz').read_bytes())
        assert sha(payload) == sample_receipt['outputs'][label]['sha256']
    footprint_receipt = json.loads((ROOT / 'askap/footprint_receipt.json').read_text())
    coverage_raw = gzip.decompress((ROOT / 'askap/host_coverage.csv.gz').read_bytes())
    assert sha(coverage_raw) == footprint_receipt['outputs']['hosts']['csv_sha256']

    lofar = [r for r in read(ROOT / 'lofar/data/vlotss_matches.csv')
             if r['target_position'] == 'on']
    v_known = {r['vlotss_name']: r for r in read(HERE / 'lofar_known.csv')}
    racs = read(ROOT / 'askap/racs_matches.csv')
    r_known = {r['radio_id']: r for r in read(HERE / 'racs_known.csv')}
    assert len(lofar) == 35 and len(v_known) == 21
    assert len(racs) == 91 and len(r_known) == 46
    assert all(r['vlotss_name'] in v_known for r in lofar)
    assert all(r['radio_id'] in r_known and r_known[r['radio_id']]['source_type'] == 'star'
               for r in racs)
    assert read(ROOT / 'askap/racs_offset_matches.csv') == []

    v_receipt = json.loads((ROOT / 'lofar/data/vlotss_receipt.json').read_text())
    r_receipt = json.loads((ROOT / 'askap/racs_control_counts.json').read_text())
    assert [r['on_pairs'] for r in v_receipt['samples']] == [1, 20, 14, 0]
    assert [r['offset_pairs'] for r in v_receipt['samples']] == [0] * 4
    assert [r['target_pairs'] for r in r_receipt['samples']] == [1, 51, 39, 0]
    assert [r['control_pairs'] for r in r_receipt['samples']] == [0] * 4

    competing = defaultdict(set)
    for row in racs:
        if row['sample'] == 'gaia_nearby':
            competing[row['radio_id']].add(row['gaia_dr3_id'])
    out = []
    for r in lofar:
        caveat = ('Published accepted V-LoTSS catalogue entry; name not cross-matched '
                  'to Yiu et al. Table A.2.' if not v_known[r['vlotss_name']]['yiu_table_A2_name']
                  else 'Published V-LoTSS catalogue; repeated sample aliases are not new epochs.')
        out.append(dict(survey='V-LoTSS DR2', sample=r['sample'], target=r['target'],
            catalogue_id=r['vlotss_name'], published_counterpart=r['published_name'],
            date_or_mjd=r['date'], separation_arcsec=r['distance_arcsec'],
            v_catalogue_value=r['stokes_v_mjy'], v_unit='integrated mJy',
            v_catalogue_error_mjy=r['error_v_mjy'],
            local_v_rms_mjy_beam=r['rms_v_mjy_beam'],
            five_rms_mjy_beam=f"{5 * float(r['rms_v_mjy_beam']):.6f}",
            disposition='published', caveat=caveat))
    for r in racs:
        note = 'Published RACS-low2 source association; not an independent discovery.'
        if len(competing[r['radio_id']]) > 1:
            note += ' Multiple Gaia entries within search circle: component identity unresolved.'
        out.append(dict(survey='RACS-low2', sample=r['sample'], target=r['target'],
            catalogue_id=r['radio_id'], published_counterpart=r['source_association'],
            date_or_mjd=r['mjd'], separation_arcsec=r['separation_arcsec'],
            v_catalogue_value=r['v_peak_mjy_beam'], v_unit='peak mJy/beam',
            v_catalogue_error_mjy='', local_v_rms_mjy_beam=r['v_rms_mjy_beam'],
            five_rms_mjy_beam=f"{5 * float(r['v_rms_mjy_beam']):.6f}",
            disposition='published', caveat=note))
    match_path = HERE / 'match_dispositions.csv'
    match_sha = save(match_path, out, FIELDS)
    host = [r for r in out if r['sample'] == 'hosts']
    assert {r['target'] for r in host} == {'GJ 1151', 'COCONUTS-2 A'}
    host_path = HERE / 'host_v_measurements.csv'
    host_sha = save(host_path, host, FIELDS)

    covered = read(ROOT / 'askap/host_coverage.csv.gz')
    assert len(covered) == sample_receipt['outputs']['hosts']['rows'] == 4862
    good_racs = sum(int(r['racs_good_sbids']) > 0 for r in covered)
    good_vast = sum(int(r['vast_good_sbids']) > 0 for r in covered)
    geom_racs = sum(int(r['racs_distinct_sbids']) > 0 for r in covered)
    geom_vast = sum(int(r['vast_distinct_sbids']) > 0 for r in covered)
    assert (good_racs, good_vast, geom_racs, geom_vast) == (4413, 1293, 4478, 1320)
    assert sum(bool(r['pixel_limit_mjy_beam']) for r in covered) == 0
    assert sum(bool(r['racs_low2_catalogue_local_5rms_mjy_beam']) for r in covered) == 1

    vl = read(ROOT / 'vlass/measurements.csv')
    vn = [r for r in vl if r['independently_counted'] == 'True']
    vo = read(ROOT / 'vlass/offsets.csv')
    vsummary = json.loads((ROOT / 'vlass/summary.json').read_text())
    assert len(vl) == 13 and len(vn) == vsummary['independent_on'] == 11
    assert sum(r['sample'] == 'hosts' for r in vl) == 0
    assert len([r for r in vo if r['independently_counted'] == 'True']) == 88
    assert max(abs(float(r['v_snr'])) for r in vn) < 5
    assert vsummary['n_5sigma_on'] == vsummary['n_5sigma_controls'] == 0

    receipt = dict(matched_catalogue_rows=len(out), dispositions=dict(Counter(r['disposition'] for r in out)),
       unique_vlotss_sources=len(v_known), unique_racs_sources=len(r_known),
       vlotss_on_pairs_by_sample=[r['on_pairs'] for r in v_receipt['samples']],
       racs_on_pairs_by_sample=[r['target_pairs'] for r in r_receipt['samples']],
       vlotss_offset_matches=0, racs_offset_matches=0,
       selected_host_good_image_polygons=dict(racs=good_racs, vast=good_vast),
       selected_host_geometric_polygons=dict(racs=geom_racs, vast=geom_vast),
       hosts_with_catalogue_V_measurement=2, hosts_with_forced_pixel_V_upper_limit=0,
       vlass_distinct_nonhost_forced_positions=11, vlass_offset_positions=88,
       vlass_five_sigma_detections=0, files_sha256=dict(
            match_dispositions_csv=match_sha, host_v_measurements_csv=host_sha))
    (HERE / 'disposition_receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

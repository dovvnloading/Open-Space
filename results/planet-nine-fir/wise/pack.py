"""Bundle acquisition receipts, local searches and exact circular-orbit unions."""
import argparse
import csv
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_fields(source, name):
    return [json.loads(p.read_text()) for p in sorted((source / name).glob('*.json'))]


def write_records(path, rows):
    data = ''.join(json.dumps(x, separators=(',', ':')) + '\n' for x in rows)
    path.write_bytes(gzip.compress(data.encode(), mtime=0))


def circular_summary(results, dest):
    with (dest / 'circular_union.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['object', 'directional_fields', 'all_orbits', 'eligible_union',
                         'blank_triggers_union', 'W2_Vega', 'recovered_union'])
        for name in ('Chen1', 'Chen2', 'RR'):
            rows = [r for r in results if r['target'].split('-h')[0] == name]
            if not rows:
                continue
            assert len({r['templates'] for r in rows}) == 1
            assert len({r['target'] for r in rows}) == len(rows)
            eligible = set().union(*(r['eligible_track_ids'] for r in rows))
            triggers = set().union(*(r['trigger_track_ids'] for r in rows))
            for mag in sorted(rows[0]['injection_recovery_ids'], key=float):
                recovered = set().union(*(r['injection_recovery_ids'][mag] for r in rows))
                assert recovered <= eligible
                writer.writerow([name, len(rows), rows[0]['templates'], len(eligible),
                                 len(triggers), mag, len(recovered)])


def pack(source, dest):
    dest.mkdir(parents=True, exist_ok=True)
    receipts = read_fields(source, 'receipts')
    results = read_fields(source, 'search_results')
    circular = read_fields(source, 'circular_results')
    for name, records in [('acquisition', receipts), ('field_results', results),
                          ('circular_fields', circular)]:
        write_records(dest / (name + '.jsonl.gz'), records)
    circular_summary(circular, dest)
    by_target = {r['target']: r for r in results}
    magnitudes = [13, 14, 14.5, 15, 16]
    fields = ['target', 'band', 'science_frames_found', 'selected_frames', 'usable_frames',
              'searched_tracks', 'eligible_tracks', 'crossings', 'eligible_controls',
              'control_crossings', 'injection_draws', *[f'rec_mag{m}' for m in magnitudes],
              'max_track_score', 'status']
    with (dest / 'field_summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for a in receipts:
            r = by_target.get(a['target'])
            used = sum(all(k + '_sha256' in p for k in ['int', 'msk']) for p in a['products'])
            row = [a['target'], a.get('band', 'W2'), a['science_count'], len(a['products']), used]
            if r:
                n = {v[0]: v[1] for v in r['injections']}
                high = r['highest_mean_scores']
                row += [r[k] for k in ['templates', 'eligible', 'threshold_crossings',
                                       'eligible_controls', 'control_crossings']]
                row += [r['injections'][0][2], *[n[m] for m in magnitudes],
                        high[0]['mean'] if high else '',
                        'searched' if r['eligible'] else 'no eligible test']
            else:
                row += [''] * (len(fields) - len(row) - 1) + ['insufficient valid images or controls']
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dest', type=Path, default=HERE / 'compact')
    pack(HERE, parser.parse_args().dest)

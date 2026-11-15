"""Screen CHIME/FRB Catalog 2 for temporally separated, in-band spectral reversals."""
import argparse
import csv
from hashlib import sha256
from itertools import combinations, groupby
from math import isfinite
from pathlib import Path
from urllib.request import urlopen

CATALOG = 'https://assets.chime-frb.ca/catalogs/catalog2/releases/20251231/data/chimefrbcat2.csv'
DIGEST = '5108ada779d279a2547d9f9e73ae25bfdd40d8496d6ba7255ec29c6629057a48'


def peak_excursion(rows):
    """Largest middle-component rise above both neighbors; only resolved triplets."""
    excursions = []
    for a, b, c in zip(rows, rows[1:], rows[2:]):
        triple = (a, b, c)
        freq = [float(r['peak_freq']) for r in triple]
        if not all(410 < f < 790 for f in freq):
            continue
        for left, right in zip(triple, triple[1:]):
            dt = (float(right['mjd_inf']) - float(left['mjd_inf'])) * 86400
            widths = [float(r['width_fitb']) for r in (left, right)]
            if not all(isfinite(w) for w in widths) or not dt > 2 * max(widths):
                break
        else:
            excursions.append(freq[1] - max(freq[0], freq[2]))
    return max(excursions) if excursions else None


def save(path, header, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', type=Path, default=Path('chimefrbcat2.csv'))
    parser.add_argument('--out', type=Path, default=Path('.'))
    args = parser.parse_args()
    if not args.catalog.exists():
        args.catalog.write_bytes(urlopen(CATALOG, timeout=60).read())
    if sha256(args.catalog.read_bytes()).hexdigest() != DIGEST:
        raise ValueError('Catalog hash differs from the pinned 20251231 release.')
    with args.catalog.open() as handle:
        raw = list(csv.DictReader(handle))
    all_events = []
    eligible = []
    for name, group in groupby(sorted(raw, key=lambda r: r['tns_name']), lambda r: r['tns_name']):
        rows = list(group)
        if any(int(r['excluded_flag']) or int(r['sidelobe_flag']) for r in rows):
            continue
        rows.sort(key=lambda r: float(r['mjd_inf']))
        score = peak_excursion(rows)
        first = rows[0]
        event = (first['repeater_name'], float(first['mjd_inf']), name, len(rows), score)
        all_events.append(event)
        if score is not None:
            eligible.append(event)
    pairs = []
    for source, events in groupby(sorted(e for e in all_events if e[0]), lambda e: e[0]):
        events = list(events)
        for one, two in zip(events, events[1:]):
            if one[4] is None or two[4] is None:
                continue
            pairs.append([source, one[2], two[2], (two[1]-one[1])*86400, one[4], two[4]])
    args.out.mkdir(exist_ok=True, parents=True)
    save(args.out / 'resolved-events.csv',
         ['source','first_mjd_inf','event','components','middle_excursion_MHz'], sorted(eligible))
    save(args.out / 'adjacent-events.csv',
         ['source','first','next','seconds','first_excursion_MHz','next_excursion_MHz'], pairs)
    target = 'FRB20220912A'
    selected = sorted(e for e in eligible if e[0] == target)
    edges = {(r[1], r[2]) for r in pairs if r[0] == target and r[3] < 1800}
    score = min(e[4] for e in selected if e[2] in ('FRB20221129A','FRB20221129B'))
    # Condition on the observing times and two largest scores in this source.
    # The cut at half an hour and source were chosen after inspection.
    null = list(combinations([e[2] for e in selected], 2))
    print(f'{len(raw)} rows; {len(all_events)} qualifying events; '
          f'{sum(e[3]>=3 for e in all_events)} with >=3 components; {len(eligible)} resolved; '
          f'{len(pairs)} adjacent resolved pairs, {sum(r[3]<1800 for r in pairs)} within 30 min.')
    print(f'{target}: {len(selected)} resolved events, {sum(e[4]>=score for e in selected)} '
          f'at >= {score:.2f} MHz; {sum(pair in edges for pair in null)}/{len(null)} '
          'placements of its two highest scores make an adjacent pair within 30 min.')


if __name__ == '__main__':
    main()

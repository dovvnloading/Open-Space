"""Select well-separated early optical peaks from Shara et al. (2023), Table 3."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

CATALOG = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/269/42/table3&-out.max=unlimited'


def search(table):
    lines = iter(line for line in table.splitlines() if line and not line.startswith('#'))
    fields = next(lines).split('\t')
    next(lines)  # units
    next(lines)  # separators
    sources = defaultdict(list)
    for row in csv.DictReader(lines, fieldnames=fields, delimiter='\t'):
        sources[int(row['Seq'])].append({key: value.strip() for key, value in row.items()})
    if len(sources) != 94:
        raise ValueError('Expected 94 catalogued novae')
    hits = []
    eligible = 0
    for source, rows in sources.items():
        early = [r for r in rows if r['DF606W'] and float(r['DF606W']) <= -25]
        eligible += bool(early)
        for r in early:
            if r['l_Vmag'] or r['l_Umag'] or not r['e_Vmag'] or not r['e_Umag']:
                continue
            if float(r['e_Vmag']) >= 0.25 or float(r['e_Umag']) >= 0.4:
                continue
            lower = float(r['Vmag']) + 1
            dips = [s['Epoch'] for s in rows if s['DF606W'] and float(r['DF606W']) < float(s['DF606W']) < -5
                    and s['Vmag'] and float(s['Vmag']) > lower
                    and (s['l_Vmag'] == '>' or (s['e_Vmag'] and float(s['e_Vmag']) < .5))]
            if len(dips) >= 2:
                hits.append((source, r['Epoch'], r['DF606W'], dips))
    return sources, eligible, hits


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, help='Local unmodified VizieR ASU TSV instead of live fetch')
    args = parser.parse_args()
    data = args.table.read_text() if args.table else urlopen(CATALOG, timeout=60).read().decode()
    sources, eligible, hits = search(data)
    print(f'{len(sources)} novae; {eligible} with early time coverage; {sum(map(len, sources.values()))} rows')
    print(hits)

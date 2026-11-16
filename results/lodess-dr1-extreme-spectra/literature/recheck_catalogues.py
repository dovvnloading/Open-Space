"""Repeat the spherical nearest-neighbour audit against published selected samples."""

import argparse
import csv
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
MIRROR = 'https://vizier.idia.ac.za/viz-bin/asu-tsv'
# Release row counts include repeated measurements where a study reports spectra.
TABLES = {
    'callingham2017_pklow': ('J/ApJ/836/174/pk72', 'GLEAM', 36),
    'callingham2017_pkfreq': ('J/ApJ/836/174/pkfreq', 'GLEAM', 1222),
    'callingham2017_gps': ('J/ApJ/836/174/gps', 'GLEAM', 261),
    'callingham2017_convex': ('J/ApJ/836/174/convex', 'GLEAM', 116),
    'slob2022': ('J/A+A/668/A186/pscat', 'LoTSS', 373),
    'ballieux2024_mps': ('J/A+A/689/A264/mps-cat', 'LoTSS', 506),
    'ballieux2024_gps': ('J/A+A/689/A264/gps-cat', 'LoTSS', 8032),
    'stein2021_mhz': ('VIII/104/mhzbreak_cand', 'Name', 173703),
    'stein2021_ghz': ('VIII/104/ghzbreak_cand', 'Name', 31179),
    'utr2_mean': ('VIII/80/mean', 'Name', 2302),
}


def load_sample(key, cache):
    table, id_col, expected = TABLES[key]
    file = cache / f'{key}.csv'
    if not file.exists():
        columns = [id_col, 'RAJ2000', 'DEJ2000']
        if key.startswith('stein'):
            columns.append('cat')
        params = {'-source': table, '-out': ','.join(columns), '-out.max': 'unlimited'}
        with urlopen(MIRROR + '?' + urlencode(params), timeout=180) as response:
            lines = response.read().decode().splitlines()
        if not any(l.startswith('#Table') for l in lines):
            raise ValueError(f'Missing table {table} in VizieR response')
        rows = [l for l in lines if l.strip() and not l.startswith('#')]
        separator = next(i for i,l in enumerate(rows) if set(l) <= {'-', '\t'})
        if len(rows) - separator - 1 != expected:
            raise ValueError(f'Unexpected number of rows in {table}')
        with file.open('w', newline='') as out:
            writer = csv.writer(out)
            writer.writerow(columns)
            writer.writerows([v.strip() for v in row.split('\t')] for row in rows[separator+1:])
    with file.open() as f:
        rows = list(csv.DictReader(f))
    if len(rows) != expected:
        raise ValueError(f'Incomplete cached table {file}')
    return [(ra_degrees(r['RAJ2000']), float(r['DEJ2000']), r[id_col], r.get('cat', '')) for r in rows]


def ra_degrees(value):
    if ' ' not in value:
        return float(value)
    hours, minutes, seconds = map(float, value.split())
    return 15 * (hours + minutes / 60 + seconds / 3600)


def sep_arcsec(ra1, dec1, ra2, dec2):
    angle = math.sin(math.radians(dec2 - dec1) / 2) ** 2
    angle += math.cos(math.radians(dec1)) * math.cos(math.radians(dec2)) * math.sin(math.radians(ra2 - ra1) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, angle)))) * 3600


def nearest(sample, ra, dec):
    # Dec difference bounds angular separation. An initial real neighbour sets
    # a strict maximum radius, inside which all possible closer rows are tested.
    start = min(sample, key=lambda row: abs(dec - row[1]))
    bound = sep_arcsec(ra, dec, start[0], start[1]) / 3600 + 1e-10
    return min((sep_arcsec(ra, dec, a, d), name, subtype)
               for a, d, name, subtype in sample if abs(dec - d) <= bound)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path, help='directory outside the checkout for full catalogues')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (HERE / 'positions.csv').open() as f:
        positions = list(csv.DictReader(f))
    with args.output.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['source', 'study', 'nearest_name', 'sep_arcsec', 'prior_subtype', 'review_radius_arcsec', 'within_radius'])
        for key in [*TABLES, 'coppejans2015_table4']:
            if key == 'coppejans2015_table4':
                with (HERE / 'coppejans2015_table4.csv').open() as study:
                    sample = [(float(r['ra_deg']), float(r['dec_deg']), r['source_id'], '') for r in csv.DictReader(study)]
            else:
                sample = load_sample(key, args.cache)
            radius = 2400 if key == 'utr2_mean' else 60
            for r in positions:
                distance, name, subtype = nearest(sample, float(r['ra']), float(r['dec']))
                writer.writerow([r['id'], key, name, round(distance, 2), subtype, radius, distance <= radius])


if __name__ == '__main__':
    main()

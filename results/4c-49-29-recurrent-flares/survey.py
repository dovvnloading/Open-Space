"""Replay the exploratory sky selection. Input names and URLs are in sources.json."""
import argparse
import csv
import gzip
import collections
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent


def unit_vectors(ra, dec):
    a, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(a) * np.cos(d), np.sin(a) * np.cos(d), np.sin(d)))


def angular_query(tree, xyz):
    distance, index = tree.query(xyz)
    return np.rad2deg(2 * np.arcsin(distance / 2)), index


def archived_rows(name):
    with gzip.open(ROOT / 'inputs' / ('sky-' + name + '.csv.gz'), 'rt') as stream:
        return list(csv.DictReader(stream))


def archived_inputs():
    flares = archived_rows('flares')
    rows = archived_rows('lat')
    catalogs = []
    for name in ('4fgl', 'fl16y'):
        selected = [r for r in rows if r['catalog'] == name]
        catalogs.append(unit_vectors([float(r['ra']) for r in selected], [float(r['dec']) for r in selected]))
    radio = [(r['name'], float(r['ra']), float(r['dec'])) for r in archived_rows('radio')]
    return flares, catalogs, radio


def live_inputs(paths):
    from astropy.io import fits

    flares = json.loads(paths['flare_list'].read_text())
    catalogs = []
    for name in ('4fgl', 'fl16y'):
        with fits.open(paths[name]) as hdus:
            rows = hdus[1].data
            catalogs.append(unit_vectors(np.array(rows['RAJ2000'], dtype=float), np.array(rows['DEJ2000'], dtype=float)))
    radio = []
    for line in paths['rfc'].read_text().splitlines():
        if line.startswith('RFC '):
            fields = line.split()
            ra = 15 * sum(float(s) / 60**i for i, s in enumerate(fields[3:6]))
            sign = -1 if fields[6].startswith('-') else 1
            dec = sign * sum(abs(float(s)) / 60**i for i, s in enumerate(fields[6:9]))
            radio.append((fields[1], ra, dec))
    return flares, catalogs, radio


def survey(flares, catalogs, radio):
    selected = [r for r in flares if r['fglassoc'] == 'None' and abs(float(r['galb'])) > 10]
    positions = unit_vectors([float(r['best_ra']) for r in selected], [float(r['best_dec']) for r in selected])
    lat_distance, _ = angular_query(cKDTree(np.concatenate(catalogs)), positions)
    radio_xyz = unit_vectors([r[1] for r in radio], [r[2] for r in radio])
    radio_distance, radio_ids = angular_query(cKDTree(radio_xyz), positions)
    matches = []
    for i in np.flatnonzero((lat_distance > 0.25) & (radio_distance < 0.2)):
        matches.append({'flare': selected[i]['flareID'], 'rfc': radio[radio_ids[i]][0],
                        'radio_deg': round(float(radio_distance[i]), 5),
                        'nearest_lat_deg': round(float(lat_distance[i]), 5)})
    counts = collections.Counter(r['rfc'] for r in matches)
    return {'all_flare_rows': len(flares), 'unassociated_high_latitude_rows': len(selected),
            'rows_outside_0.25deg_catalog_radius': int((lat_distance > 0.25).sum()),
            'rfc_positions': len(radio), 'catalog_rows_4fgl_fl16y': [len(c) for c in catalogs],
            'matched_flare_rows': len(matches), 'matched_radio_sources': len(counts),
            'repeat_sources': {k: v for k, v in counts.items() if v > 1}, 'matches': matches}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--live-cache', dest='cache', type=Path, help='directory for the four downloaded source catalogs')
    args = parser.parse_args()
    if args.cache is None:
        print(json.dumps(survey(*archived_inputs()), indent=2))
        raise SystemExit(0)
    args.cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / 'sources.json').read_text())
    paths = {}
    for key in ('flare_list', '4fgl', 'fl16y', 'rfc'):
        path = args.cache / key
        if not path.exists():
            with urlopen(manifest[key]['url'], timeout=120) as response:
                path.write_bytes(response.read())
        paths[key] = path
    print(json.dumps(survey(*live_inputs(paths)), indent=2))

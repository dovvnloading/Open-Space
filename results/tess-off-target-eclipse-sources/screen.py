"""Recreate the angular, period and primary-epoch screen from the full inputs.

Usage: python screen.py tois.csv veb.dat.gz > matches.csv
The frozen queries, hashes and published prior matches are in the README.
"""
import argparse
import csv
import gzip
import itertools
import sys

import numpy as np
from scipy.spatial import cKDTree


def direction(ra, dec):
    lon = np.radians(ra)
    lat = np.radians(dec)
    return np.column_stack((np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), np.sin(lat)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tois')
    parser.add_argument('gaia')
    args = parser.parse_args()
    with open(args.tois, newline='') as stream:
        targets = [r for r in csv.DictReader(stream) if r['tfopwg_disp']=='PC'
                   and all(r[key] for key in ('ra','dec','pl_orbper','pl_tranmid'))]
    vectors = direction([float(r['ra']) for r in targets], [float(r['dec']) for r in targets])
    tree = cKDTree(vectors)
    writer = csv.writer(sys.stdout)
    writer.writerow(['toi', 'tic', 'gaia_source_id', 'sep_arcsec', 'relative_frequency_error',
                     'absolute_primary_phase_offset'])
    with gzip.open(args.gaia, 'rb') as sources:
        while rows := list(itertools.islice(sources, 60000)):
            coordinates = direction([float(r[576:597]) for r in rows],
                                    [float(r[598:620]) for r in rows])
            neighbors = tree.query_ball_point(coordinates, 2*np.sin(np.radians(30/3600)/2))
            for i, group in enumerate(neighbors):
                for j in group:
                    r, t = rows[i], targets[j]
                    angle = 2*np.arcsin(np.linalg.norm(coordinates[i]-vectors[j])/2)
                    sep = float(np.degrees(angle)*3600)
                    freq = float(r[69:90])
                    delta = float(t['pl_orbper'])*freq - 1
                    if not (3 < sep < 30 and abs(delta) < .002):
                        continue
                    elapsed = float(t['pl_tranmid']) - (2455197.5 + float(r[50:68]))
                    center = float(r[410:422].strip() or b'nan')
                    phase = abs((elapsed*freq - center + .5) % 1 - .5)
                    writer.writerow([t['toi'], t['tid'], r[:19].decode().strip(), sep, delta, phase])


if __name__ == '__main__':
    main()

"""Check all geometrically eligible compact triples in the DESI quality sample."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from survey import C, MIN_ANGLE, angular_distance, select


def counts(ra, dec, z):
    r, d = np.deg2rad([ra, dec])
    xyz = np.stack((np.cos(r)*np.cos(d), np.sin(r)*np.cos(d), np.sin(d)), axis=1)
    # 0.00022 radians exceeds 225 proper kpc / min(DA) for this sample.
    radius = .00022
    assert radius > .225 / angular_distance(z).min()
    edges = cKDTree(xyz).query_pairs(radius, output_type='ndarray')
    sep = np.linalg.norm(xyz[edges[:, 0]]-xyz[edges[:, 1]], axis=1)
    adj = defaultdict(set)
    for a, b in edges[sep > MIN_ANGLE]:
        adj[a].add(b)
    seq = np.array([(a,b,c) for a, neighbors in adj.items() for b in neighbors
                    for c in (neighbors & adj.get(b, set()))])
    diam = np.maximum.reduce([np.linalg.norm(xyz[seq[:,i]] - xyz[seq[:,j]], axis=1)
                              for i,j in combinations(range(3), 2)])
    regions = (np.floor(ra/3)*80 + np.floor((np.sin(np.deg2rad(dec))+1)/.035)).astype(int)
    order = np.argsort(regions, kind='stable')
    bounds = np.r_[0, np.flatnonzero(np.diff(regions[order]))+1, len(order)]
    pools = [order[a:b] for a,b in zip(bounds[:-1], bounds[1:])]

    def eligible(values):
        mean = values.mean(axis=1)
        return (np.ptp(values,axis=1)*C/(1+mean) <= 1000) & (diam*angular_distance(mean) <= .225)

    anchors = []
    for j in range(3):
        # Keep every qualifying real pair where the triangle fits the diameter cutoff.
        for i in np.flatnonzero(eligible(z[np.delete(seq, j, axis=1)])):
            anchors.append((i, j, seq[i,j]))
    print('Geometric triangles, observed, anchors:', len(seq), eligible(z[seq]).sum(), len(anchors))
    random = np.random.default_rng(101339)
    total, conditional = Counter(), Counter()
    for _ in range(1000):
        zz = z.copy()
        for pool in pools:
            zz[pool] = random.permutation(z[pool])
        total[int(eligible(zz[seq]).sum())] += 1
        matches = set()
        for i,j,target in anchors:
            values = z[seq[i]].copy()
            values[j] = zz[target]
            mean = values.mean()
            if np.ptp(values)*C/(1+mean) <= 1000 and diam[i]*angular_distance(mean) <= .225:
                matches.add(i)
        conditional[len(matches)] += 1
    print('All labels:', sorted(total.items()))
    print('Preserve pairs:', sorted(conditional.items()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('catalog')
    args = parser.parse_args()
    with fits.open(args.catalog) as hdus:
        table = hdus[1].data
        ix = select(table)
        counts(*[table[key][ix].astype('d') for key in ['TARGET_RA','TARGET_DEC','Z']])

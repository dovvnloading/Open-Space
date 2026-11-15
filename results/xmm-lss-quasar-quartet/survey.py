"""Replay the DESI parent-sample search and both redshift randomizations.

Usage: python survey.py /path/to/QSO_cat_iron_cumulative_v0.fits
Requires numpy, scipy, astropy. Typical peak memory is a few GB.
"""
import argparse
import itertools
from collections import Counter, defaultdict

import numpy as np
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

C = 299792.458
GRID = np.linspace(0, 5, 5001)
DISTANCE = cumulative_trapezoid(C / 70 / np.sqrt(.3 * (1 + GRID)**3 + .7), GRID, initial=0)
MIN_ANGLE = np.deg2rad(8 / 3600)


def angular_distance(z):
    return np.interp(z, GRID, DISTANCE) / (1 + z)


def select(table):
    z = table['Z']
    keep = ((table['ZWARN'] == 0) & (table['COADD_FIBERSTATUS'] == 0) &
            (np.char.strip(table['SPECTYPE']) == 'QSO') & (z > .6) & (z < 4.4) &
            (abs(table['Z_RR'] - table['Z_QN']) / (1 + z) < .015))
    lines = ['LYA', 'CIV', 'CIII', 'MgII', 'Hbeta', 'Halpha']
    count = sum((table['C_' + name] > .8).astype('u1') for name in lines)
    ix = np.flatnonzero(keep & (count >= 2))
    ix = ix[np.argsort(-table['TSNR2_QSO'][ix], kind='stable')]
    return ix[np.unique(table['TARGETID'][ix], return_index=True)[1]]


def cliques(edges):
    adjacent = defaultdict(set)
    for a, b in edges:
        adjacent[a].add(b)
    out = []
    for a, neighbors in adjacent.items():
        for b in neighbors:
            common = neighbors & adjacent.get(b, set())
            for c in common:
                out.extend((a, b, c, d) for d in common & adjacent.get(c, set()))
    return np.asarray(out, dtype=np.int32)


def broad_search(xyz, z, da, targetids):
    logz = np.log1p(z)
    width = .013
    edges = []
    for n in range(int(logz.min() / width), int(logz.max() / width) + 1):
        ids = np.flatnonzero((logz >= n*width) & (logz < (n+2)*width))
        if not len(ids):
            continue
        pairs = ids[cKDTree(xyz[ids]).query_pairs(.002, output_type='ndarray')]
        a, b = pairs.T
        good = ((abs(logz[a] - logz[b]) < 1400 / C) & (np.minimum(logz[a], logz[b]) < (n+1)*width))
        a, b = a[good], b[good]
        angle = np.linalg.norm(xyz[a] - xyz[b], axis=1)
        keep = (angle > MIN_ANGLE) & (angle * .5 * (da[a] + da[b]) < 1.8)
        edges.extend(zip(a[keep], b[keep]))
    groups = []
    for quartet in cliques(edges):
        coords = xyz[quartet]
        diameter = max(np.linalg.norm(coords[i] - coords[j]) for i, j in itertools.combinations(range(4), 2))
        groups.append((diameter * np.mean(da[quartet]) * 1000, np.ptp(logz[quartet])*C, targetids[quartet].tolist()))
    for diameter, velocity, ids in sorted(groups):
        print(f'{diameter:.1f} kpc {velocity:.1f} km/s', *ids)


def controls(xyz, z, da, ra, dec):
    # Enumerate every quartet geometrically able to satisfy the final thresholds.
    edges = cKDTree(xyz).query_pairs(1 / da.min(), output_type='ndarray')
    theta = np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1)
    seq = cliques(edges[theta >= MIN_ANGLE])
    diam = np.maximum.reduce([np.linalg.norm(xyz[seq[:, a]] - xyz[seq[:, b]], axis=1)
                              for a, b in itertools.combinations(range(4), 2)])

    def qualifying(redshifts):
        vals = redshifts[seq]
        keep = np.flatnonzero(np.ptp(vals, axis=1) * C / (1 + vals.mean(axis=1)) <= 500)
        return np.count_nonzero(diam[keep] * angular_distance(vals[keep].mean(axis=1)) <= 1)

    regions = (np.floor(ra / 3)*80 + np.floor((np.sin(np.deg2rad(dec)) + 1)/.035)).astype(int)
    order = np.argsort(regions, kind='stable')
    bounds = np.r_[0, np.flatnonzero(np.diff(regions[order]))+1, len(order)]
    blocks = [order[a:b] for a, b in zip(bounds[:-1], bounds[1:])]
    observed = qualifying(z)
    print('Geometric quartets, observed qualifying quartets:', len(seq), observed)
    rng = np.random.default_rng(11251919)
    random_counts = []
    zz = z.copy()
    for _ in range(500):
        for ids in blocks:
            zz[ids] = rng.permutation(z[ids])
        random_counts.append(qualifying(zz))
    print('All-label shuffle:', sorted(Counter(random_counts).items()))

    # Preserve each observed coherent triplet; only the candidate fourth is randomized.
    anchors = []
    for j in range(4):
        companions = z[np.delete(seq, j, axis=1)]
        mean = companions.mean(axis=1)
        # Geometry uses the four-point diameter even while identifying the anchors.
        keep = (np.ptp(companions, axis=1)*C/(1+mean) <= 500) & (diam*angular_distance(mean) <= 1)
        for i in np.flatnonzero(keep):
            anchors.append((i, j, seq[i, j]))
    targetids = np.unique([a[2] for a in anchors])
    pools = {regions[block[0]]: block for block in blocks}
    regions_used = np.unique(regions[targetids])
    rng = np.random.default_rng(290201)
    counts = []
    for _ in range(10000):
        replacements = {}
        for region in regions_used:
            targets = targetids[regions[targetids] == region]
            draws = rng.choice(pools[region], size=len(targets), replace=False)
            replacements.update(zip(targets, z[draws]))
        matches = set()
        for i, j, target in anchors:
            values = z[seq[i]].copy()
            values[j] = replacements[target]
            mean = values.mean()
            if np.ptp(values)*C/(1+mean) <= 500 and diam[i]*angular_distance(mean) <= 1:
                matches.add(i)
        counts.append(len(matches))
    print('Triplet anchors, distinct geometries:', len(anchors), len(set(a[0] for a in anchors)))
    print('Fixed-triplet completions:', sorted(Counter(counts).items()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('qso_fits')
    options = parser.parse_args()
    with fits.open(options.qso_fits, memmap=True) as hdus:
        table = hdus[1].data
        ix = select(table)
        ra, dec, z = [table[key][ix].astype('d') for key in ['TARGET_RA', 'TARGET_DEC', 'Z']]
        ids = table['TARGETID'][ix].copy()
    r, d = np.deg2rad([ra, dec])
    xyz = np.stack((np.cos(r)*np.cos(d), np.sin(r)*np.cos(d), np.sin(d)), axis=1)
    da = angular_distance(z)
    print('Unique quality-screened quasars:', len(z))
    broad_search(xyz, z, da, ids)
    controls(xyz, z, da, ra, dec)

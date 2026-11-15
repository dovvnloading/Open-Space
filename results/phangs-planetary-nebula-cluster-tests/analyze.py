"""Search the catalog cross-products. Usage: python analyze.py INPUT_DIRECTORY"""
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from scipy.spatial import cKDTree

from prepare import digest


def read_table(path, index=0):
    lines = path.read_text().splitlines()
    start = [i for i, line in enumerate(lines) if line.startswith("---")][index]
    records = itertools.takewhile(bool, lines[start + 1:])
    reader = csv.DictReader(records, fieldnames=lines[start - 2].split("\t"), delimiter="\t")
    return [{k: v.strip() for k, v in row.items()} for row in reader]


def vectors(ra, dec):
    ra, dec = np.radians([ra, dec])
    return np.stack([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)], axis=-1)


def coordinates(rows, sexagesimal):
    ra = [x['RAJ2000'] for x in rows]
    dec = [x['DEJ2000'] for x in rows]
    coord = SkyCoord(ra, dec, unit=(u.hourangle if sexagesimal else u.deg, u.deg))
    return np.column_stack((coord.ra.deg, coord.dec.deg))


def m31(rows):
    xyz = vectors(*coordinates(rows, True).T)
    ii, jj = cKDTree(xyz).query_pairs(2*np.sin(150/(2*780000)), output_type='ndarray').T
    hits = []
    for i, j in zip(ii, jj):
        a, b = rows[i], rows[j]
        if any(not row[k] for row in (a, b) for k in ('LogM', 'RVel')):
            continue
        distance = 2*780000*np.arcsin(np.linalg.norm(xyz[i]-xyz[j])/2)
        escape = np.sqrt(2*0.00430091*(10**float(a['LogM']) + 10**float(b['LogM']))/distance)
        dv = abs(float(a['RVel']) - float(b['RVel']))
        hits.append(dict(names=[a['Name'], b['Name']], separation_pc=round(float(distance), 3),
                         delta_velocity_kms=round(dv, 3), escape_kms=round(float(escape), 3)))
    return dict(total_rows=len(rows), total_pairs=len(rows)*(len(rows)-1)//2,
                under_150pc=len(ii), complete=len(hits),
                nominally_compatible=sum(x['delta_velocity_kms']<x['escape_kms'] for x in hits),
                closest=min(hits, key=lambda x:x['separation_pc']))


def clusters(directory):
    records = {}
    for filename in ('machine.fits', 'human.fits'):
        with fits.open(directory / filename) as hdus:
            for x in hdus[1].data:
                pointing = x['PHANGS_GALAXY'].strip()
                ident = int(x['ID_PHANGS_CLUSTER'])
                records[(pointing, ident)] = (pointing, ident, x['PHANGS_RA'], x['PHANGS_DEC'],
                    x['PHANGS_DR5_THILKER25_AGE'], bool(x['PHANGS_GLOBULAR_FLOYD24']))
    return list(records.values())


def split_sources(primary, alternative):
    orig = coordinates(primary, True)
    alt = coordinates(alternative, False)
    tree = cKDTree(vectors(*orig.T))
    d, _ = tree.query(vectors(*alt.T))
    sep = np.degrees(2*np.arcsin(np.minimum(1,d/2)))*3600
    good = np.array([x['c']=='pn' for x in alternative])
    confirmed_dist, _ = cKDTree(vectors(*alt[good].T)).query(vectors(*orig.T))
    agreed = np.degrees(2*np.arcsin(np.minimum(1,confirmed_dist/2)))*3600 <= .5
    result = {}
    for name, mask in [('primary', [x['Type']=='PN' for x in primary]), ('remnant', [x['Type']=='SNR' for x in primary]),
                       ('agreed', [x['Type']=='PN' and a for x,a in zip(primary,agreed)])]:
        result[name] = [(row['Galaxy'], xy) for row,xy,ok in zip(primary, orig, mask) if ok]
    result['additional'] = [(row['Galaxy'].upper(), xy) for row,xy,ok in zip(alternative,alt,good&(sep>1.0)) if ok]
    return result


def counts(rows, catalog, seed, draws):
    rng = np.random.default_rng(seed)
    radii = np.array([.1, .2, .3, .45])
    observed = np.zeros(len(radii), int)
    null = np.zeros((draws, len(radii)), int)
    for galaxy in sorted({g for g,_ in rows}):
        sky = np.array([xy for g,xy in rows if g==galaxy])
        center = sky.mean(axis=0)
        scale = np.array([3600*np.cos(np.radians(center[1])), 3600])
        pointings = galaxy.lower().replace('ngc0', 'ngc')
        local = np.array([x[2:4] for x in catalog if x[0].startswith(pointings)])
        if not local.size:
            continue
        tree = cKDTree((local-center)*scale)
        sky = (sky-center)*scale
        observed += (tree.query(sky)[0][:,None]<radii).sum(0)
        angles = rng.uniform(0, 2*np.pi, draws)
        radii_shift = np.sqrt(rng.uniform(9, 100, draws))
        shifts = radii_shift[:,None]*np.column_stack((np.cos(angles), np.sin(angles)))
        dist = tree.query((sky+shifts[:,None,:]).reshape(-1,2))[0].reshape(draws,-1)
        null += np.stack([(dist<r).sum(1) for r in radii],axis=1)
    return dict(observed=observed.tolist(), mean_shifted=null.mean(0).round(3).tolist(),
                p_marginal=((1+(null>=observed).sum(0))/(draws+1)).tolist())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', type=Path)
    parser.add_argument('--draws', type=int, default=4999)
    args = parser.parse_args()
    d = args.inputs
    for fname, expected in json.loads(Path(__file__).with_name('inputs.json').read_text()).items():
        if digest(d/fname) != expected['sha256']:
            raise ValueError(f'Input changed: {fname}')
    cat = clusters(d)
    src = split_sources(read_table(d/'pne.tsv',1), read_table(d/'congiu.tsv'))
    out = {'m31':m31(read_table(d/'m31.tsv')), 'source_counts':{k:len(v) for k,v in src.items()},
           'radii_arcsec':[.1,.2,.3,.45], 'draws':args.draws, 'seed':71412,
           'hst_union_size':len(cat), 'controls':{}}
    for name, rows in src.items():
        for cut, selected in [('all',cat),('age_30_to_3000_myr',[x for x in cat if 30<=x[4]<=3000]),
                              ('globular_flag',[x for x in cat if x[5]])]:
            out['controls'][name+'_'+cut] = counts(rows, selected, 71412, args.draws)
    print(json.dumps(out,indent=2,allow_nan=False))

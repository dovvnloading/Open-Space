"""Recompute the HD 72223 association checks from the frozen public records."""
import csv
import hashlib
import json
import math
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / 'inputs'
PRIMARY = '5707071630039205248'
OUTER = '5707071698757576576'
PARAMETERS = ('parallax', 'pmra', 'pmdec')


def table(path):
    root = ElementTree.parse(INPUTS / path)
    assert all(info.get('value') == 'OK' for info in root.findall('.//{*}INFO[@name="QUERY_STATUS"]'))
    data = root.find('.//{*}TABLE')
    columns = [item.get('name') for item in data.findall('{*}FIELD')]
    return [dict(zip(columns, [cell.text for cell in row])) for row in data.findall('.//{*}TR')]


def vizier(path):
    lines = [line for line in (INPUTS / path).read_text().splitlines() if line and not line.startswith('#')]
    names = lines[0].split('\t')
    return [dict(zip(names, [v.strip() for v in line.split('\t')])) for line in lines[3:]]


def correlation(row):
    if 'corr_vec' in row:
        # NSS DR3 upper triangle, ordered by column; ra, dec, parallax, pmra, pmdec come first.
        assert int(row['bit_index']) & 31 == 31
        entries = list(map(float, row['corr_vec'].split()))
        a, b, c = [entries[index] for index in (5, 8, 9)]
    else:
        a, b, c = [float(row[name]) for name in ('parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr')]
    sigma = np.array([float(row[key + '_error']) for key in PARAMETERS])
    return np.outer(sigma, sigma) * np.array([[1, a, b], [a, 1, c], [b, c, 1]])


def separation(primary, outer):
    ra1, dec1 = [math.radians(float(primary[k])) for k in ('ra', 'dec')]
    ra2, dec2 = [math.radians(float(outer[k])) for k in ('ra', 'dec')]
    h = math.sin((dec1-dec2)/2)**2 + math.cos(dec1)*math.cos(dec2)*math.sin((ra1-ra2)/2)**2
    return math.degrees(2*math.asin(math.sqrt(h))) * 3600


def main():
    for filename, wanted in json.loads((ROOT / 'hashes.json').read_text()).items():
        actual = hashlib.sha256((INPUTS / filename).read_bytes()).hexdigest()
        if actual != wanted:
            raise ValueError(f'Input checksum failed for {filename}')
    cohort = {r['source_id']:r['table'] for r in csv.DictReader((INPUTS / 'rekhi_ids.csv').read_text().splitlines())}
    gaia = {r['source_id']:r for r in table('gaia_sources.xml')}
    nss = {r['source_id']:r for r in table('orbital_solutions.xml') if r['parallax'] is not None}
    results = []
    for catalog in ('primary_members.tsv', 'secondary_members.tsv'):
        for row in vizier(catalog):
            ids = [row['Source1'], row['Source2']]
            primary, = [i for i in ids if i in cohort]
            outer, = [i for i in ids if i != primary]
            a, b = nss[primary], gaia[outer]
            diff = np.array([float(b[k])-float(a[k]) for k in PARAMETERS])
            cov = correlation(a) + correlation(b)
            angle = separation(a,b)
            factor = 4.74047 / float(a['parallax'])
            results.append(dict(primary=primary, outer=outer, table=cohort[primary], binary_type=row['BinType'],
                                R=float(row['R']), sep_arcsec=angle, projected_au=angle*1000/float(a['parallax']),
                                identical_kinematics_chi2=float(diff @ np.linalg.solve(cov, diff)),
                                delta_parallax_sigma=diff[0]/np.sqrt(cov[0,0]), vtan_relative_km_s=factor*np.linalg.norm(diff[1:]),
                                delta_outer_parallax_pmra_pmdec=diff.tolist()))
    a, b = nss[PRIMARY], gaia[OUTER]
    detections = list(csv.DictReader((INPUTS / 'ps1_detections.csv').read_text().splitlines()))
    motions = []
    for visit in detections:
        if visit['objID'] != '85131277522736571':
            continue
        dt = 2000 + (float(visit['obsTime'])-51544.5)/365.25 - 2016
        offset = np.array([(float(visit['ra'])-float(b['ra']))*math.cos(math.radians(float(b['dec']))), float(visit['dec'])-float(b['dec'])])*3600
        motion = np.array([float(a['pmra']),float(a['pmdec'])])*dt/1000
        motions.append((float(visit['obsTime']),float(np.linalg.norm(offset-motion)),float(np.linalg.norm(offset)),int(visit['infoFlag2'])))
    checks = {}
    for name, group in [('all',motions),('infoFlag2_zero',[r for r in motions if r[3]==0])]:
        checks[name] = dict(n=len(group),mjd_range=[min(r[0] for r in group),max(r[0] for r in group)],
                            median_comoving_arcsec=float(np.median([r[1] for r in group])),
                            median_stationary_arcsec=float(np.median([r[2] for r in group])),
                            comoving_closer=sum(r[1] < r[2] for r in group))
    host_result = next(r for r in results if r['primary'] == PRIMARY)
    return dict(n_targets=len(cohort), wide_catalog_pairs=results, ps1_astrometric_checks=checks,
                catalog_wd_probability=float(vizier('white_dwarf.tsv')[0]['Pwd']),
                field_source_ids=[r['Source'] for r in vizier('local_field.tsv')],
                inner_period_d=float(a['period']),
                absolute_G_at_orbit_distance=float(b['phot_g_mean_mag'])+5*math.log10(float(a['parallax']))-10,
                escape_at_projected_sep_km_s=42.1219*math.sqrt((.93+.60+.64)/host_result['projected_au']))


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))

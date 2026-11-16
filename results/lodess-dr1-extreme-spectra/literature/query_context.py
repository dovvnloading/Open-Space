"""Recover the archival cone queries and name-resolved bibliographies as raw tables."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree as ET

HERE = Path(__file__).resolve().parent
VIZIER = 'https://vizier.idia.ac.za/viz-bin/asu-tsv'
NED = 'https://ned.ipac.caltech.edu/NED::API/'
ASTRON = 'https://vo.astron.nl/tap/sync'
SIMBAD = 'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
CONTEXT = ['LOL1J112022.1+502936', 'LOL1J142947.2+503033', 'LOL1J144343.3+503432']
NED_NAMES = {'LOL1J112022.1+502936': 'NVSS J112022+502936',
             'LOL1J142947.2+503033': 'NVSS J142947+503033',
             'NVSS J100005+812702': 'NVSS J100005+812702'}


def request(base, params):
    with urlopen(base + '?' + urlencode(params), timeout=90) as response:
        return response.read()


def adql(base, sql):
    data = request(base, {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': sql})
    if data.lstrip().startswith(b'<'):
        raise ValueError(data[:350])
    return data.decode()


def query(target, out):
    name, ra, dec = target['id'], float(target['ra']), float(target['dec'])
    basename = name.replace(' ', '_')
    point = f"POINT('ICRS',{ra},{dec})"
    circle = f"CIRCLE('ICRS',{ra},{dec},0.016666667)"
    objects = adql(SIMBAD, f"""SELECT b.oid, b.main_id, b.otype, b.rvz_redshift,
        DISTANCE(POINT('ICRS',b.ra,b.dec),{point})*3600 AS sep_as
        FROM basic AS b WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),{circle}) ORDER BY sep_as""")
    (out / f'{basename}_simbad.csv').write_text(objects)
    for r in csv.DictReader(objects.splitlines()):
        if float(r['sep_as']) > 25:
            continue
        oid = int(r['oid'])
        for suffix, sql in {
                'ids': f'SELECT id FROM ident WHERE oidref={oid}',
                'refs': f"""SELECT r.bibcode, r.title FROM ref AS r JOIN has_ref AS h
                         ON h.oidbibref=r.oidbib WHERE h.oidref={oid}"""}.items():
            (out / f'{basename}_{oid}_{suffix}.csv').write_text(adql(SIMBAD, sql))
    if name in NED_NAMES:
        for action in ['References', 'Photometry', 'Crossids']:
            data = request(NED + action + 'OfObject', {'TARGET': NED_NAMES[name]})
            if not any(p.get('value') == 'OK' for p in ET.fromstring(data).findall(".//{*}PARAM[@name='QUERY_STATUS']")):
                raise ValueError(f'NED query failed for {name} {action}')
            (out / f'{basename}_ned_{action}.xml').write_bytes(data)
    if name in CONTEXT:
        catalogues = ['J/A+A/678/A151/catalog', 'VII/289/dr16q', 'VIII/104/spectra', 'J/ApJ/793/82']
        observations = adql(ASTRON, f"""SELECT obsid,raj2000,decj2000,obsdate,obstime,
            baselinequal,obsqual,baselinequalft,obsflag FROM lbcs.main
            WHERE 1=CONTAINS(POINT('ICRS',raj2000,decj2000),{circle})""")
        (out / f'{basename}_lbcs.csv').write_text(observations)
    elif name == 'NVSS J100005+812702':
        catalogues = ['VIII/104/spectra', 'J/MNRAS/515/1736', 'VIII/44/table2',
                      'J/A+A/370/78', 'J/A+A/463/97/table4', 'J/A+A/544/A25/sources']
        field = request(VIZIER, {'-source': 'VIII/44/table3', 'Field': '153', '-out': '**'})
        (out / 'miyun_field153.tsv').write_bytes(field)
    else:
        return
    for table in catalogues:
        data = request(VIZIER, {'-source': table, '-c': f'{ra} {dec:+.9f}',
                               '-c.rs': '30', '-out': '**', '-out.max': '1000'})
        if b'#Table' not in data or b'#INFO\tError=' in data:
            raise ValueError(f'VizieR query failed for {name} {table}')
        (out / f'{basename}_{table.replace("/", "_")}.tsv').write_bytes(data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw_dir', type=Path)
    args = parser.parse_args()
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    with (HERE / 'positions.csv').open() as f:
        positions = list(csv.DictReader(f))
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda pos: query(pos, args.raw_dir), positions))

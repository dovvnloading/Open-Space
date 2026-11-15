"""Rebuild compact input tables from the published Zenodo archive and VizieR."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import tarfile
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

DEST = Path(__file__).resolve().parent
RECORD = 'https://zenodo.org/api/records/4422053/files/TTVposteriors.tar.gz/content'
ARCHIVE_SHA = '9df25b6589d64d7cad07ef00cc4085293d18857680efe022a44e61db854c7af6'
CATALOG = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
IDS = {2092: ('2092.02', '2092.01', '2092.03'), 2086: ('2086.01', '2086.02', '2086.03')}


def export(archive, timing_dir=None):
    if hashlib.sha256(archive).hexdigest() != ARCHIVE_SHA:
        raise ValueError('Zenodo archive SHA-256 differs from the published snapshot')
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:gz') as data, (DEST/'posterior.csv').open('w') as fout:
        writer = csv.writer(fout, lineterminator='\n')
        writer.writerow(['koi', 'row_0based', *[f'{p}_{field}' for p in 'bcd' for field in ('mass_ratio_earth_solar', 'period_d', 'esinw', 'ecosw', 'transit_bjd_minus_2454900')]])
        for koi in IDS:
            member = data.extractfile(f'TTVposteriors/KOI-{koi:04d}.dat')
            rows = member.read().decode().splitlines()
            for i in np.random.default_rng(20261115).choice(len(rows), 128, replace=False):
                writer.writerow([koi, i, *rows[i].split()])
    with (DEST/'transits.csv').open('w') as f:
        writer=csv.writer(f, lineterminator='\n')
        writer.writerow(['koi', 'n', 'linear_bjd_minus_2454900', 'oc_min', 'error_min', 'out', 'over'])
        for koi in IDS:
            if timing_dir:
                lines=(timing_dir/f'ttv_{koi}.tsv').read_text().splitlines()
            else:
                query={'-source':'J/ApJS/225/9/table3','KOI':f'{koi}..{koi+0.99}','-out':'KOI,N,tn,O-C,e_O-C,Out,Over','-out.max':'unlimited'}
                lines=urlopen(CATALOG+'?'+urlencode(query),timeout=90).read().decode().splitlines()
            rows=[x.split() for x in lines if x.lstrip().startswith(f'{koi}.')]
            for row in sorted(rows,key=lambda r:(r[0],int(r[1]))):
                writer.writerow(row)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--timing-dir',type=Path)
    args=parser.parse_args()
    export(args.archive.read_bytes() if args.archive else urlopen(RECORD,timeout=180).read(),args.timing_dir)

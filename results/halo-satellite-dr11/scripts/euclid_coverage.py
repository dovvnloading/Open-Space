"""Intersect official Euclid DR1 Foundation MER MOCs with frozen ICRS rectangles."""
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import urllib.request

import fitsio
import hpgeom
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://nasa-irsa-euclid-dr1.s3.amazonaws.com/dr1/MOC/'

def main():
    fields = json.loads((ROOT/'catalog/footprint.json').read_text())['fields']
    rows = []
    for mode in ['wide_survey','wide_mode','deep_survey','deep_mode']:
        url = BASE + 'euclid_dr1_mer_catalogue_' + mode + '.moc.fits'
        with urllib.request.urlopen(url, timeout=90) as r:
            raw = r.read()
        with os.fdopen(os.memfd_create('moc'),'w+b') as f:
            f.write(raw); f.flush()
            with fitsio.FITS(f'/proc/self/fd/{f.fileno()}') as h:
                header = h[1].read_header()
                if header['ORDERING']!='NUNIQ' or header['COORDSYS']!='C':
                    raise ValueError('Unexpected MOC encoding or coordinate frame')
                uniq = h[1]['PIXEL'][:].astype(np.int64)
                max_order = int(header['MOCORDER'])
        for field, box in fields.items():
            # For each MOC order compare only its native cells to an inclusive
            # (conservative) rectangle cover. A count of zero proves disjointness.
            matches = 0
            for order in range(max_order+1):
                base = 4 * 4**order
                native = uniq[(uniq>=base)&(uniq<4*base)] - base
                if len(native):
                    rect = hpgeom.query_box(2**order,*box['ra'],*box['dec'],inclusive=True,nest=True)
                    matches += len(np.intersect1d(rect,native,assume_unique=True))
            rows.append(dict(field=field, mode=mode, intersecting_moc_cells=int(matches),
                             moc_order=max_order, input_sha256=hashlib.sha256(raw).hexdigest(),
                             input_size=len(raw), url=url, retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    (ROOT/'comparison/euclid_moc.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__ == '__main__':
    main()

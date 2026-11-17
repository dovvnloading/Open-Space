"""Verify the published receipt ledger against local byte streams without re-downloading them."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from counterparts import HERE


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);cache=ap.parse_args().cache
    tables=json.loads((HERE/'tables.json').read_text());rows=list(csv.DictReader((HERE/'photometry.csv').open()))
    catalog=list(csv.DictReader((HERE/'catalog.csv').open()))
    assert len(rows)==len(catalog)==2429
    assert {r['event'] for r in rows}=={r['event'] for r in catalog} and len({r['event'] for r in rows})==len(rows)
    for table in tables:
        assert sha(cache/table['year']/'lenses.par')==table['sha256']
    for row in rows:
        y,_,n=row['event'].split('-');path=cache/y/f'blg-{n}.dat'
        if row['status']!='ok' or sha(path)!=row['sha256'] or path.stat().st_size!=int(row['bytes']):
            raise ValueError(f"raw input mismatch: {row['event']}")
    ans=dict(verified_primary_tables=len(tables),verified_curves=len(rows),unique_bulge_ids=len(catalog),
             input_bytes=sum(int(r['bytes']) for r in rows),epochs=sum(int(r['n']) for r in rows),
             sha256_catalog_csv=sha(HERE/'catalog.csv'),sha256_photometry_csv=sha(HERE/'photometry.csv'),
             ancillary_2024=[])
    for event in ['2024-BLG-0519','2024-BLG-0680','2024-BLG-0825','2024-BLG-0944']:
        y,_,n=event.split('-');p=cache/y/f'blg-{n}.dat'
        if p.exists():
            ans['ancillary_2024'].append(dict(event=event,url=f'https://www.astrouw.edu.pl/ogle/ogle4/ews/{y}/blg-{n}/phot.dat',
                                              bytes=p.stat().st_size,sha256=sha(p)))
    print(json.dumps(ans,indent=2))

if __name__=='__main__':main()

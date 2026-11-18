"""Check public HST/JWST and Euclid Q2 image footprints at KMT lead positions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'other'))
from imaging import search, CACHE
from crossmatch import coords


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('events',nargs='+')
    p.add_argument('--out',type=Path,default=HERE/'archive-queries.json')
    a=p.parse_args()
    targets={}
    for cat in (HERE/'kmt').glob('catalog-*.tsv'):
        with cat.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                targets[row['Event']]=coords(row['RA'],row['Dec'])
    with (HERE/'ogle/catalog.csv').open() as f:
        for row in csv.DictReader(f):
            targets['OGLE-'+row['event']]=coords(row['ra_hms'],row['dec_dms'])
    (CACHE/'imaging').mkdir(parents=True,exist_ok=True)
    jobs=[(name,*targets[name],svc) for name in a.events for svc in ('MAST','IRSA')]
    with ThreadPoolExecutor(4) as pool: results=list(pool.map(search,jobs))
    a.out.write_text(json.dumps([{'query':q,'products':v} for q,v in results],indent=2)+'\n')


if __name__=='__main__':main()

"""Freeze SIMBAD cone identifications; associate only within 3 arcsec, retaining alternatives."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np

BASE = Path(__file__).resolve().parent
ENDPOINT = 'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
QUERY = """SELECT b.oid,b.main_id,b.ra,b.dec,b.otype,b.sp_type,b.rvz_redshift,
b.rvz_type,b.rvz_bibcode,b.nbref
FROM basic AS b
WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',270,66.56,0.52))
ORDER BY oid"""

def fetch(base,ra,dec):
    query=QUERY.replace('270,66.56',f'{ra},{dec}')
    url = ENDPOINT + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=10000, QUERY=query))
    with urlopen(url, timeout=120) as response:
        data = response.read()
        status = response.status
    rows = list(csv.DictReader(io.StringIO(data.decode('utf8'))))
    if not rows or len(rows) >= 10000 or 'rvz_type' not in rows[0]:
        raise ValueError('Empty, truncated or unexpected SIMBAD response')
    (base/'simbad.csv').write_bytes(data)
    receipt = dict(utc=dt.datetime.now(dt.timezone.utc).isoformat(),url=url,query=query,status=status,
                   rows=len(rows),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    (base/'simbad_query.json').write_text(json.dumps(receipt,indent=2)+'\n')

def read(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))

def associate(base):
    src = read(base.parent/'sources.csv')
    cat = read(base/'simbad.csv')
    ra,dec = np.radians([[float(o[k]) for o in cat] for k in ['ra','dec']])
    with (base/'matches.csv').open('w',newline='') as out:
        fields=['source_id','n_within_3arcsec','offset_arcsec','oid','main_id','otype','sp_type','rvz_redshift','rvz_type','rvz_bibcode','nbref']
        writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader()
        for row in src:
            a,d = np.radians([float(row['ra']),float(row['dec'])])
            h = np.sin((dec-d)/2)**2+np.cos(d)*np.cos(dec)*np.sin((ra-a)/2)**2
            dist = np.degrees(2*np.arcsin(np.sqrt(np.clip(h,0,1))))*3600
            inds=np.flatnonzero(dist<=3)
            for i in sorted(inds,key=lambda i:dist[i]):
                writer.writerow(dict(source_id=row['source_id'],n_within_3arcsec=len(inds),offset_arcsec=round(dist[i],4),
                                     **{k:cat[i][k] for k in fields[3:]}))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--fetch',action='store_true');p.add_argument('--match',action='store_true')
    p.add_argument('--wide',action='store_true')
    args=p.parse_args()
    base=BASE.parent/'wide'/'catalogs' if args.wide else BASE
    base.mkdir(parents=True,exist_ok=True)
    ra,dec=(209.53125,30) if args.wide else (270,66.56)
    if args.fetch:fetch(base,ra,dec)
    if args.match:associate(base)

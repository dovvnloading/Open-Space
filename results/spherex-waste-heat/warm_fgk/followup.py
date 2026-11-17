"""Check the published 100–600 K FGK fits for catalog counterparts and optical dips."""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
sys.path.insert(0, str(BASE / 'optical'))
import ztf
from fetch import position, separation

def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def output(name, keys, rows):
    with (ROOT / name).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def catalog(targets):
    cache = ROOT/'cache'
    near = []
    for t in targets:
        ra, dec = position(t, 2000)
        near.append(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.00416667))")
    query = 'SELECT main_id,otype,ra,dec FROM basic WHERE ' + ' OR '.join(near)
    failures = []
    for endpoint in ('https://simbad.cds.unistra.fr/simbad/sim-tap/sync',
                     'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'):
        url = endpoint + '?' + urlencode(dict(REQUEST='doQuery',LANG='ADQL', FORMAT='csv', MAXREC=1000,QUERY=query))
        try:
            with urlopen(url,timeout=80) as response:
                data = response.read()
            lines = csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
            if not {'ra','dec','otype'}.issubset(lines.fieldnames or []):
                raise ValueError(data[:180])
            rows = list(lines)
            if len(rows) >= 1000:
                raise ValueError('Result truncated')
            (cache/'simbad.csv.gz').write_bytes(gzip.compress(data,mtime=0))
            receipt = dict(url=url, query=query, rows=len(rows),sha256=hashlib.sha256(data).hexdigest(),
                           utc=dt.datetime.now(dt.timezone.utc).isoformat(),earlier_errors=failures)
            break
        except Exception as exc:
            failures.append(dict(url=url,error=str(exc)))
    else:
        receipt = dict(query=query,errors=failures)
        rows = []
    (ROOT/'simbad_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    matches = []
    for t in targets:
        ra, dec = position(t,2000)
        for s in rows:
            dist = separation(float(s['ra']),float(s['dec']),ra,dec)
            if dist <= 15:
                matches.append(dict(source_id=t['source_id'],main_id=s['main_id'],otype=s['otype'],
                                    separation_arcsec=round(dist,3),ra=s['ra'],dec=s['dec']))
    output('simbad_neighbors.csv',['source_id','main_id','otype','separation_arcsec','ra','dec'],matches)
    return receipt, matches

def main():
    fits = {r['source_id']:r for r in read(BASE/'literature'/'contardo_hogg.csv') if 100<=float(r['bb_temperature_k'])<=600}
    targets = [dict(row,label=row['source_id']) for row in read(BASE/'literature'/'fgk_astrometry.csv') if row['source_id'] in fits]
    assert len(targets)==16
    (ROOT/'cache').mkdir(exist_ok=True)
    ztf.ROOT = ROOT
    with ThreadPoolExecutor(5) as pool:
        cat = pool.submit(catalog,targets)
        optics = list(pool.map(ztf.inspect,targets))
        receipt, matches = cat.result()
    (ROOT/'ztf_receipts.json').write_text(json.dumps([r for r,_ in optics],indent=2)+'\n')
    keys='label oid points nights mjd_min mjd_max median_r p95_minus_median faintest_minus_median nights_fainter_by_0p2'.split()
    all_curves = [v for _,curves in optics for v in curves]
    output('r_lightcurves.csv',keys,all_curves)
    summaries = []
    for t,(q,curves) in zip(targets,optics):
        sid = t['source_id']
        identified = [m for m in matches if m['source_id']==sid]
        status = 'error' if 'error' in q else ('retained' if curves else 'no qualifying curve')
        summaries.append(dict(source_id=sid,T_published_K=fits[sid]['bb_temperature_k'], Lir_over_Lstar=fits[sid]['lir_over_lstar'],
                              distance_pc=round(1000/float(t['parallax']),1),ruwe=t['ruwe'],
                              simbad_status='query failed' if 'errors' in receipt else 'queried',
                              simbad_within_15arcsec=len(identified) if 'errors' not in receipt else '',
                              ztf_status=status,ztf_objects=len(curves),
                              ztf_objects_with_0p2mag_faint_night=sum(v['nights_fainter_by_0p2']>0 for v in curves)))
    output('summary.csv',summaries[0].keys(),summaries)
    with ZipFile(ROOT/'inputs.zip','w',compression=ZIP_DEFLATED) as archive:
        for item in sorted((ROOT/'cache').glob('*.gz')):
            archive.write(item,item.name)

if __name__=='__main__':
    main()

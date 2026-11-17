"""Inspect the archival RACS-high catalogue within 10 arcsec."""
import argparse
import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent


def query(item):
    oid, ra, dec = item
    query = ("SELECT name,ra,dec_corr,total_flux FROM AS110.racs_high_sources_v01 "
             "WHERE 1=CONTAINS(POINT('ICRS',ra,dec_corr),"
             f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},{10/3600:.9f}))")
    params = {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query}
    url = 'https://casda.csiro.au/casda_vo_tools/tap/sync?' + urlencode(params)
    receipt = [oid, datetime.now(UTC).isoformat(), url]
    try:
        with urlopen(url, timeout=45) as response:
            raw = response.read()
            if response.status != 200 or not raw.startswith(b'name,ra,dec_corr,total_flux'):
                raise ValueError('Unexpected CASDA TAP response')
        matches = list(csv.DictReader(raw.decode().splitlines()))
        return receipt + ['ok', hashlib.sha256(raw).hexdigest(), len(raw), len(matches), repr(matches)]
    except Exception as exc:
        return receipt + ['error', '', '', '', str(exc)]


def main(rubin, ztf):
    objects = {'control_RACS_0000+37_21':(4.14197827,36.19873653)}
    for path in (rubin, ztf):
        with path.open() as inp:
            for r in csv.DictReader(inp):
                objects.setdefault(r['oid'], (float(r.get('ra',r.get('meanra'))),float(r.get('dec',r.get('meandec')))))
    with ThreadPoolExecutor(4) as pool:
        result = list(pool.map(query, [(k, *v) for k,v in objects.items()]))
    with (HERE / 'archival_radio.csv').open('w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['oid','request_utc','url','status','sha256_response','bytes','n','matches_or_error'])
        w.writerows(result)
    if result[0][6] != 1:
        raise ValueError(f'RACS reference source not recovered: {result[0]}')
    print('successful', sum(r[3]=='ok' for r in result[1:]), 'science objects',len(result)-1)
    print('sources with archival counterparts', [r[0] for r in result[1:] if r[6]])


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('rubin',type=Path)
    p.add_argument('ztf',type=Path)
    a = p.parse_args()
    main(a.rubin, a.ztf)

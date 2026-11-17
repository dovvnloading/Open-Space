"""Query DESI DR1 for all modelled Q1 catalogue positions, then apply a 2 arcsec radial match.

python desi.py MODEL_TABLE_DIR RAW_CACHE OUT_DIR
"""
import argparse
import concurrent.futures
import hashlib
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from infer import read, write, sep


def collect(data, cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    cat = {r['id_str']: r for r in read(data / 'q1_discovery_engine_lens_catalog.csv')}
    ids = sorted({r['id_str'] for r in read(data / 'modeling_lens_mass.csv')} & set(cat))
    columns = 'targetid,mean_fiber_ra,mean_fiber_dec,z,zerr,zwarn,spectype,deltachi2,zcat_primary'
    endpoint = 'https://datalab.noirlab.edu/tap/sync?'
    def batch(pair):
        idx, keys = pair
        filename = cache / f'{idx:03}.csv'
        windows = []
        for k in keys:
            r = cat[k]
            ra, dec = float(r['right_ascension']), float(r['declination'])
            dd = 2.01 / 3600
            dr = dd / np.cos(np.deg2rad(dec))
            windows.append('(mean_fiber_ra BETWEEN %.8f AND %.8f AND mean_fiber_dec BETWEEN %.8f AND %.8f)'
                           % (ra-dr,ra+dr,dec-dd,dec+dd))
        adql = 'SELECT '+columns+' FROM desi_dr1.zpix WHERE '+' OR '.join(windows)+' ORDER BY targetid'
        params = dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=adql)
        if not filename.exists():
            for attempt in range(3):
                try:
                    response = urllib.request.urlopen(endpoint + urllib.parse.urlencode(params), timeout=120).read()
                    if not response.startswith(b'targetid,'):
                        raise ValueError(response[:320])
                    filename.write_bytes(response)
                    break
                except (OSError, ValueError):
                    if attempt == 2:
                        raise
                    time.sleep(2)
        rows = []
        for target in read(filename):
            for key in keys:
                r = cat[key]
                dist = sep(float(r['right_ascension']), float(r['declination']),
                           float(target['mean_fiber_ra']), float(target['mean_fiber_dec']))
                if dist <= 2:
                    adoptable = (dist <= 1 and target['zwarn'] == '0' and target['zcat_primary'] == 't'
                                 and target['spectype'] == 'GALAXY' and float(target['deltachi2']) >= 25
                                 and float(target['z']) > 0)
                    rows.append(dict(id_str=key, sep_arcsec=dist, **target, adoptable=int(adoptable)))
        return rows, dict(batch=idx, n_keys=len(keys), n_raw=len(read(filename)), n_matches=len(rows),
                          query_sha256=hashlib.sha256(adql.encode()).hexdigest(),
                          response_sha256=hashlib.sha256(filename.read_bytes()).hexdigest())
    chunks = [(i,ids[j:j+20]) for i,j in enumerate(range(0,len(ids),20))]
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        result = list(pool.map(batch,chunks))
    write(out/'desi_matches.csv', [r for rows,_ in result for r in rows])
    write(out/'desi_manifest.csv', [m for _,m in result])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for field in ['model_dir','cache','out']:
        p.add_argument(field,type=Path)
    a = p.parse_args()
    collect(a.model_dir, a.cache, a.out)

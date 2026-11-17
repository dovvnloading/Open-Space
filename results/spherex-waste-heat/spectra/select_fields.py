"""Retrieve Gaia astrometry and select comparison stars inside each SPHEREx tile."""
import argparse
import csv
import hashlib
import io
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / 'neowise-nuclear-flares'))
from spherex_tools import healpix_nest

ENDPOINT = 'https://gea.esac.esa.int/tap-server/tap/sync'
COLUMNS = ('source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,parallax_over_error,'
           'phot_g_mean_mag,bp_rp,ruwe')
RECEIPTS = []

def gaia(query):
    url = ENDPOINT + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=query))
    with urlopen(url, timeout=120) as response:
        raw = response.read()
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    if b'source_id,' not in raw[:100]:
        raise RuntimeError(raw[:300])
    RECEIPTS.append(dict(url=url, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(), rows=len(rows)))
    return rows

def pos(row):
    # Year 1 representative epoch. Annual parallax here is subpixel (few mas).
    years = 2025.9 - float(row['ref_epoch'])
    ra, dec = float(row['ra']), float(row['dec'])
    return ((ra + float(row['pmra'])*years/3.6e6/math.cos(math.radians(dec))) % 360,
            dec + float(row['pmdec'])*years/3.6e6)

def field(row):
    ra,dec = pos(row)
    tile = healpix_nest(5,ra,dec)
    g, c = map(float,[row['phot_g_mean_mag'],row['bp_rp']])
    minimum_parallax = .15 if row['label'].startswith('FGK') else 3.33
    query = (f'SELECT {COLUMNS} FROM gaiadr3.gaia_source WHERE '
             f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.7f},{dec:.7f},0.85)) "
             f'AND parallax_over_error>10 AND ruwe<1.4 AND parallax>{minimum_parallax:.4f} '
             f'AND bp_rp BETWEEN {c-.5:.5f} AND {c+.5:.5f} '
             f'AND phot_g_mean_mag BETWEEN {g-1.5:.5f} AND {g+1.5:.5f}')
    mg = g + 5*math.log10(float(row['parallax'])) - 10
    pool = []
    for other in gaia(query):
        if other['source_id'] == row['source_id'] or healpix_nest(5,*pos(other)) != tile:
            continue
        og,oc = map(float,[other['phot_g_mean_mag'],other['bp_rp']])
        omg = og + 5*math.log10(float(other['parallax'])) - 10
        if abs(omg-mg) > 1: continue
        score=((og-g)/1.0)**2+((oc-c)/.2)**2+((omg-mg)/.6)**2
        pool.append((score,other['source_id'],other))
    ans=[]
    for idx, item in enumerate([row]+[p[2] for p in sorted(pool)[:2]]):
        nra,ndec = pos(item)
        ans.append(dict(label=row['label']+(f'_cmp{idx}' if idx else ''), **{k:item[k] for k in COLUMNS.split(',')},
                        ra_2025_9=nra,dec_2025_9=ndec,tile=tile))
    return ans

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fgk',action='store_true',help='First three Contardo & Hogg Table 2 sources plus all entries within 300 pc')
    opts=parser.parse_args()
    suffix='_fgk' if opts.fgk else ''
    labels={r['source_id']:r['label'] for r in csv.DictReader((ROOT.parent/'literature/candidates.csv').open())}
    if opts.fgk:
        fgk=list(csv.DictReader((ROOT.parent/'literature/fgk_astrometry.csv').open()))
        subset=fgk[:3] + [r for r in fgk[3:] if float(r['parallax'])>1000/300]
        labels={row['source_id']:f'FGK{index+1}' for index,row in enumerate(subset)}
    rows=gaia(f"SELECT {COLUMNS} FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(labels)})")
    if set(labels)!={r['source_id'] for r in rows}: raise ValueError('incomplete Gaia ID lookup')
    for row in rows: row['label']=labels[row['source_id']]
    with ThreadPoolExecutor(max_workers=4) as pool: fields=sum(pool.map(field,rows),[])
    with (ROOT/f'positions{suffix}.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields[0]); w.writeheader(); w.writerows(sorted(fields,key=lambda r:r['label']))
    (ROOT/f'gaia_queries{suffix}.json').write_text(json.dumps(sorted(RECEIPTS,key=lambda r:r['url']),indent=2)+'\n')
if __name__=='__main__': main()

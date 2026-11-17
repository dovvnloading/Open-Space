"""Extract planet-centered calibrated detector regions in complementary NIRCam channels."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from jwst_cutouts import make_cut


def write(path, fields, rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(fields);w.writerows(rows)


def selected(r):
    name=r['productFilename']
    if not name.endswith('_cal.fits') or '_nrc' not in name:return False
    if name.split('_')[-2].endswith('long'):return True
    if r['proposal_id']=='6379' and r['filters'] in ('F070W','F115W'):return True
    return r['proposal_id']=='8975' and int(name[7:10]) in (2,4,6,8)


def run(cache, uranus_f210m=False):
    cache.mkdir(exist_ok=True,parents=True)
    rows=list(csv.DictReader((ROOT.parent/'manifests/products.csv').open()))
    if uranus_f210m:
        rows=[r for r in rows if r['planet']=='Uranus' and r['proposal_id']=='2739' and r['filters']=='F210M' and r['productFilename'].endswith('_cal.fits')]
    else:
        rows=[r for r in rows if selected(r)]
    prefix='u2739_' if uranus_f210m else ''
    def get(r):
        base,spans=make_cut(r,166 if 'long' in r['productFilename'] else 232,cache)
        return (base,spans)
    with ThreadPoolExecutor(12) as pool:result=list(pool.map(get,rows))
    write(ROOT/(prefix+'cutouts.csv'),'source filter visit midpoint_mjd_utc planet_x planet_y cutout_x0 cutout_y0 filename output_bytes output_sha256'.split(),[r[0] for r in result])
    write(ROOT/(prefix+'ranges.csv'),'source offset length sha256'.split(),[s for r in result for s in r[1]])
    print(len(rows),sum(bool(r[0][8]) for r in result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path)
    p.add_argument('--uranus-f210m',action='store_true')
    a=p.parse_args();run(a.cache,a.uranus_f210m)

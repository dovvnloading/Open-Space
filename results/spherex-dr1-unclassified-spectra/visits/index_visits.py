"""Index held-out QR3 detector footprints around the frozen field; cache TAP replies."""
import argparse
import csv
import hashlib
import io
import json
import struct
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

TAP = 'https://irsa.ipac.caltech.edu/TAP/sync'
IRSA = 'https://irsa.ipac.caltech.edu/'
FIELD = (270., 66.56)
HERE = Path(__file__).resolve().parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def cut_url(uri, ra, dec):
    return IRSA + uri + '?' + urlencode({'center': f'{ra:.8f},{dec:.8f}d', 'size': '0.065'})

def corners(poly):
    """Four DSS/IRSA detector corners, ordered (low x low y, high low, high high, low high)."""
    wkb = bytes.fromhex(poly)
    if len(wkb) != 97 or struct.unpack_from('<III', wkb, 5) != (4326, 1, 5):
        raise ValueError('unexpected archive footprint WKB')
    return np.asarray(struct.unpack_from('<10d', wkb, 17)).reshape(5, 2)[:4]

def projected(points, ra, dec):
    """Gnomonic coordinates around target, omitting the unused angular scale."""
    a, d = np.deg2rad(points.T)
    da = a - np.deg2rad(ra)
    d0 = np.deg2rad(dec)
    cosc = np.sin(d) * np.sin(d0) + np.cos(d) * np.cos(d0) * np.cos(da)
    return np.stack((np.cos(d) * np.sin(da) / cosc,
                     (np.sin(d)*np.cos(d0) - np.cos(d)*np.sin(d0)*np.cos(da)) / cosc), axis=-1)

def approximate_xy(poly, ra, dec):
    points = projected(corners(poly), ra, dec)
    # The inverse of an affine corner approximation is adequate for choosing images.
    # Final positions and wavelengths are always recomputed using the exposure WCS.
    xaxis = ((points[1]-points[0])+(points[2]-points[3]))/2
    yaxis = ((points[3]-points[0])+(points[2]-points[1]))/2
    return (np.linalg.solve(np.stack((xaxis,yaxis),axis=-1), -points.mean(axis=0))+.5)*2040-.5

def read_url(url, timeout=150):
    error = None
    for _ in range(3):
        try:
            with urlopen(url,timeout=timeout) as resp:
                return resp.read()
        except OSError as exc:
            error = exc
    raise error

def build(detector, cache, out, ra, dec):
    query = f"""SELECT a.uri,p.time_bounds_lower,p.poly
FROM spherex.artifact a JOIN spherex.plane p ON a.planeid=p.planeid
WHERE 1=INTERSECTS(p.poly,CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.5))
AND p.energy_bandpassname='SPHEREx-D{detector}'
ORDER BY p.time_bounds_lower,a.uri"""
    url = TAP+'?'+urlencode(dict(QUERY=query, LANG='ADQL', FORMAT='csv', REQUEST='doQuery', MAXREC=100000))
    dest = cache/f'footprints-D{detector}.csv'
    if not dest.exists() or not dest.read_bytes().startswith(b'uri,'):
        dest.write_bytes(read_url(url))
    data = dest.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if not rows or 'uri' not in rows[0] or len(rows)>=100000:
        raise ValueError(f'invalid or capped TAP response: {data[:300]}')
    qr3_count=sum('/qr3/' in r['uri'] for r in rows)
    rows = sorted((r for r in rows if '/qr3/' in r['uri'] or '/dr1/' in r['uri']), key=lambda r:(float(r['time_bounds_lower']),r['uri']))
    for row in rows:
        x,y=approximate_xy(row['poly'],ra,dec)
        if 300<x<1750 and 300<y<1750:
            break
    else:
        raise ValueError('cannot find a central WCS reference')
    url_cut = cut_url(row['uri'],ra,dec)
    raw=read_url(url_cut)
    with fits.open(io.BytesIO(raw)) as hd:
        grid = hd['WCS-WAVE'].data[0]
        result=dict(detector=detector, reference_cutout=url_cut,reference_sha256=sha(raw),
                    x=grid['X'].tolist(), y=grid['Y'].tolist(), wave_um=grid['VALUES'][:,:,0].tolist(),
                    fwhm_um=grid['VALUES'][:,:,1].tolist())
    (out/f'calibration_D{detector}.json').write_text(json.dumps(result,sort_keys=True,indent=1)+'\n')
    (out/f'index_D{detector}.json').write_text(json.dumps(dict(query=query,endpoint=TAP,sha256=sha(data),
                    qr3_rows=qr3_count,n_rows=len(data.splitlines())-1,first_mjd=rows[0]['time_bounds_lower'],last_mjd=rows[-1]['time_bounds_lower']),
                    indent=2)+'\n')
    print(detector,len(rows),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache',type=Path)
    parser.add_argument('out',type=Path)
    parser.add_argument('ra',type=float)
    parser.add_argument('dec',type=float)
    parser.add_argument('--detectors',type=int,nargs='+',default=list(range(1,7)))
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.out.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda d:build(d,args.cache,args.out,args.ra,args.dec),args.detectors))

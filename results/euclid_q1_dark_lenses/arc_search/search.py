"""Search fixed strips of native Euclid Q1 VIS data. Python 3.10+, NumPy, SciPy."""
import argparse
import hashlib
import json
import math
import os
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
RADII = (4, 6, 9, 12)  # 0.2 arcsec binned pixels
MARGIN = 24
NATIVE_CHUNK = 1200
N_AZ = 8
SECTOR_MIN = 2.2
MEAN_MIN = 3.4
CENTER_MAX = 2.5


def digest(data):
    return hashlib.sha256(data).hexdigest()


def get_range(url, a, b):
    for attempt in range(3):
        try:
            req = Request(f'{url}?range_start={a}&range_end={b}',
                          headers={'Range': f'bytes={a}-{b}', 'Accept-Encoding':'identity'})
            with urlopen(req, timeout=120) as f:
                if f.status != 206 or f.headers['Content-Range'].split('/')[0] != f'bytes {a}-{b}':
                    raise ValueError('archive did not return requested range')
                data = f.read(b-a+2)
                if len(data) != b-a+1:
                    raise ValueError('incomplete range')
                return data
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def decode_header(raw):
    cards = [raw[k:k+80].decode('ascii') for k in range(0,len(raw),80)]
    end = next(j for j, c in enumerate(cards) if c[:8].strip() == 'END')
    h = {c[:8].strip():c[10:].split('/')[0].strip().strip("'").strip()
         for c in cards[:end] if c[8:10] == '= '}
    assert h['BITPIX'] == '-32' and h['NAXIS1'] == h['NAXIS2'] == '19200'
    assert h['CTYPE1'] == 'RA---TAN' and h['CTYPE2'] == 'DEC--TAN'
    assert h['FILTER'] == 'VIS' and h['DATASETR'].startswith('Q1')
    return h, ((end+1)*80+2879)//2880*2880


def world(h, x, y):
    dx, dy = x+1-float(h['CRPIX1']), y+1-float(h['CRPIX2'])
    xi, eta = np.deg2rad([float(h['CD1_1'])*dx+float(h['CD1_2'])*dy,
                           float(h['CD2_1'])*dx+float(h['CD2_2'])*dy])
    d0 = math.radians(float(h['CRVAL2']))
    den = math.cos(d0)-eta*math.sin(d0)
    ra = (float(h['CRVAL1']) + math.degrees(math.atan2(xi,den))) % 360
    dec = math.degrees(math.atan2(math.sin(d0)+eta*math.cos(d0), math.hypot(den,xi)))
    return round(ra,8), round(dec,8)


def sigma(values):
    med = np.median(values)
    return float(med), float(1.4826*np.median(np.abs(values-med)))


def find_rings(arr, good):
    # Kernels are measured against the actual mosaic noise, including correlations.
    # The margin contains the full exterior reference ring and the Gaussian support.
    eligible = ndi.minimum_filter(good.astype('u1'), size=2*MARGIN+1, mode='constant') > 0
    if eligible.sum() < 300:
        return [], 0
    sm = ndi.gaussian_filter(arr, 1., mode='nearest', truncate=3)
    med, sc = sigma(sm[eligible])
    if sc <= 0:
        return [], 0
    center = (sm-med)/sc
    hits = []
    for radius in RADII:
        sectors = []
        for k in range(N_AZ):
            # Sample three azimuths per sector; compare to the local exterior sky
            # in the same direction to reduce false responses next to gradients.
            sector = np.zeros_like(sm)
            for delta in (-0.28, 0., 0.28):
                ang = (k+delta)*2*np.pi/N_AZ
                def offset(r):
                    return (-round(r*np.sin(ang)), -round(r*np.cos(ang)))
                sector += np.roll(sm,offset(radius),(0,1)) - np.roll(sm,offset(1.65*radius),(0,1))
            base, noise = sigma(sector[eligible][::4])
            if noise <= 0:
                break
            sectors.append((sector-base)/noise)
        if len(sectors) != N_AZ:
            continue
        sectors = np.stack(sectors)
        mean = sectors.mean(axis=0)
        nsector = (sectors >= SECTOR_MIN).sum(axis=0)
        ok = eligible & (center <= CENTER_MAX) & (nsector >= 6) & (mean >= MEAN_MIN)
        maxima = ok & (mean == ndi.maximum_filter(np.where(ok,mean,-np.inf),size=2*radius+1))
        for y,x in zip(*np.nonzero(maxima)):
            hits.append((float(mean[y,x]),int(x),int(y),int(radius),int(nsector[y,x]),float(center[y,x])))
    keep = []
    for point in sorted(hits,reverse=True):
        if all(math.hypot(point[1]-other[1],point[2]-other[2]) > max(point[3],other[3]) for other in keep):
            keep.append(point)
    return keep, int(eligible.sum())


def run_strip(item, y0, height, h, headsize):
    width = int(h['NAXIS1'])
    start = headsize + 4*width*y0
    # For transport each 384-row strip is split into 8 requests of 3.69 MB.
    step = 48*width*4
    end = start + height*width*4
    with ThreadPoolExecutor(4) as pool:
        parts = list(pool.map(lambda a: get_range(item['url'], a, min(a+step,end)-1), range(start,end,step)))
    raw = b''.join(parts)
    del parts
    im = np.frombuffer(raw,'>f4').reshape(height,width)
    report = {'tile':item['tile'], 'row':y0,'height':height,'first_byte':start,'last_byte':end-1,
              'sha256':digest(raw),'eligible_binned_centers':0,'triggers':[],'injections':[]}
    for x0 in range(0,width,NATIVE_CHUNK):
        a = np.array(im[:,x0:x0+NATIVE_CHUNK],dtype='f4')
        valid = np.isfinite(a) & (a != 0) & (np.abs(a)<1e5)
        # Rare native zero-valued pixels are conservatively masked as potential no-data.
        a[~valid] = 0.
        a = a.reshape(height//2,2,NATIVE_CHUNK//2,2).mean(axis=(1,3))
        valid = valid.reshape(height//2,2,NATIVE_CHUNK//2,2).all(axis=(1,3))
        found, n = find_rings(a,valid)
        report['eligible_binned_centers'] += n
        for avg,x,y,r,nsector,c in found:
            xp,yp = x0+2*x+0.5,y0+2*y+0.5
            ra,dec = world(h,xp,yp)
            report['triggers'].append({'x0':xp,'y0':yp,'ra':ra,'dec':dec,'radius_arcsec':r*.2,
                                      'sector_mean':round(avg,3),'sectors':nsector,'center_z':round(c,3)})
        # One fixed injection experiment per strip, distinct from the sky triggers.
        if x0 == 8*NATIVE_CHUNK and n:
            for j,(r,amplitude) in enumerate([(r,s) for r in RADII for s in (2,5)]):
                cx,cy = 75+(j%4)*150,48+(j//4)*96
                if not valid[cy-MARGIN:cy+MARGIN+1,cx-MARGIN:cx+MARGIN+1].all():
                    continue
                sy,sx = slice(cy-32,cy+33),slice(cx-32,cx+33)
                _,noise = sigma(a[valid])
                yy,xx = np.ogrid[-32:33,-32:33]
                added = amplitude*noise*np.exp(-.5*((np.hypot(xx,yy)-r)/1.0)**2)
                ai = a.copy()
                ai[sy,sx] += added
                # A pre-existing matched detection cannot count as an injection recovery.
                def matched(peaks):
                    return any(np.hypot(f[1]-cx,f[2]-cy)<=2 and f[3]==r for f in peaks)
                fi,_ = find_rings(ai,valid)
                report['injections'].append({'r':r*.2,'peak_sigma':amplitude,'preexisting':matched(found),'recovered':matched(fi)})
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=HERE/'runs')
    p.add_argument('--tile')
    p.add_argument('--rows',type=int,help='limit rows per field for a transport smoke test')
    opts=p.parse_args()
    opts.out.mkdir(exist_ok=True,parents=True)
    for item in json.loads((HERE/'fields.json').read_text()):
        if opts.tile and item['tile'] != opts.tile:
            continue
        raw=get_range(item['url'],0,5759)
        h,headsize=decode_header(raw)
        for y0 in item['rows_zero_based'][:opts.rows]:
            path=opts.out/f"{item['tile']}_{y0}.json"
            if path.exists():
                continue
            try:
                report=run_strip(item,y0,item['height'],h,headsize)
                report['header_sha256']=digest(raw)
                report['script_sha256']=digest(Path(__file__).read_bytes())
                tmp=path.with_suffix('.tmp')
                tmp.write_text(json.dumps(report,indent=2)+'\n')
                os.replace(tmp,path)
                print(item['tile'],y0,report['eligible_binned_centers'],len(report['triggers']),flush=True)
            except (OSError,ValueError) as exc:
                print('FAILED',item['tile'],y0,str(exc),flush=True)
if __name__=='__main__':
    main()

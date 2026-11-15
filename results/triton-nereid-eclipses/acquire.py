"""Build the local event-state table from JPL Horizons and the cited SPK files."""
import csv
import io
import struct
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from functools import lru_cache
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from numpy.polynomial.chebyshev import chebval, chebder

DAY = 86400.0
J2000 = 2451545.0
JPL = 'https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/'
SAI = 'https://ftp.imcce.fr/pub/ephem/satel/SAI/NEPTUNE/'
KERNELS = {
    'nep097': JPL + 'nep097.bsp',
    'nep105': JPL + 'nep105.bsp',
    'sai801': SAI + '2015/SAI-2015-Triton-v2.bsp',
    'sai802': SAI + '2011/SAI-2011-Nereid.bsp',
}
EVENTS = [('2024', '2024-02-07 12:00'), ('2026', '2026-01-30 19:00')]
BODIES = [8, 10, 399, 801, 802, 899]


def horizons(body, start):
    stop = (datetime.fromisoformat(start) + timedelta(hours=10)).isoformat(' ', timespec='minutes')
    params = dict(format='text', COMMAND=str(body), CENTER='500@0', EPHEM_TYPE='VECTORS',
                  START_TIME=start, STOP_TIME=stop, STEP_SIZE='30m', OUT_UNITS='KM-S',
                  REF_PLANE='FRAME', TIME_TYPE='TDB', VEC_TABLE='2', VEC_CORR='NONE',
                  CSV_FORMAT='YES', OBJ_DATA='NO')
    params = {k: v if k == 'format' else "'"+v+"'" for k,v in params.items()}
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urlencode(params)
    raw = urlopen(url, timeout=90).read().decode()
    lines = raw.split('$$SOE')[1].split('$$EOE')[0]
    return np.loadtxt(io.StringIO(lines), delimiter=',', usecols=(0,2,3,4,5,6,7))


@lru_cache(maxsize=256)
def fetch_range(url, start, size):
    end = start + size - 1
    req = Request(url+f'?window={start}-{end}', headers={'Range': f'bytes={start}-{end}'})
    with urlopen(req, timeout=65) as rsp:
        if rsp.status != 206:
            raise IOError(f'Expected partial response, received {rsp.status}')
        data = rsp.read()
    if len(data) != size:
        raise IOError(f'Expected {size} bytes, received {len(data)}')
    return data


def segment(url, target, et):
    header = fetch_range(url, 0, 1024)
    endian = {'LTL-IEEE':'<', 'BIG-IEEE':'>'}[header[88:96].decode().strip()]
    nd, ni = struct.unpack_from(endian+'ii', header, 8)
    assert (nd, ni) == (2,6)
    record = struct.unpack_from(endian+'i', header, 76)[0]
    found = None
    while record:
        summaries = fetch_range(url, (record-1)*1024, 1024)
        nxt, _, count = struct.unpack_from(endian+'3d', summaries)
        for i in range(int(count)):
            lo, hi, dest, origin, frame, typ, first, last = struct.unpack_from(endian+'2d6i', summaries, 24+i*40)
            if dest == target and lo <= et <= hi:
                found = (origin, frame, typ, first, last)
        record = int(nxt)
    if found is None:
        raise ValueError((url,target,et))
    return endian, found


def states(kernel, target, times, base):
    url = KERNELS[kernel]
    etimes = (times-J2000)*DAY
    endian, (origin, frame, kind, first, last) = segment(url, target, etimes.mean())
    if frame != 1 or kind not in (2,3):
        raise ValueError((url, frame, kind))
    # SPK uses one-based double-word addresses. Read only the needed coefficient records.
    init, interval, rec_size, _ = struct.unpack(endian+'4d', fetch_range(url, (last-4)*8, 32))
    rec_size = int(rec_size)
    ids = np.floor((etimes-init)/interval).astype(int)
    begin, end = ids.min(), ids.max()
    raw = fetch_range(url, (first-1+rec_size*begin)*8, rec_size*(end-begin+1)*8)
    recs = np.frombuffer(raw, dtype=endian+'f8').reshape(-1,rec_size)
    values = []
    for epoch, idx in zip(etimes, ids):
        mid, radius, *coefs = recs[idx-begin]
        c = np.reshape(coefs, (3 if kind==2 else 6,-1))
        z = (epoch-mid)/radius
        pos = [chebval(z, k) for k in c[:3]]
        vel = ([chebval(z,chebder(k))/radius for k in c[:3]] if kind==2
               else [chebval(z,k) for k in c[3:]])
        values.append(pos+vel)
    return np.column_stack([times, base[origin][:,1:]+np.asarray(values)])


def build(output, include_sai=True):
    with open(output,'w',newline='') as f:
        w = csv.writer(f,lineterminator='\n')
        w.writerow(('event','solution','body','jd_tdb','x_km','y_km','z_km','vx_km_s','vy_km_s','vz_km_s'))
        for label,start in EVENTS:
            with ThreadPoolExecutor(max_workers=2) as pool:
                base = dict(zip(BODIES,pool.map(lambda b: horizons(b,start),BODIES)))
            tables = {('horizons',b):a for b,a in base.items()}
            cases = [('nep097',801),('nep097',802),('nep105',802)]
            if include_sai:
                cases.extend([('sai801',801),('sai802',802)])
            for kernel,target in cases:
                tab = states(kernel,target,base[target][:,0],base)
                tables[kernel,target] = tab
            for (name,body),a in tables.items():
                for row in a:
                    w.writerow([label,name,body,*[format(v,'.16g') for v in row]])
            f.flush()


if __name__ == '__main__':
    build(sys.argv[1] if len(sys.argv)>1 else 'states.csv')

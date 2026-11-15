"""Re-query the 174 frozen encounters. Cache original responses outside the repository."""
import argparse
import csv
import hashlib
import json
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from measure import DAY, C, ROOT, TT_UTC, angles, covered, jd

HORIZONS = 'https://ssd.jpl.nasa.gov/api/horizons.api'
MPC = 'https://data.minorplanetcenter.net/api/get-obs'
CAD = 'https://ssd-api.jpl.nasa.gov/cad.api?' + urlencode({
    'date-min': '2018-01-01', 'date-max': '2050-12-31', 'dist-max': '0.002', 'fullname': 'true', 'sort': 'dist'})


def get(cache, key, url, body=None):
    path = cache / (key + '.json')
    if not path.exists():
        req = Request(url, data=json.dumps(body).encode() if body else None, method='GET')
        if body:
            req.add_header('Content-Type', 'application/json')
        for attempt in range(3):
            try:
                with urlopen(req, timeout=90) as response:
                    raw = response.read()
                json.loads(raw)
                path.write_bytes(raw)
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(2)
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def horizons_url(parameters):
    return HORIZONS + '?' + urlencode(dict(format='json', **{k: f"'{v}'" for k, v in parameters.items()}))


def vector_api(cache, key, name, start, end, step):
    url = horizons_url(dict(COMMAND=name, CENTER='500@399', EPHEM_TYPE='VECTORS', START_TIME=start,
                       STOP_TIME=end, STEP_SIZE=step, REF_PLANE='FRAME', VEC_TABLE='2', VEC_CORR='NONE',
                       OUT_UNITS='KM-S', CSV_FORMAT='YES', OBJ_DATA='NO'))
    data, digest = get(cache, key, url)
    table = data['result'].split('$$SOE')[1].split('$$EOE')[0]
    rows = [line.split(',') for line in table.strip().splitlines()]
    t = np.array([float(r[0]) for r in rows])
    rv = np.array([[float(v) for v in r[2:8]] for r in rows])
    return CubicHermiteSpline(t, rv[:, :3], rv[:, 3:] * DAY), digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--only', help='One frozen target, e.g. "2024 XS16"')
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    sun, _ = vector_api(args.cache, 'sun', '10', '2017-12-01', '2051-02-01', '5 d')
    with (ROOT / 'encounters.csv').open() as handle:
        cases = list(csv.DictReader(handle))
    for case in cases:
        n = case['target']
        if args.only and n != args.only:
            continue
        t = float(case['ca_jd_tdb'])
        name = (n if n.isdecimal() else 'DES=' + n) + ';'
        key = n.replace(' ', '_') + '_' + str(int(t))
        body, digest = vector_api(args.cache, key, name, f'JD{t-1:.8f}', f'JD{t+1:.8f}', '5 m')
        ts = np.linspace(t-1, t+1, 17281)
        sep, e, s = angles(body, sun, ts, 6371.)
        # At these distances Earth is larger than the Sun near deep eclipse.
        # Rank with normalized signed limb distance; only deeply eclipsing
        # targets advance to the MPC exposure-time screen.
        i = np.argmin((sep-e)/s)
        frac = covered(sep[i], e[i], s[i])
        hits = 0
        if frac > .05 and case['future_at_search'] == '0':
            obj, _ = get(args.cache, 'mpc_' + n.replace(' ', '_'), MPC, {'desigs':[n], 'output_format':['ADES_DF']})
            for obs in obj[0]['ADES_DF']:
                if obs['Obstype'] != 'optical' or not obs['mag'] or not obs['obstime']:
                    continue
                rec = jd(datetime.fromisoformat(obs['obstime'].replace('Z', '+00:00'))) + TT_UTC/DAY
                if abs(rec-t) > 1:
                    continue
                emit = rec - np.linalg.norm(body(rec))/C/DAY
                hits += covered(*angles(body, sun, emit, 6371.)) > .05
        print(n, f'{frac:.6f}', hits, digest)
    # Horizons response headers include a request timestamp. Full live-response
    # hashes can differ even when every element of the returned state agrees.


if __name__ == '__main__':
    main()

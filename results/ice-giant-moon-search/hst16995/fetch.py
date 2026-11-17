"""Retrieve frozen HST FLTs and spacecraft-relative Horizons ephemerides."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
SPECIAL = {75051: '2023U1', 85051: 'S2002_N5', 85052: 'S2021_N1'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(cache, export):
    cache.mkdir(exist_ok=True, parents=True)
    selected = list(csv.DictReader((ROOT/'selection.csv').open()))
    hashes, coords, receipts = [], [], []
    for r in selected:
        name = r['rootname'] + '_flt.fits'
        path = cache / name
        url = f"https://stpubdata.s3.amazonaws.com/hst/public/{name[:4]}/{r['rootname']}/{name}"
        if not path.exists():
            with urlopen(url, timeout=120) as response:
                raw = response.read()
            if not raw.startswith(b'SIMPLE  =') or len(raw) % 2880:
                raise ValueError(f'Incomplete FITS: {name}')
            path.write_bytes(raw)
        raw = path.read_bytes()
        hashes.append(dict(file=name, bytes=len(raw), sha256=digest(raw), url=url))

    for target, codes in [('NEPTUNE', [899, *range(801, 815), 85051, 85052]),
                          ('URANUS', [799, *range(701, 728), 75051])]:
        rows = [r for r in selected if r['target'] == target]
        jd = [(float(r['mjd_start'])+float(r['mjd_end']))/2 + 2400000.5 for r in rows]
        times = '\n'.join('%.9f' % t for t in jd)
        for code in codes:
            params = dict(format='json', COMMAND=f"'{code}'", OBJ_DATA='NO', MAKE_EPHEM='YES',
                          EPHEM_TYPE='OBSERVER', CENTER="'500@-48'", TLIST="'"+times+"'",
                          TLIST_TYPE='JD', QUANTITIES="'1,19,20,21,24'", CSV_FORMAT='YES',
                          ANG_FORMAT='DEG', EXTRA_PREC='YES', TIME_DIGITS='FRACSEC')
            url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urlencode(params)
            path = cache / f'{code}.json'
            if not path.exists():
                with urlopen(url, timeout=120) as response:
                    path.write_bytes(response.read())
            raw = path.read_bytes()
            result = json.loads(raw)['result']
            ident = next(line for line in result.splitlines() if line.startswith('Target body name:'))
            if f'({code})' not in ident or SPECIAL.get(code, '') not in ident or '_merged' not in ident:
                raise ValueError(f'Wrong target for {code}: {ident}')
            table = result.split('$$SOE')[1].split('$$EOE')[0].strip()
            from astropy.time import Time
            from datetime import datetime, timezone
            for i, (r, values) in enumerate(zip(rows, csv.reader(table.splitlines()), strict=True)):
                timestamp = datetime.strptime(values[0].strip(), '%Y-%b-%d %H:%M:%S.%f').replace(tzinfo=timezone.utc)
                if abs(Time(timestamp).jd - jd[i])*86400 > .003:
                    raise ValueError(f'Epoch mismatch: {r["rootname"]}, {code}')
                coords.append(dict(rootname=r['rootname'], body=code, ra=values[3].strip(),
                                   dec=values[4].strip(), r_au=values[5].strip(), delta_au=values[7].strip(),
                                   light_min=values[9].strip(), phase_deg=values[10].strip()))
            receipts.append(dict(body=code, identity=ident.strip(), sha256=digest(raw), url=url))
    for name, rows in [('inputs.csv', hashes), ('ephemerides.csv', coords), ('queries.csv', receipts)]:
        path = ROOT/name
        if export:
            with path.open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        else:
            frozen = list(csv.DictReader(path.open()))
            if frozen != [{k: str(v) for k,v in r.items()} for r in rows]:
                raise ValueError(f'{name} changed; inspect archive recalibration or ephemeris revision')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--freeze', action='store_true', help='write the receipts for a new run')
    args = parser.parse_args()
    run(args.cache, args.freeze)

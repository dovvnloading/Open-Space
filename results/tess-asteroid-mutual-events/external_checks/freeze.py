"""Build the small public-data snapshot consumed by check.py."""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MPC = 'https://data.minorplanetcenter.net/api/get-obs'
HORIZONS = 'https://ssd.jpl.nasa.gov/api/horizons.api'
TSSYS = 'https://archive.konkoly.hu/pub/tssys/dr1/lightcurves_spectra/'


def julian_date(field):
    year, month, day = field.split()
    date = datetime.date(int(year), int(month), int(float(day)))
    return date.toordinal() + 1721424.5 + float(day) % 1


def get(request, cache):
    return cache.read_bytes() if cache else urllib.request.urlopen(request, timeout=90).read()


def freeze(cache, local_tess):
    manifest, photometry, spin_rows, geometry = [], [], [], []
    for number in (29489, 3177):
        payload = {'desigs': [str(number)], 'output_format': ['OBS80']}
        request = urllib.request.Request(MPC, data=json.dumps(payload).encode(),
                    headers={'Content-Type': 'application/json'}, method='GET')
        inputs = {}
        data = get(request, cache/f'{number}_obs.json' if cache else None)
        inputs['mpc'] = data
        for line in json.loads(data)[0]['OBS80'].splitlines():
            if len(line) < 80 or line[14] in ('s', 'v', 'r', 'R'):
                continue
            try:
                jd, mag = julian_date(line[15:32]), float(line[65:70])
            except ValueError:
                continue
            station = line[77:80]
            if station == 'C57' or (2458356 < jd < 2458436 and station != 'I41'):
                photometry.append([number, f'{jd:.8f}', mag, line[70].strip(), station])
        params = dict(format='text', COMMAND=f"'{number};'", MAKE_EPHEM="'YES'",
             EPHEM_TYPE="'OBSERVER'", CENTER="'500@399'", START_TIME="'2018-07-01'",
             STOP_TIME="'2018-12-15'", STEP_SIZE="'1d'", QUANTITIES="'19,20,24'",
             CSV_FORMAT="'YES'", CAL_FORMAT="'JD'")
        horizon_url = HORIZONS+'?'+urllib.parse.urlencode(params)
        inputs['horizons'] = get(horizon_url, cache/f'{number}_horizons.txt' if cache else None)
        ephem = inputs['horizons'].decode().split('$$SOE')[1].split('$$EOE')[0].strip()
        for line in ephem.splitlines():
            row = line.split(',')
            geometry.append([number] + [row[j].strip() for j in (0, 3, 5, 7)])
        tess_url = TSSYS+str(number)+'.lc'
        inputs['tssys'] = get(tess_url, local_tess/f'{number}.lc' if local_tess else None)
        a = np.loadtxt(io.BytesIO(inputs['tssys']), usecols=(1, 4, 5, 9, 10, 11, 12))
        for t, m, err, flag, delta, r, phase in a:
            if flag == 0 and err < .15:
                spin_rows.append([number, t, m, delta, r, phase])
        manifest.append({'mpc': number, 'MPC_request': {'url': MPC, 'GET_json': payload},
             'Horizons_url': horizon_url, 'TSSYS_url': tess_url, 'source_sha256':
             {k:hashlib.sha256(v).hexdigest() for k,v in inputs.items()}})
    for name,header,data in [('observations','mpc,jd,mag,band,station',photometry),
        ('geometry','mpc,jd,r_au,delta_au,phase_deg',geometry),
        ('tess_spin_input','mpc,jd,mag,delta_au,r_au,phase_deg',spin_rows)]:
        text = io.StringIO()
        writer = csv.writer(text,lineterminator='\n')
        writer.writerow(header.split(','))
        writer.writerows(data)
        (HERE/(name+'.csv.gz')).write_bytes(gzip.compress(text.getvalue().encode(),mtime=0))
    (HERE/'sources.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-cache', type=Path)
    parser.add_argument('--tssys-directory', type=Path)
    args = parser.parse_args()
    freeze(args.raw_cache, args.tssys_directory)

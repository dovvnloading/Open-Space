"""Refine each preselected encounter against spacecraft-centred JPL Horizons."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar


def identifier(name):
    parts = name.split()
    if parts[0].isdigit() and (len(parts) > 2 or not 1900 <= int(parts[0]) <= 2026):
        return parts[0]
    return name


def refine(row, cache):
    query_id = identifier(row['name'])
    day = int(float(row['screen_jd_tdb']))
    query = {'format':'json','COMMAND':repr(query_id+';'),'CENTER':"'500@-28'",
             'EPHEM_TYPE':"'VECTORS'",'START_TIME':repr(f'JD{day-6.5}'),
             'STOP_TIME':repr(f'JD{day+6.5}'),'STEP_SIZE':"'3 h'",'OUT_UNITS':"'KM-S'",
             'VEC_TABLE':"'2'",'CSV_FORMAT':"'YES'"}
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urlencode(query)
    path = cache / f"-28_{query_id.replace(' ','_')}_{day}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url, timeout=90) as response:
            raw = response.read()
        if b'$$SOE' not in raw:
            raise ValueError(f'Horizons returned no states for {query_id}')
        path.write_bytes(raw)
    lines = json.loads(raw)['result'].split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()
    data = np.array([[float(r[0]), *map(float, r[2:8])] for r in csv.reader(lines)])
    t = (data[:,0]-day)*86400
    orbit = CubicHermiteSpline(t, data[:,1:4], data[:,4:7], axis=0)
    nearest = np.argmin(np.linalg.norm(data[:,1:4],axis=1))
    a, b = t[max(0,nearest-1)], t[min(len(t)-1,nearest+1)]
    opt = minimize_scalar(lambda dt: np.linalg.norm(orbit(dt)), bounds=(a,b),method='bounded')
    # Include endpoints in case the closest point lies outside the refinement interval.
    dt = min([a,opt.x,b],key=lambda x: np.linalg.norm(orbit(x)))
    return {**row,'query_id':query_id,'closest_jd_tdb':day+dt/86400,
            'distance_km':np.linalg.norm(orbit(dt)),'speed_km_s':np.linalg.norm(orbit(dt,1)),
            'at_window_edge':bool(dt in (t[0],t[-1])),
            'source_sha256':hashlib.sha256(raw).hexdigest()}


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('screen_csv',type=Path)
    p.add_argument('cache',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with args.screen_csv.open() as source:
        results = [refine(row,args.cache) for row in csv.DictReader(source)]
    with args.output.open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

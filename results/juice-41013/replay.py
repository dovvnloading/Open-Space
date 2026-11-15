"""Calculate geometric encounters and apparent observing conditions from frozen data."""
import csv
from datetime import datetime, timedelta
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

DATA = Path(__file__).resolve().parent / 'data'
REF_JD = 2461735.0  # 2027 November 25 12:00 TDB
DAY = 86400
AU = 149597870.7


def load_source(name):
    raw = gzip.decompress((DATA / (name+'.gz')).read_bytes())
    receipt = json.loads((DATA/'sources.json').read_text())[name]
    if hashlib.sha256(raw).hexdigest() != receipt['sha256']:
        raise ValueError(f'Hash mismatch for {name}')
    return json.loads(raw)


def horizons_rows(name):
    result = load_source(name)['result']
    return list(csv.reader(result.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))


def spline(jd, pos, vel):
    return CubicHermiteSpline((jd - REF_JD)*DAY, pos, vel, axis=0, extrapolate=False)


def horizons_orbit(name):
    rows = horizons_rows(name)
    state = np.array([[float(r[0]), *map(float, r[2:8])] for r in rows])
    return spline(state[:,0], state[:,1:4], state[:,4:])


def miriade_orbit(name):
    rows = load_source(name+'.json')['data']
    time = np.array([r['Date'] for r in rows])  # TT, treated as TDB (<2 ms)
    pos = np.array([[r[k] for k in ('px','py','pz')] for r in rows])*AU
    vel = np.array([[r[k] for k in ('vx','vy','vz')] for r in rows])*AU/DAY
    return spline(time, pos, vel)


def encounter(asteroid, spacecraft):
    def distance(t):
        return np.linalg.norm(asteroid(t)-spacecraft(t))
    t = minimize_scalar(distance, bounds=(-2*3600, 18*3600), method='bounded').x
    speed = np.linalg.norm(asteroid(t,1)-spacecraft(t,1))
    clock = datetime(2027,11,25,12) + timedelta(seconds=t)
    return {'tdb':clock.isoformat(timespec='seconds'),'distance_km':round(distance(t),2),
            'relative_km_s':round(speed,5)}, t


def measure():
    asteroid = horizons_orbit('aster-fine.json')
    craft = horizons_orbit('craft-fine.json')
    nominal, t = encounter(asteroid, craft)
    models = {key:encounter(miriade_orbit(key),craft)[0] for key in ('astorb','mpcorb')}
    with gzip.open(DATA/'esa-states.csv.gz','rt') as f:
        revisions = list(csv.DictReader(f))
    kernels = {}
    for name in dict.fromkeys(r['revision'] for r in revisions):
        rows = [r for r in revisions if r['revision']==name]
        state = np.array([[float(r[k]) for k in ('jd_tdb','x_km','y_km','z_km','vx_km_s','vy_km_s','vz_km_s')] for r in rows])
        kernels[name] = encounter(asteroid,spline(state[:,0],state[:,1:4],state[:,4:]))[0]

    obs = horizons_rows('obs-fine.json')
    brightest = min(obs,key=lambda r:float(r[7]))
    nearby = next(r for r in obs if r[0].strip().endswith('19:20'))
    obs_speed = np.hypot(float(nearby[5]),float(nearby[6]))/3600  # arcsec/s
    widths = {}
    for albedo in (0.04,0.14,0.40):
        diameter = 1329*10**(-16.21/5)/np.sqrt(albedo)
        widths[albedo] = {'diameter_km':round(diameter,3),
                          'samples_across':round(diameter/(nominal['distance_km']*15e-6),2)}
    def threshold(x):
        return np.linalg.norm(asteroid(x)-craft(x))-100000
    hours = (brentq(threshold,t,min(asteroid.x[-1],craft.x[-1]))-brentq(threshold,max(asteroid.x[0],craft.x[0]),t))/3600
    with (DATA/'encounters.csv').open() as f:
        candidates = sorted(csv.DictReader(f),key=lambda r:float(r['distance_km']))
    return {'nominal':nominal,'hours_inside_100000_km':round(hours,3),'asteroid_solution_checks':models,
            'spacecraft_history':kernels, 'apparent_V_brightest_sample':float(brightest[7]),
            'brightest_sample_utc':brightest[0].strip(),'sun_elong_deg_at_brightest':float(brightest[12]),
            'phase_deg_at_brightest':float(brightest[14]),'sky_motion_arcsec_s_1920utc':round(obs_speed,2),
            'three_sigma_asteroid_range_km_1920utc':float(nearby[22]),
            'sampled_diameter_by_assumed_albedo':widths, 'refined_count':len(candidates),
            'two_closest':[{'name':r['name'],'distance_km':round(float(r['distance_km']))} for r in candidates[:2]]}


if __name__ == '__main__':
    print(json.dumps(measure(),indent=2))

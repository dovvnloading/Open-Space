"""Replay the frozen nomination, or scan cached/downloadable orbital catalogues."""
import argparse
import csv
from datetime import date
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
K = 0.01720209895
AU_KM = 149597870.7
GRID = 2460676.5 + 8 * np.arange(183)


def basis(inclination, perihelion, node):
    inc, w, om = np.deg2rad([inclination, perihelion, node])
    p = np.stack((np.cos(w)*np.cos(om) - np.sin(w)*np.sin(om)*np.cos(inc),
                  np.cos(w)*np.sin(om) + np.sin(w)*np.cos(om)*np.cos(inc),
                  np.sin(w)*np.sin(inc)), axis=-1)
    q = np.stack((-np.sin(w)*np.cos(om) - np.cos(w)*np.sin(om)*np.cos(inc),
                  -np.sin(w)*np.sin(om) + np.cos(w)*np.cos(om)*np.cos(inc),
                  np.cos(w)*np.sin(inc)), axis=-1)
    return p, q


def numbered_elements(lines):
    rows = np.array([[int(line[:6]), *map(float, line[25:].split()[:8])] for line in lines])
    num, epoch, a, e, inc, w, om, mean, h = rows.T
    mask = (epoch > 59000) & (inc < 180) & (a < 10) & (e < 0.98)
    num, epoch, a, e, inc, w, om, mean, h = rows[mask].T
    p, q = basis(inc, w, om)
    return num.astype(int), epoch, np.deg2rad(mean), e, p*a[:, None], q*(a*np.sqrt(1-e*e))[:, None], K/a**1.5


def numbered_state(elements, t):
    _, epoch, mean, e, p, q, motion = elements
    m = mean + motion*(t - 2400000.5 - epoch)
    eccentric = m + e*np.sin(m)*(1 + e*np.cos(m))
    for _ in range(3):
        eccentric -= (eccentric - e*np.sin(eccentric) - m)/(1 - e*np.cos(eccentric))
    c, s = np.cos(eccentric), np.sin(eccentric)
    position = (c-e)[:, None]*p + s[:, None]*q
    velocity = (-s[:, None]*p + c[:, None]*q)*(motion/(1-e*c))[:, None]
    return position, velocity


def comet_track(line):
    name = line[102:159].strip()
    if not name or not line[81:89].strip():
        return None
    epoch = date.fromisoformat(line[81:89]).toordinal() - date(1858, 11, 17).toordinal()
    q, e, w, node, inc = [float(line[x:y]) for x, y in ((30,39),(41,49),(51,59),(61,69),(71,79))]
    if epoch < 58800 or q > 4.5 or e == 1 or (name[0].isdigit() and '-' in name.split('/')[0]):
        return None
    tp = date(int(line[14:18]), int(line[19:21]), 1).toordinal() + 1721424.5 + float(line[22:29]) - 1
    axis = q/abs(1-e)
    n = K/axis**1.5
    m = n*(GRID-tp)
    if e < 1:
        m = (m+np.pi) % (2*np.pi) - np.pi
        u = m + np.sign(m)*0.75*e
        for _ in range(15):
            u -= (u-e*np.sin(u)-m)/(1-e*np.cos(u))
        c, s = np.cos(u), np.sin(u)
        x, y = axis*(c-e), axis*np.sqrt(1-e*e)*s
        udot = n/(1-e*c)
        vx, vy = -axis*s*udot, axis*np.sqrt(1-e*e)*c*udot
    else:
        u = np.arcsinh(m/e) + np.sign(m)*0.75
        for _ in range(25):
            u -= (e*np.sinh(u)-u-m)/(e*np.cosh(u)-1)
        c, s = np.cosh(u), np.sinh(u)
        x, y = axis*(e-c), axis*np.sqrt(e*e-1)*s
        udot = n/(e*c-1)
        vx, vy = -axis*s*udot, axis*np.sqrt(e*e-1)*c*udot
    p, v = basis(inc, w, node)
    name = name.split(' (')[0] if name.startswith(('C/','P/','A/','I/')) else name.split('/')[0]
    return name, (x[:, None]*p+y[:, None]*v, vx[:, None]*p+vy[:, None]*v)


def cached(cache, filename, url):
    path = cache / filename
    if not path.exists():
        with urlopen(url, timeout=180) as response:
            payload = response.read()
        path.write_bytes(payload)
    return path


def population(cache, inputs):
    cache.mkdir(exist_ok=True, parents=True)
    paths = []
    for k, name in [('numbered_asteroids', 'elements.gz'), ('mpc_comets', 'mpccom.txt')]:
        path = cached(cache, name, inputs[k]['url'])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != inputs[k]['sha256']:
            raise ValueError(f'{name} has changed: {digest}. Supply the recorded snapshot to replay.')
        paths.append(path)
    with gzip.open(paths[0], 'rt') as stream:
        next(stream)
        next(stream)
        asteroids = numbered_elements(stream)
    tracks = dict(track for line in paths[1].read_text().splitlines() if (track := comet_track(line)))
    for row in csv.DictReader((ROOT/'comet_grid_sources.csv').open()):
        target = row['comet']
        path = cached(cache, f'{target}.json', row['query'])
        text = json.loads(path.read_text())['result']
        raw = csv.reader(text.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines())
        vectors = np.array([[float(r[0]), *map(float, r[2:8])] for r in raw])
        if not np.array_equal(vectors[:, 0], GRID):
            raise ValueError(f'Unexpected epochs for {target}')
        tracks[target] = vectors[:, 1:4], vectors[:, 4:7]
    return asteroids, tracks


def screen(elements, tracks):
    from scipy.spatial import cKDTree
    best = {}
    for j, t in enumerate(GRID):
        r, v = numbered_state(elements, t)
        tree = cKDTree(r)
        for name, (rc, vc) in tracks.items():
            match = tree.query_ball_point(rc[j], 0.19)
            if not match:
                continue
            dr, dv = r[match] - rc[j], v[match] - vc[j]
            with np.errstate(divide='ignore', invalid='ignore'):
                dt = -np.sum(dr*dv, axis=1)/np.sum(dv*dv, axis=1)
            distance = np.linalg.norm(dr + dt[:,None]*dv, axis=1)
            for i in np.flatnonzero((abs(dt) < 4.05) & (distance < .006)):
                key = name, int(elements[0][match[i]])
                candidate = (distance[i]*AU_KM, t+dt[i], name, key[1], np.linalg.norm(dv[i])*1731.456,
                             np.linalg.norm(rc[j]))
                if key not in best or candidate[0] < best[key][0]:
                    best[key] = candidate
    return sorted(best.values())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--population', type=Path, metavar='CACHE', help='Scan both full catalogue snapshots')
    parser.add_argument('--output', type=Path, help='Optional candidate CSV path')
    args = parser.parse_args()
    inputs = json.loads((ROOT/'ephemerides.json').read_text())['search_inputs']
    if args.population:
        elements, tracks = population(args.population, inputs)
    else:
        subset = inputs['element_records']
        elements = numbered_elements([subset['asteroid']])
        tracks = dict([comet_track(subset['comet'])])
    hits = screen(elements, tracks)
    follow = [r for r in hits if r[0] < 110000 and r[4] > .4 and r[5] < 3.8]
    print('asteroids, comet tracks, coarse nominations, follow-up eligible',
          len(elements[0]), len(tracks), len(hits), len(follow))
    for distance, jd, comet, asteroid, speed, rh in hits:
        if comet == '407P' and asteroid == 75547:
            print(f'407P / (75547): JD {jd:.8f}, {distance:.1f} km, {speed:.4f} km/s')
    if args.output:
        with args.output.open('w', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(['screen_km','jd_tdb','comet','asteroid','speed_km_s','comet_r_at_grid_au'])
            writer.writerows(hits)

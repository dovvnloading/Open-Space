"""Recreate the Mars moon penumbra screen and the October 2027 shape checks."""
import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import spiceypy as spice
from PIL import Image, ImageDraw
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
SUN_RADIUS = 695700.0  # km, IAU nominal
PHOBOS_RADIUS = 11.08  # km, approximate sphere for initial contact only
SHAPE = 'phobos_2014_09_22.bds'
NOE = 'NOE-4-2020.bsp'
SITES = [(-4.82602331, 137.388886), (18.41996953, 77.22012875)]  # planetocentric lat, east lon
EVENT_DAYS = ['2027-10-08T17:29', '2027-10-09T23:49', '2027-10-11T06:13']


def kernels(cache, download):
    manifest = json.loads((HERE / 'kernels.json').read_text())
    cache.mkdir(parents=True, exist_ok=True)
    for name, item in manifest.items():
        path = cache / name
        if download and not path.exists():
            part = path.with_suffix('.part')
            with part.open('wb') as out:
                if name == NOE:
                    out.truncate(item['size'])
                    for span in item['ranges']:
                        offset, n = span['offset'], span['length']
                        endpoint = item['url'] + f'?start={offset}'
                        with urlopen(Request(endpoint, headers={'Range': f'bytes={offset}-{offset+n-1}'}), timeout=60) as reply:
                            if reply.headers.get('Content-Range') != f'bytes {offset}-{offset+n-1}/{item["size"]}':
                                raise ValueError('Server returned an unexpected byte range')
                            body = reply.read(n+1)
                        if hashlib.sha256(body).hexdigest() != span['sha256']:
                            raise ValueError(f'Incorrect NOE bytes at {offset}')
                        out.seek(offset)
                        out.write(body)
                else:
                    with urlopen(item['url'], timeout=180) as reply:
                        shutil.copyfileobj(reply, out)
            part.replace(path)
        with path.open('rb') as f:
            if name == NOE:
                if path.stat().st_size != item['size']:
                    raise ValueError('Invalid NOE kernel size')
                for span in item['ranges']:
                    f.seek(span['offset'])
                    if hashlib.sha256(f.read(span['length'])).hexdigest() != span['sha256']:
                        raise ValueError(f'Incorrect NOE bytes at {span["offset"]}')
            elif hashlib.file_digest(f, 'sha256').hexdigest() != item['sha256']:
                raise ValueError(f'Unexpected checksum for {name}')
        if name not in (SHAPE, NOE):
            spice.furnsh(str(path))


def norm(v):
    return np.linalg.norm(v, axis=-1)


def unit(v):
    return v / norm(v)[..., None]


def separation(v, w):
    return np.arctan2(norm(np.cross(v, w)), np.sum(v * w, axis=-1))


def pos(target, t, observer, correction='CN+S'):
    return np.asarray(spice.spkpos(str(target), t, 'J2000', correction, str(observer))[0])


def geometry(t, receiver=402, blocker=401):
    sun = pos(10, t, receiver)
    moon = pos(blocker, t, receiver)
    mars = pos(499, t, receiver)
    distance = norm(moon)
    r = PHOBOS_RADIUS if blocker == 401 else 6.2
    return (separation(sun, moon), math.asin(SUN_RADIUS/norm(sun)), math.asin(r/distance),
            separation(sun, mars) - math.asin(3396.2/norm(mars)), distance)


def altitude(t, body, lat, lon):
    site = spice.srfrec(499, np.radians(lon), np.radians(lat))
    ray, _ = spice.spkcpo(str(body), t, 'IAU_MARS', 'OBSERVER', 'CN+S', site, 'MARS', 'IAU_MARS')
    return np.degrees(math.asin(np.dot(unit(site), unit(ray[:3]))))


def closest(t, receiver, blocker, span):
    return t + minimize_scalar(lambda dt: geometry(t+dt, receiver, blocker)[0],
                               bounds=(-span, span), method='bounded').x


def screen(step):
    start, end = (spice.str2et(d) for d in ['2026-11-15T00:00', '2033-01-01T00:00'])
    events = []
    for lo in np.arange(start, end, 86400*60):
        hi = min(end, lo+86400*60)
        ts = np.arange(lo-step, hi+2*step, step)
        phobos = pos(401, ts, 499, 'NONE')
        deimos = pos(402, ts, 499, 'NONE')
        sun = unit(pos(10, ts, 499, 'NONE'))
        relative = phobos - deimos
        impact = norm(np.cross(relative, sun))
        indices = np.where((impact[1:-1] <= impact[:-2]) & (impact[1:-1] < impact[2:]) & (impact[1:-1]<750))[0]+1
        for i in indices:
            receiver, blocker = (402,401) if np.dot(relative[i],sun[i])>0 else (401,402)
            t = closest(ts[i], receiver, blocker, step)
            sep, solar, moon, clearance, distance = geometry(t, receiver, blocker)
            if lo <= t < hi and sep < solar+moon and clearance > solar:
                events.append([utc(t), receiver, np.degrees(sep)*60, np.degrees(solar)*60,
                               np.degrees(moon)*60, np.degrees(clearance), distance])
    return events


def utc(t):
    return spice.et2utc(t, 'ISOC', 3)


def projected_area(t, vertices, plates, resolution):
    moon, lt = spice.spkpos('PHOBOS', t, 'J2000', 'CN', 'DEIMOS')
    z = unit(moon)
    x = unit(np.cross([0,0,1], z))
    basis = np.vstack([x, np.cross(z,x), z])
    rotation = spice.pxform('IAU_PHOBOS', 'J2000', t-lt)
    shape = vertices @ (rotation.T @ basis.T)
    shape[:,2] += norm(moon)
    coordinates = shape[:,:2]/shape[:,2:]
    sun = pos(10, t, 402, 'CN') @ basis.T
    sun_radius = math.tan(math.asin(SUN_RADIUS/norm(sun)))
    lower = coordinates.min(axis=0)
    pitch = (coordinates.max(axis=0)-lower).max()/(resolution-4)
    pixels = np.rint((coordinates-lower)/pitch).astype(int)+2
    image = Image.new('1',(resolution,resolution))
    draw = ImageDraw.Draw(image)
    for tri in pixels[plates].reshape(-1,6):
        draw.polygon(tuple(tri),fill=1)
    yy,xx = np.nonzero(np.asarray(image))
    xy = np.column_stack((xx-2,yy-2))*pitch+lower
    covered = np.count_nonzero(norm(xy-sun[:2]/sun[2]) < sun_radius)
    return 100*covered*pitch**2/(np.pi*sun_radius**2)


def detail(model, cache):
    handle = spice.dasopr(str(cache / SHAPE))
    segment = spice.dlabfs(handle)
    nv, nt = spice.dskz02(handle, segment)
    vertices = spice.dskv02(handle, segment, 1, nv)
    plates = spice.dskp02(handle, segment, 1, nt)-1
    spice.dascls(handle)
    phobos_bound = norm(vertices).max()+0.001
    rows = []
    for date in EVENT_DAYS:
        t = closest(spice.str2et(date),402,401,240)
        sep, sun, _, clearance, distance = geometry(t)

        def spherical_gap(dt):
            a, b, c, _, _ = geometry(t+dt)
            return a-b-c

        def all_sites_margin(dt):
            a, b, _, _, d = geometry(t+dt)
            # Conservative geometric parallax for the 7.8-km bounding ellipsoid.
            # Use a conservative lower bound of 200 million km for this season's Sun distance.
            assert norm(pos(10,t+dt,402)) > 2e8
            return b-a-math.asin(phobos_bound/(d-7.8))-math.asin(7.8/d)-math.asin(7.8/2e8)

        ingress, egress = [brentq(spherical_gap,*bracket) for bracket in [(-900,0),(0,900)]]
        first, last = [brentq(all_sites_margin,*bracket) for bracket in [(-900,0),(0,900)]]
        skies = [altitude(t,b,*site) for site in SITES for b in (402,10)]
        areas = [projected_area(t,vertices,plates,resolution) for resolution in (512,1024)]
        rows.append([model,utc(t),np.degrees(sep)*60,distance,np.degrees(clearance),*skies,utc(t+ingress),
                     utc(t+egress),egress-ingress,utc(t+first),utc(t+last),last-first,*areas])
    return rows


def write_csv(path, header, rows):
    with path.open('w',newline='') as f:
        writer = csv.writer(f,lineterminator='\n')
        writer.writerow(header.split())
        for row in rows:
            writer.writerow([f'{x:.6f}' if isinstance(x,float) else x for x in row])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--download',action='store_true')
    p.add_argument('--screen',action='store_true')
    p.add_argument('--step',type=int,choices=[150,300],default=300)
    args = p.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    kernels(args.cache,args.download)
    if args.screen:
        write_csv(args.output/'screen.csv','target_center_utc shadowed_body_naif sep_arcmin solar_radius_arcmin blocker_radius_arcmin mars_limb_clearance_deg range_km',screen(args.step))
    rows = detail('MAR099',args.cache)
    spice.furnsh(str(args.cache/NOE))
    rows += detail('NOE-4-2020',args.cache)
    write_csv(args.output/'contacts.csv','orbit target_center_utc sep_arcmin distance_km mars_limb_clearance_deg gale_deimos_alt_deg gale_sun_alt_deg jezero_deimos_alt_deg jezero_sun_alt_deg sphere_first_utc sphere_last_utc sphere_duration_s full_all_first_utc full_all_last_utc full_all_duration_s area_pct_512 area_pct_1024',rows)


if __name__ == '__main__':
    main()

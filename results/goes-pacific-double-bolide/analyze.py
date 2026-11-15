"""Reconstruct the Pacific GLM doublet and audit temporal coincidences in the stereo catalog."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime
from itertools import combinations, product
from pathlib import Path
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
A = 6378.137  # WGS84 equatorial radius, km
B = A * (1 - 1 / 298.257223563)
SAT_LONGITUDES = [-75.1999969, -137.0]
NOMINAL_HEIGHT = 35786.0234


def geodetic_point(lon, lat, major, minor):
    lon, lat = np.radians([lon, lat])
    normal = major**2 / np.hypot(major * np.cos(lat), minor * np.sin(lat))
    return normal * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                              (minor / major)**2 * np.sin(lat)])


def to_geodetic(point):
    horizontal = np.hypot(*point[:2])
    eccentricity_sq = 1 - (B / A)**2
    latitude = np.arctan2(point[2], horizontal * (1 - eccentricity_sq))
    for _ in range(8):
        normal = A / np.sqrt(1 - eccentricity_sq * np.sin(latitude)**2)
        latitude = np.arctan2(point[2] + eccentricity_sq * normal * np.sin(latitude), horizontal)
    return [np.degrees(np.arctan2(point[1], point[0])), np.degrees(latitude),
            horizontal / np.cos(latitude) - normal]


def triangulate(east, west):
    # Ray origins use the satellite positions recorded in both NASA products.
    stations = [geodetic_point(lon, 0, A + NOMINAL_HEIGHT, B + NOMINAL_HEIGHT)
                for lon in SAT_LONGITUDES]
    sight = [geodetic_point(*coord, A + 14, B + 6) - station
             for coord, station in zip([east, west], stations)]
    scale = np.linalg.lstsq(np.array([sight[0], -sight[1]]).T,
                            stations[1] - stations[0], rcond=None)[0]
    intersections = np.asarray(stations) + scale[:, None] * np.asarray(sight)
    return intersections.mean(axis=0), np.linalg.norm(intersections[1] - intersections[0])


def read_flash(index, sat):
    path = HERE / f'flash-{index}-G{sat}.csv'
    return np.loadtxt(path, delimiter=',', skiprows=11)


def bins(measurements):
    overlap = max(data[0, 0] for data in measurements), min(data[-1, 0] for data in measurements)
    for t in range(int(overlap[0]), int(overlap[1]) + 1, 20):
        groups = [data[(data[:, 0] >= t) & (data[:, 0] < t + 20) &
                       (data[:, 0] <= overlap[1])] for data in measurements]
        if not all(len(group) for group in groups):
            continue
        centroids = [np.average(group[:, 1:3], axis=0, weights=group[:, 3]) for group in groups]
        point, miss = triangulate(*centroids)
        strength = math.sqrt(np.prod([group[:, 3].sum() for group in groups]))
        yield t, point, miss, strength, centroids


def error_envelope(centroids):
    # Each satellite's reported footprint is offset independently in both axes.
    # Evaluate every corner and centre/edge of two +/-10 km coordinate boxes.
    points = []
    for dx1, dy1, dx2, dy2 in product([-10., 0., 10.], repeat=4):
        perturbed = []
        for (lon, lat), dx, dy in zip(centroids, [dx1, dx2], [dy1, dy2]):
            perturbed.append([lon + dx / (111.32 * np.cos(np.radians(lat))), lat + dy / 111.32])
        points.append(triangulate(*perturbed)[0])
    return points


def reconstruct():
    peaks, envelopes, bounds = [], [], []
    for index in (1, 2):
        measurements = [read_flash(index, sat) for sat in (16, 18)]
        series = list(bins(measurements))
        time, point, miss, _, centroids = max(series, key=lambda row:row[3])
        lon, lat, height = to_geodetic(point)
        peaks.append(point)
        envelopes.append(error_envelope(centroids))
        bounds.append([min(d[0, 0] for d in measurements), max(d[-1, 0] for d in measurements)])
        heights = [to_geodetic(p)[2] for p in envelopes[-1]]
        print(f'Flash {index}, peak stereo bin {time}, lon {lon:.3f}, lat {lat:.3f}, '
              f'height {height:.1f} km, ray miss {miss:.1f} km. '
              f'Offset grid heights {min(heights):.1f} to {max(heights):.1f} km.')
        print('  Stereo heights by 20 ms bin:', [round(float(to_geodetic(row[1])[2]), 1) for row in series])
    distances = [np.linalg.norm(a-b) for a,b in product(*envelopes)]
    dt = (bounds[1][1] - bounds[0][0]) / 1000
    print(f'Peak locations {np.linalg.norm(peaks[1]-peaks[0]):.1f} km apart. '
          f'Grid separation {min(distances):.1f} to {max(distances):.1f} km. '
          f'Full duration {dt:.3f} s, min grid separation / full duration {min(distances)/dt:.1f} km/s.')


def load_catalog():
    with (HERE / 'stereo-catalog.csv').open() as handle:
        rows = sorted(csv.DictReader(handle), key=lambda r:int(r['utc_onset_ms']))
    times = np.array([int(r['utc_onset_ms']) for r in rows]) / 1000
    lat, lon = np.radians([[float(r[k]) for r in rows] for k in ('latitude_deg', 'longitude_deg')])
    xyz = np.array([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]).T
    return rows, times, xyz


def arc_dist(xyz, i, j):
    return 2 * 6371 * np.arcsin(np.minimum(1, np.linalg.norm(xyz[i] - xyz[j], axis=-1) / 2))


def audit():
    rows, times, xyz = load_catalog()
    close = [(i,j) for j,t in enumerate(times) for i in range(np.searchsorted(times,t-10),j)]
    print(f'{len(rows)} stereo events. All times separated by at most 10 s:')
    for i,j in close:
        print(f'  {rows[i]["event_id"]}  {rows[j]["event_id"]}  {times[j]-times[i]:.3f} s, '
              f'{arc_dist(xyz,i,j):.0f} km (catalog positions).')
    # Conditional expectation over independent permutations of position labels in each block.
    # Keep the exact time sequence, local position distribution, and storm-day event counts.
    for width in [600, 1800, 3600, 7200, 86400]:
        for offset in [0, width/2]:
            groups = defaultdict(list)
            blocks = np.floor((times+offset)/width).astype(int)
            for i,block in enumerate(blocks):
                groups[block].append(i)
            expected = 0.
            for i,j in close:
                if times[j]-times[i] > 2.5:
                    continue
                li,lj = groups[blocks[i]],groups[blocks[j]]
                pairs = list(combinations(li,2)) if blocks[i]==blocks[j] else list(product(li,lj))
                left,right = np.array(pairs).T
                expected += np.mean(arc_dist(xyz,left,right) <= 350)
            print(f'  {width:g}s blocks, {offset:g}s phase: conditional expected count {expected:.5f}')


def verify_noaa(meta, directory):
    from netCDF4 import Dataset
    quality_filtered = defaultdict(list)
    for sat, product_meta in meta['noaa_l2'].items():
        path = directory / product_meta['url'].split('/')[-1]
        if not path.exists():
            path.write_bytes(urlopen(product_meta['url'],timeout=60).read())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == product_meta['sha256']
        with Dataset(path) as image:
            reference = datetime.fromisoformat(image.time_coverage_start.replace('Z','+00:00')).timestamp()
            optical = np.array([image['group_time_offset'][:],image['group_lon'][:],
                                image['group_lat'][:],image['group_energy'][:]]).astype(float).T
            optical[:,0] = np.floor(optical[:,0]*1000) + round(reference*1000)
            for index in (1,2):
                data = read_flash(index,sat)
                # Four columns must agree for each optical group from the corresponding telescope.
                selected = np.ones(len(optical), bool)
                for col, margin in [(0,1.),(1,1.e-5),(2,1.e-5)]:
                    selected &= (optical[:,col]>=min(data[:,col])-margin)&(optical[:,col]<=max(data[:,col])+margin)
                original = optical[selected]
                assert original.shape == data.shape
                assert np.allclose(data[:,0],original[:,0],rtol=0,atol=1.)
                assert np.allclose(data[:,1:3],original[:,1:3],rtol=0,atol=1.e-5)
                assert np.allclose(data[:,3],original[:,3],rtol=1.e-5,atol=0)
                good = np.asarray(image['group_quality_flag'][:][selected]) == 0
                quality_filtered[index].append(data[good])
                print(f'NOAA G{sat}, flash {index}: {len(data)} groups verified, {good.sum()} good quality.')
    for index in (1, 2):
        full = [read_flash(index, sat) for sat in meta['noaa_l2']]
        peak = max(bins(full), key=lambda row: row[3])
        clean = max(bins(quality_filtered[index]), key=lambda row: row[3])
        assert peak[0] == clean[0] and np.linalg.norm(peak[1] - clean[1]) < 1.e-6



def plot(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    figure,axes=plt.subplots(2,2,figsize=(8.6,5.4),sharex='col',layout='constrained')
    origin=1700583076000
    for col,index in enumerate((1,2)):
        data=[read_flash(index,sat) for sat in (16,18)]
        top,bottom=axes[:,col]
        for samples,sat,color in zip(data,(16,18),('#22669C','#BA5527')):
            # Group energy per 10 ms bin, normalized separately for each satellite and flash.
            edges=np.arange(np.floor(samples[0,0]/10)*10,np.ceil(samples[-1,0]/10)*10+11,10.)
            energy,_=np.histogram(samples[:,0],edges,weights=samples[:,3])
            top.stairs(energy/energy.max(),(edges-origin)/1000,color=color,label=f'GOES-{sat}')
        track=list(bins(data))
        bottom.plot([(v[0]+10-origin)/1000 for v in track],[to_geodetic(v[1])[2] for v in track],
                    'o-',ms=4,color='#233C51',lw=1.2)
        top.set_title(f'Flash {index}',loc='left',fontweight='bold')
        top.set_ylim(-0.03,1.12)
        top.legend(loc='upper left',frameon=False,fontsize=8)
        bottom.set_ylim(78,101)
        bottom.set_xlabel('Seconds after 16:11:16 UTC')
        bottom.xaxis.set_major_locator(MultipleLocator(.1))
        for axis in (top,bottom):
            axis.grid(axis='y',color='#999999',alpha=.25,lw=.7)
    axes[0,0].set_ylabel('Normalized energy / 10 ms')
    axes[1,0].set_ylabel('Stereo height (km)')
    figure.suptitle('Two Pacific bolides on 21 November 2023',fontweight='bold')
    figure.savefig(output)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--noaa-cache',type=Path,help='verify the underlying NOAA granules, requires netCDF4')
    parser.add_argument('--plot',type=Path,help='draw time series and stereo heights, requires matplotlib')
    arguments=parser.parse_args()
    metadata=json.loads((HERE/'sources.json').read_text())
    for i,event in enumerate(metadata['events'],1):
        for sat,resource in event['products'].items():
            assert hashlib.sha256((HERE/f'flash-{i}-G{sat}.csv').read_bytes()).hexdigest()==resource['sha256']
    reconstruct()
    audit()
    if arguments.plot:
        plot(arguments.plot)
    if arguments.noaa_cache:
        arguments.noaa_cache.mkdir(parents=True,exist_ok=True)
        verify_noaa(metadata,arguments.noaa_cache)

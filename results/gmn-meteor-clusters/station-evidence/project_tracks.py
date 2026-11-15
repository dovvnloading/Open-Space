"""Project the published endpoints into a 320x180 thumbnail using a simple equidistant lens."""
import argparse
import csv
from pathlib import Path

import numpy as np


def ecef(lon, lat, height_m):
    lon, lat = np.radians([lon, lat])
    eccentricity_squared = 0.00669437999
    radius = 6378137 / np.sqrt(1 - eccentricity_squared * np.sin(lat)**2)
    return np.array([(radius + height_m) * np.cos(lat) * np.cos(lon),
                     (radius + height_m) * np.cos(lat) * np.sin(lon),
                     (radius * (1 - eccentricity_squared) + height_m) * np.sin(lat)])


def project(camera, track, suffix):
    lon, lat, height, az, el, roll, width, height_fov = [
        float(camera[col]) for col in ('longitude_deg', 'latitude_deg', 'height_m', 'azimuth_deg',
                                       'altitude_deg', 'roll_h_deg', 'fov_h_deg', 'fov_v_deg')]
    origin = ecef(lon, lat, height)
    site_lon, site_lat = np.radians([lon, lat])
    east = [-np.sin(site_lon), np.cos(site_lon), 0]
    north = [-np.sin(site_lat)*np.cos(site_lon), -np.sin(site_lat)*np.sin(site_lon), np.cos(site_lat)]
    up = [np.cos(site_lat)*np.cos(site_lon), np.cos(site_lat)*np.sin(site_lon), np.sin(site_lat)]
    obj = ecef(float(track['lon'+suffix]), float(track['lat'+suffix]), 1000 * float(track['h'+suffix]))
    line = np.array([east, north, up]) @ (obj - origin)
    line /= np.linalg.norm(line)
    az, el, roll, width, height_fov = np.radians([az, el, roll, width, height_fov])
    center = [np.cos(el)*np.sin(az), np.cos(el)*np.cos(az), np.sin(el)]
    right = [np.cos(az), -np.sin(az), 0]
    dx, dy = line @ right, line @ np.cross(right, center)
    angle = np.arccos(np.clip(line @ center, -1, 1))
    scale = angle / np.sin(angle) if angle > 1e-8 else 1
    x = 160 + (np.cos(roll)*dx + np.sin(roll)*dy) * scale * 320 / width
    y = 90 - (-np.sin(roll)*dx + np.cos(roll)*dy) * scale * 180 / height_fov
    return np.array([x, y])


def main():
    here = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tracks', type=Path, default=here.parent / 'tracks.csv')
    args = parser.parse_args()
    with args.tracks.open(newline='') as stream:
        tracks = {r['id']: r for r in csv.DictReader(stream)}
    with (here / 'camera-tiles.csv').open(newline='') as stream:
        tiles = {r['station']: r for r in csv.DictReader(stream)}
    with (here / 'pointing.csv').open(newline='') as stream:
        cameras = list(csv.DictReader(stream))
    print('camera,projected_start_x,projected_start_y,projected_end_x,projected_end_y,feature_to_segment_px')
    for camera in cameras:
        station = camera['station']
        track = tracks[camera['trajectory_id']]
        start, end = [project(camera, track, suffix) for suffix in ('', '2')]
        tile = tiles[station]
        distance = ''
        if tile['feature_x_px']:
            point = np.array([float(tile['feature_x_px']), float(tile['feature_y_px'])])
            segment = end - start
            fraction = np.clip((point - start) @ segment / (segment @ segment), 0, 1)
            distance = f'{np.linalg.norm(point - start - fraction * segment):.1f}'
        values = ','.join(f'{x:.1f}' for x in [*start, *end])
        print(f'{station},{values},{distance}')


if __name__ == '__main__':
    main()

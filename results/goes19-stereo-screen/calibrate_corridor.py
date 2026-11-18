"""Check nominal angular residuals on the two already published 2023 bolides."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np

from corridor import altitude_arc, angular_residual
from geometry import A, B, surface


INPUT_HASHES = {
    'flash-1-G16.csv': '9399c867cb4b50f358360c59f831b9b720b5886b310c7eabf48d6533faa52b40',
    'flash-1-G18.csv': '73581981648c2737a7ba76ffa77192c0045fe24928bc68cd7a76863aea1be675',
    'flash-2-G16.csv': 'ff2a3f2043bf28f10625624499ce4d56c5f827f223b6f1ea0900cf546d3c2034',
    'flash-2-G18.csv': '38599671a42419fbadc21a9e97038e9e28ca2c7c010964d38afa0db447aa82e5',
}
BIN_MS = 20
HEIGHTS = (40., 120.)
STATION_HEIGHT = 35786.0234


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_controls(directory):
    controls = {}
    for name, expected in INPUT_HASHES.items():
        raw = (directory / name).read_bytes()
        if sha(raw) != expected:
            raise ValueError(f'control digest mismatch: {name}')
        groups = np.loadtxt(io.BytesIO(raw), delimiter=',', skiprows=11, ndmin=2)
        if groups.shape[1] != 4 or not np.isfinite(groups).all() or np.any(groups[:, 3] <= 0):
            raise ValueError(f'invalid control groups: {name}')
        controls[name] = groups
    return controls


def evaluate(east, west):
    east_station = surface(-75.1999969, 0, A + STATION_HEIGHT, B + STATION_HEIGHT)
    west_station = surface(-137.0, 0, A + STATION_HEIGHT, B + STATION_HEIGHT)
    east_bin = east[:, 0].astype(np.int64) // BIN_MS
    bins, inverse = np.unique(east_bin, return_inverse=True)
    selected = int(bins[np.argmax(np.bincount(inverse, weights=east[:, 3]))])
    start = selected * BIN_MS
    first = east[east_bin == selected]
    second = west[(west[:, 0] >= start) & (west[:, 0] < start + BIN_MS)]
    east_centroid = np.average(first[:, 1:3], weights=first[:, 3], axis=0)
    arc = altitude_arc(east_station, west_station, *east_centroid, HEIGHTS)
    result = {'bin_start_unix_ms': start, 'east_groups': len(first), 'west_groups': len(second),
              'east_centroid_deg': east_centroid.tolist(), 'west_arc_endpoints_ecef_unit': [p.tolist() for p in arc]}
    if not len(second):
        return dict(result, status='no_west_groups_in_bin')
    west_centroid = np.average(second[:, 1:3], weights=second[:, 3], axis=0)
    direction = surface(*west_centroid, A + 14, B + 6) - west_station
    return dict(result, status='measured', west_centroid_deg=west_centroid.tolist(),
                arc_length_arcsec=3600 * np.degrees(angular_residual((arc[0], arc[0]), arc[1])),
                nearest_arc_residual_arcsec=3600 * np.degrees(angular_residual(arc, direction)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    control_path, output = args.controls.resolve(), args.output.resolve()
    if any(path.drive.upper() != 'D:' for path in (control_path, output)) or output.exists():
        raise ValueError('use a D: control path and a new D: output file')
    controls = read_controls(control_path)
    code = Path(__file__).resolve().parent
    result = {
        'scope': 'Retrospective calibration on 2023 published bolide groups. No 2026 values read.',
        'runtime': {'python': sys.version, 'numpy': np.__version__},
        'input_sha256': INPUT_HASHES,
        'code_sha256': {name: sha((code / name).read_bytes()) for name in ('calibrate_corridor.py', 'corridor.py', 'geometry.py')},
        'configuration': {'bin_width_ms': BIN_MS, 'altitude_shell_offsets_km': HEIGHTS,
                          'bin_selection': 'Max GOES-16 energy on Unix UTC-aligned bins, earliest on ties. Compare GOES-18 in same bin.'},
        'controls': {i: evaluate(controls[f'flash-{i}-G16.csv'], controls[f'flash-{i}-G18.csv']) for i in (1, 2)},
    }
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(output)


if __name__ == '__main__':
    main()

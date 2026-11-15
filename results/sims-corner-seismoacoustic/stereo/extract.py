"""Extract GOES lightning groups and verify the frozen NOAA granules."""
import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen

import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
COLS = ['event', 'sat', 'id', 'seconds', 'lat', 'lon', 'joules', 'quality']


def extract(cache, output):
    manifest = json.loads((HERE / 'sources.json').read_text())
    with output.open('w', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(COLS)
        for row in manifest['granules']:
            path = cache / row['url'].rsplit('/', 1)[1]
            if not path.is_file():
                cache.mkdir(parents=True, exist_ok=True)
                raw = urlopen(row['url'], timeout=90).read()
                if hashlib.sha256(raw).hexdigest() != row['sha256']:
                    raise ValueError('Source bytes changed: ' + row['url'])
                path.write_bytes(raw)
            if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('Granule checksum mismatch: ' + str(path))
            with netCDF4.Dataset(path) as data:
                lat = data['group_lat'][:]
                lon = data['group_lon'][:]
                # Time in seconds relative to the integer UTC second in CNEOS.
                clock = data['group_time_offset']
                epoch = datetime.fromisoformat(clock.units.split(' since ')[1])
                times = (epoch - datetime.fromisoformat(row['event'])).total_seconds() + clock[:].astype(float)
                mask = (abs(lat - row['lat']) <= 1.8) & (abs(lon - row['lon']) <= 2.4) & (abs(times) <= 3.)
                ids = data['group_id'][:]
                energies = data['group_energy'][:]
                flags = data['group_quality_flag'][:]
                for i in np.flatnonzero(mask):
                    writer.writerow([row['event'], row['sat'], int(ids[i]), f'{times[i]:.8f}',
                                     f'{lat[i]:.7f}', f'{lon[i]:.7f}', f'{energies[i]:.10e}', int(flags[i])])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE / 'groups.csv')
    args = parser.parse_args()
    extract(args.cache, args.output)

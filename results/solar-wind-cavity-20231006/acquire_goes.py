"""Retrieve GOES-R Level 2 one-minute data and store scalar field and quality."""
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from netCDF4 import Dataset, num2date

ROOT = Path(__file__).resolve().parent


def build(reader):
    sources = {}
    rows = {}
    for satellite in [16, 18]:
        filename = f'dn_magn-l2-avg1m_g{satellite}_d20231006_v2-0-2.nc'
        url = (f'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{satellite}/'
               f'l2/data/magn-l2-avg1m/2023/10/{filename}')
        raw = reader(satellite, url)
        sources[str(satellite)] = {'url': url, 'sha256': hashlib.sha256(raw).hexdigest()}
        with Dataset(filename, memory=raw) as dataset:
            times = num2date(dataset['time'][:], units=dataset['time'].units)
            variables = ['b_total', 'b_quality', 'DQF']
            values = [dataset[key][:] for key in variables]
            longitudes = dataset['orbit_llr_geo'][:, 1]
            for t, field, quality, dqf, lon in zip(times, *values, longitudes):
                if 20 <= t.hour <= 23:
                    rows.setdefault(t.strftime('%Y-%m-%dT%H:%MZ'), []).extend(
                        [round(float(field), 4), int(quality), int(dqf), round(float(lon), 3)])
    with (ROOT/'goes.csv').open('w', newline='') as handle:
        cols = [f'g{sat}.{key}' for sat in [16,18] for key in ['B','quality','DQF','lon_deg']]
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(['time_utc', *cols])
        writer.writerows([time, *row] for time,row in sorted(rows.items()))
    (ROOT/'goes_sources.json').write_text(json.dumps(sources, indent=2)+'\n')


if __name__ == '__main__':
    def fetch(satellite, url):
        with urlopen(url, timeout=120) as response:
            return response.read()
    build(fetch)

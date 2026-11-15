"""Rebuild the aligned one-minute medians from the instrument scalar products."""
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCES = json.loads((ROOT / 'sources.json').read_text())
START = '2023-10-06T12:00:00Z'
STOP = '2023-10-07T06:00:00Z'


def reduce_product(name, spec, raw):
    columns = [col for field in spec['fields'] for col in field['columns']]
    frame = pd.read_csv(io.BytesIO(raw), names=['time', *columns])
    frame.index = pd.to_datetime(frame.pop('time'), utc=True)
    frame = frame[(frame.index >= pd.Timestamp(START)) & (frame.index < pd.Timestamp(STOP))]
    if frame.empty:
        raise ValueError(f'No observations in requested date range for {name}')
    for field in spec['fields']:
        cols = field['columns']
        frame[cols] = frame[cols].replace(field['fill'], float('nan'))
    if name == 'dscovr':
        frame = frame[frame['flag'] == 0].drop(columns='flag')
    return frame.resample('1min').median().add_prefix(name + '.')


def build(raw_reader):
    frames = []
    provenance = {}
    for name, spec in SOURCES.items():
        query = {'id': spec['id'], 'parameters': ','.join(f['parameter'] for f in spec['fields']),
                 'time.min': START, 'time.max': STOP, 'format': 'csv'}
        url = 'https://cdaweb.gsfc.nasa.gov/hapi/data?' + urlencode(query)
        raw = raw_reader(name, url)
        provenance[name] = {'url': url, 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
        frames.append(reduce_product(name, spec, raw))
    output = pd.concat(frames, axis=1)
    output.index.name = 'time_utc'
    text = output.to_csv(date_format='%Y-%m-%dT%H:%MZ', float_format='%.8g')
    (ROOT / 'observations.csv.gz').write_bytes(gzip.compress(text.encode(), mtime=0))
    (ROOT / 'retrieval.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    def fetch(name, url):
        with urlopen(url, timeout=120) as response:
            return response.read()
    build(fetch)

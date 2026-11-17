"""Measure known stellar pulsations away from the isolated short-cadence trigger windows."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

from wd_search import flatten


def measure(visit, period_limits, cache):
    path = cache / visit['filename']
    if not path.is_file():
        with urlopen(visit['url'], timeout=150) as response:
            path.write_bytes(response.read())
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha != visit['sha256']:
        raise ValueError('Light curve hash changed: ' + path.name)
    with fits.open(path) as hdus:
        data = hdus[1].data
        good = (data['QUALITY'] == 0) & np.isfinite(data['TIME'] + data['SAP_FLUX'])
        t, f = data['TIME'][good], data['SAP_FLUX'][good].astype(float)
    trend = flatten(t, f)
    good = np.isfinite(trend)
    for event in visit['isolated']:
        # Exclude the same central windows and the shoulders used by the detector.
        good &= abs(t-event['epoch_btjd']) > 1.8*event['duration_min']/1440
    t, f = t[good], f[good]/trend[good]-1
    # Also reject extreme unflagged points from the stellar-frequency fit.
    median = np.median(f)
    keep = abs(f-median) <= 8*1.4826*np.median(abs(f-median))
    t, f = t[keep], f[keep]
    ls = LombScargle(t-t[0], f)
    freq, power = ls.autopower(minimum_frequency=86400./period_limits[1],
                               maximum_frequency=86400./period_limits[0], samples_per_peak=5)
    best = freq[np.argmax(power)]
    _, sine, cosine = ls.model_parameters(best)
    return dict(tic=visit['tic'], sector=visit['sector'], cadence_s=visit['cadence_s'],
                file=path.name, url=visit['url'], sha256=sha, fit_points=len(t),
                search_period_s=period_limits, strongest_period_s=float(86400./best),
                semi_amplitude=float(np.hypot(sine, cosine)), isolated_boxes=len(visit['isolated']))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('scan', type=Path, help='directory of fast_roster.py product JSONs, or their gzip JSON array')
    p.add_argument('cache', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    if args.scan.is_dir():
        records = [json.loads(f.read_text()) for f in args.scan.glob('*.json')]
    else:
        with gzip.open(args.scan, 'rt') as stream:
            records = json.load(stream)
    args.cache.mkdir(exist_ok=True, parents=True)
    # Ranges surround the known pulsation modes of each target (see SHORT_LITERATURE.md).
    bands = {'8988748': (700, 1400), '199666369': (100, 120)}
    selected = [r for r in records if r['tic'] in bands and r.get('isolated') and r['cadence_s']==20]
    out = [measure(r, bands[r['tic']], args.cache) for r in sorted(selected,key=lambda r:(r['tic'],r['sector']))]
    args.output.write_text(json.dumps(out, indent=2)+'\n')

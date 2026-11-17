"""Check three isolated PRIME peaks for known moving objects and archival imaging."""
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from crossmatch import coords
from fetch import OUT, write_csv
from imaging import search
from independent import skybot
from prime_candidates import NAMES, CACHE


def moving(task):
    name, ra, dec, peak = task
    data = np.load(Path('.survey-cache/curves') / (name + '.npz'))['data']
    valid = np.isfinite(data).all(axis=1) & (data[:, 2] > 0)
    times = data[valid, 0]
    jd = round(float(times[np.argmin(abs(times - peak))]), 6)
    # K95 is a neighboring instrument at SAAO Sutherland. Use a wide 120"
    # radius; neither a stationary centroid nor uncataloged objects are tested.
    r = skybot(name, ra, dec, jd, 'K95', CACHE / (name + '-skybot.xml'))
    r.update(observed_jd=jd, model_jd=peak, radius_arcsec=120,
             observer='SAAO Sutherland K95 site proxy')
    return r


if __name__ == '__main__':
    rows = list(csv.DictReader((OUT / 'flag_alternatives.csv').open()))
    rows += list(csv.DictReader((OUT / 'calibrated_alternatives.csv').open()))
    t0 = {r['event']: float(r['fs_t0']) for r in rows if r['event'] in NAMES}
    targets = [(r['event'], *coords(r['ra'], r['dec']), t0[r['event']])
               for r in csv.DictReader((OUT / 'index.csv').open()) if r['event'] in t0]
    with ThreadPoolExecutor(3) as pool:
        ephemerides = list(pool.map(moving, targets))
    write_csv(OUT / 'prime_skybot.csv', ephemerides)
    with ThreadPoolExecutor(6) as pool:
        observations = list(pool.map(search, [(n, ra, dec, a) for n, ra, dec, _ in targets
                                              for a in ('MAST', 'IRSA')]))
    write_csv(OUT / 'prime_imaging_queries.csv', [r for r, _ in observations])
    write_csv(OUT / 'prime_imaging_products.csv', [p for _, rows in observations for p in rows])

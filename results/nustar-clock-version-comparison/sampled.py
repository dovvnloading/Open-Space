"""Sample differences between frozen NuSTAR fine-clock solutions."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time


INPUTS = Path(os.environ['NUSTAR_CLOCK_INPUTS'])
OUTPUT = Path(os.environ['NUSTAR_CLOCK_OUTPUT'])
EXPECTED = {
    234: 'fe872bf396be44c6e4ee20551e4fde89c9cff0809e16057b5a35bcc2b0015e8e',
    235: '24a440096719c2f8c992d6c03e3c3ac769bc7f7638826b1a422d64a1a9d75d19',
    236: 'f30744159e4dcff43a525d15318bd4fb88e561c1fda9fcc508e0db752bb30f96',
}
VERSIONS = (234, 235, 236)
PAIRS = ((234, 235), (235, 236), (234, 236))

OUTPUT.mkdir(parents=True, exist_ok=True)
for version, expected_digest in EXPECTED.items():
    path = INPUTS / f'v{version}.fits.gz'
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_digest:
        raise ValueError(f'input digest mismatch {path.name}')

start = Time('2026-10-10T02:00:01', scale='utc')
stop = Time('2026-10-10T07:50:00', scale='utc')
duration = (stop - start).to_value('sec')
seconds = np.r_[np.arange(0, duration, 60), duration]
curves = {}
evidence = {}

for version in VERSIONS:
    path = INPUTS / f'v{version}.fits.gz'
    with fits.open(path, memmap=False) as hdus:
        header = hdus['NU_FINE_CLOCK'].header
        rows = hdus['NU_FINE_CLOCK'].data
        reference = Time(
            header['MJDREFI'], header['MJDREFF'],
            format='mjd', scale=header['TIMESYS'].lower(),
        )
        timestamps = (start.tt - reference).to_value('sec') + seconds
        knots = np.asarray(rows['TIME'], dtype=float)
        offset = np.asarray(rows['CLOCK_OFF_CORR'], dtype=float)
        frequency = np.asarray(rows['CLOCK_FREQ_CORR'], dtype=float)
        indices = np.searchsorted(knots, timestamps, side='right') - 1
        covered = (indices >= 0) & (indices < len(knots) - 1)
        values = np.full(len(timestamps), np.nan)
        left = indices[covered]
        gap = knots[left + 1] - knots[left]
        x = (timestamps[covered] - knots[left]) / gap

        # FITS describes the offset and its derivative at each control point.
        values[covered] = (
            (2 * x**3 - 3 * x**2 + 1) * offset[left]
            + (x**3 - 2 * x**2 + x) * gap * frequency[left]
            + (-2 * x**3 + 3 * x**2) * offset[left + 1]
            + (x**3 - x**2) * gap * frequency[left + 1]
        )
        curves[version] = values
        evidence[str(version)] = {
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'header_tstart': header['TSTART'],
            'header_tstop': header['TSTOP'],
            'min_knot': float(knots.min()),
            'max_knot': float(knots.max()),
            'covered_samples': int(covered.sum()),
            'uncovered_samples': int((~covered).sum()),
            'max_bracket_seconds': float(gap.max()) if len(gap) else None,
            'monotone_knots': bool(np.all(np.diff(knots) > 0)),
        }

result = {
    'duration_s': duration,
    'samples': len(seconds),
    'versions': evidence,
    'pairs': {},
}
for earlier, later in PAIRS:
    difference = (curves[later] - curves[earlier]) * 1e6
    difference = difference[np.isfinite(difference)]
    result['pairs'][f'{later}-{earlier}'] = {
        'common_samples': len(difference),
        'min_us': float(difference.min()) if len(difference) else None,
        'max_us': float(difference.max()) if len(difference) else None,
        'max_abs_us': float(np.max(abs(difference))) if len(difference) else None,
        'max_abs_change_from_first_us': (
            float(np.max(abs(difference - difference[0]))) if len(difference) else None
        ),
    }

(OUTPUT / 'sampled.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))

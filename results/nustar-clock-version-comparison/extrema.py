"""Compare three published NuSTAR clock tables across a fixed scheduling envelope."""

import hashlib
import json
import os
from pathlib import Path

from astropy.io import fits
from astropy.time import Time
import astropy
import numpy as np
import scipy
from scipy.interpolate import CubicHermiteSpline


SOURCE = Path(os.environ['NUSTAR_CLOCK_INPUTS'])
RUN = Path(os.environ['NUSTAR_CLOCK_OUTPUT'])
RUN.mkdir(parents=True, exist_ok=True)
START, END = '2026-10-10T02:00:01', '2026-10-10T07:50:00'
HASHES = {
    234: 'fe872bf396be44c6e4ee20551e4fde89c9cff0809e16057b5a35bcc2b0015e8e',
    235: '24a440096719c2f8c992d6c03e3c3ac769bc7f7638826b1a422d64a1a9d75d19',
    236: 'f30744159e4dcff43a525d15318bd4fb88e561c1fda9fcc508e0db752bb30f96',
}


def audit():
    model = {}
    coverage = {}
    epochs = Time([START, END], scale='utc').tt
    clock_reference = None
    for version, digest in HASHES.items():
        path = SOURCE / f'v{version}.fits.gz'
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'input digest mismatch {path.name}')
        with fits.open(path) as hdus:
            hdu = hdus['NU_FINE_CLOCK']
            header = hdu.header
            key = tuple(header[name] for name in ('TIMESYS', 'TIMEUNIT', 'MJDREFI', 'MJDREFF'))
            if key[:2] != ('TT', 's') or (clock_reference is not None and key != clock_reference):
                raise ValueError('inconsistent time reference')
            clock_reference = key
            start, end = (epochs - Time(*key[2:], scale='tt', format='mjd')).to_value('sec')
            if start < header['TSTART'] or end > header['TSTOP']:
                raise ValueError('observation interval outside stated clock validity')
            time = np.asarray(hdu.data['TIME'], dtype=float)
            if np.any(np.diff(time) <= 0) or start < time[0] or end > time[-1]:
                raise ValueError('clock knots do not bracket the interval')
            first = np.searchsorted(time, start, side='right') - 1
            last = np.searchsorted(time, end, side='left') + 1
            selection = slice(first, last)
            knots = time[selection].copy() - start
            offset = np.asarray(hdu.data['CLOCK_OFF_CORR'][selection], dtype=float)
            frequency = np.asarray(hdu.data['CLOCK_FREQ_CORR'][selection], dtype=float)
            error = np.asarray(hdu.data['CLOCK_ERR_CORR'][selection], dtype=float)
            model[version] = CubicHermiteSpline(knots, offset, frequency, extrapolate=False)
            coverage[version] = {'bracketing_knots': len(knots), 'max_gap_seconds': float(np.diff(knots).max()),
                                 'clock_error_column_range': [float(error.min()), float(error.max())]}
    span = float(end - start)
    result = {'scope': 'Published clock model differences over the as-flown schedule, including gaps. No photon or phase measurement.',
              'utc_envelope': [START, END], 'input_sha256': HASHES, 'coverage': coverage,
              'runtime': {'numpy': np.__version__, 'scipy': scipy.__version__, 'astropy': astropy.__version__},
              'pairs': {}}
    for old, new in ((234, 235), (235, 236), (234, 236)):
        a, b = model[old], model[new]
        breaks = np.unique(np.clip(np.r_[0., span, a.x, b.x], 0., span))
        # Each difference between successive union knots is itself a cubic.
        difference = CubicHermiteSpline(breaks, b(breaks) - a(breaks), b(breaks, 1) - a(breaks, 1),
                                        extrapolate=False)
        roots = difference.derivative().roots(extrapolate=False)
        points = np.r_[breaks, roots[np.isfinite(roots)]]
        microseconds = 1e6 * difference(points)
        lower, upper = int(microseconds.argmin()), int(microseconds.argmax())
        result['pairs'][f'{new}-{old}'] = {
            'union_segments': len(breaks) - 1, 'finite_derivative_roots': int(np.isfinite(roots).sum()),
            'min_us': float(microseconds[lower]), 'max_us': float(microseconds[upper]),
            'envelope_offset_at_min_seconds': float(points[lower]),
            'envelope_offset_at_max_seconds': float(points[upper]),
            'max_abs_change_from_start_us': float(np.max(np.abs(microseconds - microseconds[0]))),
        }
    result['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    with (RUN / 'extrema.json').open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(audit(), stream, indent=2, allow_nan=False)
        stream.write('\n')

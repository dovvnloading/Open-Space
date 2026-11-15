"""Profile J0720 eclipse curvature in the two frozen DR24 light curves."""
import argparse
import math
from pathlib import Path

import numpy as np
from astropy import units as units
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from scipy.optimize import brentq, least_squares, minimize_scalar

from reproduce import chirp_rate, DAY, YEAR, render_table

HERE = Path(__file__).resolve().parent
PERIOD = .0314116865 * DAY
EPOCH = 2458796.922867  # BJD_TDB


def integrated_eclipse(seconds, inner, outer, exposure):
    """Boxcar-averaged symmetric trapezoid, unit central depth."""
    def ramp(x):
        return np.maximum(x, 0)**2/(2*(outer-inner))

    def primitive(x):
        return ramp(x+outer) - ramp(x+inner) - ramp(x-inner) + ramp(x-outer)

    return (primitive(seconds+exposure/2) - primitive(seconds-exposure/2))/exposure


def fit_curvature(phase, years, red, mag, error, exposure, predicted):
    def residuals(params, curvature):
        inner, ingress, g0, r0, dg, dr, epoch_offset, slope = params
        eclipse = integrated_eclipse(phase - epoch_offset - slope*years - curvature*years**2,
                                     inner, inner+ingress, exposure)
        depth = dg+(dr-dg)*red
        model = g0+(r0-g0)*red-2.5*np.log10(1-depth*eclipse)
        return (mag-model)/np.hypot(error, .025)

    start = [10, 35, 18.92, 19.21, .48, .48, 0, 0]
    low = [0, 3, 18.6, 18.9, .04, .04, -40, -30]
    high = [35, 65, 19.3, 19.6, .95, .95, 40, 30]

    def profile(curvature):
        fitted = least_squares(residuals, start, args=(curvature,), bounds=(low, high))
        return float(fitted.fun @ fitted.fun)

    result = minimize_scalar(profile, bounds=(-2,2), method='bounded')
    if not result.success:
        raise RuntimeError('Curvature minimization failed')
    def confidence_edge(value):
        return profile(value) - result.fun - 3.841459
    lower = brentq(confidence_edge, -6, result.x)
    upper = brentq(confidence_edge, result.x, 6)
    factor = 2*PERIOD/YEAR**2
    return dict(best_pdot=result.x*factor, low_pdot=lower*factor, high_pdot=upper*factor,
                delta_chi2_zero=profile(0)-result.fun, delta_chi2_GR=profile(predicted)-result.fun)


def calculate():
    data = np.genfromtxt(HERE/'ztf_j0720.csv', names=True, delimiter=',', encoding='utf-8', dtype=None)
    data = data[data['catflags']==0]
    # HMJD is heliocentric UTC at midpoint. Change origin to SSB, and UTC to TDB.
    ra, dec = np.deg2rad([110.0124128, 64.6631933])
    direction = np.array([math.cos(ra)*math.cos(dec), math.sin(ra)*math.cos(dec), math.sin(dec)])
    dates = Time(data['hmjd_utc'], format='mjd', scale='utc')
    sun_xyz = get_body_barycentric('sun', dates.tdb).xyz.to_value(units.au)
    correction_s = (direction @ sun_xyz)*499.00478383615643
    bjd = dates.tdb.jd + correction_s/DAY
    phase = ((bjd-EPOCH)*DAY + PERIOD/2) % PERIOD - PERIOD/2
    mjd = bjd-2400000.5
    years = (mjd-58796)/365.25
    rate = chirp_rate(PERIOD, .31, .274)
    predicted = .5*rate*YEAR**2/PERIOD
    results = []
    for exposure, window, omit_dense in [(30,180,False),(20,180,False),(40,180,False),(30,120,False),(30,180,True)]:
        mask = abs(phase)<window
        if omit_dense:
            mask &= np.floor(mjd)!=58796
        row = dict(exposure_s=exposure, half_window_s=window, omit_mjd58796=int(omit_dense),
                   selected=len(mjd[mask]), nights=len(set(np.floor(mjd[mask]))))
        row.update(fit_curvature(phase[mask], years[mask], (data['band'][mask]=='r').astype(int),
                                 data['mag'][mask], data['magerr'][mask], exposure, predicted))
        results.append(row)
    return render_table(results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = HERE/'ztf_timing.csv'
    calculated = calculate()
    if args.check:
        import io
        output = np.genfromtxt(io.StringIO(calculated), delimiter=',', names=True)
        frozen = np.genfromtxt(path, delimiter=',', names=True)
        for column in output.dtype.names:
            # Float-optimizer tolerances across SciPy/BLAS versions.
            if not np.allclose(output[column], frozen[column], rtol=3e-3, atol=1e-17):
                raise SystemExit(f'ZTF result differs in {column}')
    else:
        path.write_text(calculated)

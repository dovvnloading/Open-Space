"""Quantify the influence of one old exposure on a joint PS1–ZTF curvature fit."""
import argparse
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, get_body_barycentric
from astropy.time import Time
from astropy.utils import iers
from scipy.optimize import least_squares, minimize_scalar

from reproduce import DAY, YEAR, chirp_rate, render_table
from ztf_timing import EPOCH, PERIOD, integrated_eclipse

HERE = Path(__file__).resolve().parent


def calculate():
    iers.conf.auto_download = False
    source = SkyCoord(110.0124128*u.deg, 64.6631933*u.deg)
    z = np.genfromtxt(HERE/'ztf_j0720.csv', names=True, delimiter=',', dtype=None, encoding='utf8')
    z = z[z['catflags']==0]
    time = Time(z['hmjd_utc'], format='mjd', scale='utc')
    sun = get_body_barycentric('sun', time.tdb).xyz.to_value(u.au)
    z_bjd = time.tdb.jd + (source.cartesian.xyz.value @ sun)*499.00478383615643/DAY
    ps = np.genfromtxt(HERE/'ps1_j0720.csv', names=True, delimiter=',', dtype=None, encoding='utf8')
    ps = ps[ps['filterID']<=2]
    time = Time(ps['obsTime'],format='mjd',scale='utc', location=EarthLocation.from_geodetic(-156.2572,20.7074,3055))
    ps_bjd = (time.tdb + time.light_travel_time(source)).jd
    mjd = np.r_[z_bjd, ps_bjd]-2400000.5
    phase = ((mjd-(EPOCH-2400000.5))*DAY + PERIOD/2)%PERIOD-PERIOD/2
    years = (mjd-58796)/365.25
    mag = np.r_[z['mag'], -2.5*np.log10(ps['psfFlux']/3631)]
    error = np.r_[z['magerr'], 1.0857*ps['psfFluxErr']/ps['psfFlux']]
    band = np.r_[np.where(z['band']=='r',1,0), 1+ps['filterID']]
    exposure = np.r_[np.full(len(z),30.),ps['expTime']]
    ps_baseline = [np.median(mag[(band==b)&(abs(phase)>180)]) for b in (2,3)]
    predicted = .5*chirp_rate(PERIOD,.31,.274)*YEAR**2/PERIOD
    results = []
    for omitted in (0,55589,55918,55955,56293):
        keep = (abs(phase)<180)&(np.floor(mjd)!=omitted)
        selected_band = band[keep]
        def profile(curvature):
            def residuals(x):
                inner, ingress = x[:2]
                baseline, depth = x[2:6], x[6:10]
                shift = x[10] + x[11]*years[keep] + curvature*years[keep]**2
                profile_flux = integrated_eclipse(phase[keep]-shift, inner, inner+ingress, exposure[keep])
                model = baseline[selected_band]-2.5*np.log10(1-depth[selected_band]*profile_flux)
                # Each PS1 baseline is anchored to off-window observations, not to the dips.
                return np.r_[(mag[keep]-model)/np.hypot(error[keep],.025), (baseline[2:]-ps_baseline)/.015]
            start = [13,33,18.9,19.20,18.91,19.15,.39,.39,.39,.39,0,0]
            low = [0,3,18.4,18.8,18.4,18.8,.04,.04,.05,.05,-40,-30]
            high = [35,65,19.3,19.6,19.4,19.7,.9,.9,.9,.9,40,30]
            fit = least_squares(residuals,start,bounds=(low,high))
            return float(fit.fun @ fit.fun)
        result = minimize_scalar(profile,bounds=(-.65,.65),method='bounded')
        results.append(dict(omit_mjd=omitted,selected=int(keep.sum()),
                            best_pdot=round(result.x*2*PERIOD/YEAR**2,14),
                            delta_chi2_zero=round(profile(0)-result.fun,1),
                            delta_chi2_GR=round(profile(predicted)-result.fun,1)))
    return render_table(results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    dest = HERE/'joint_stress.csv'
    text = calculate()
    if args.check:
        if dest.read_text()!=text:
            raise SystemExit('Joint timing sensitivity changed')
    else:
        dest.write_text(text)

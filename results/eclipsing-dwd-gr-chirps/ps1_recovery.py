"""Recover independent PS1 primary-eclipse exposures using the published linear clock."""
import argparse
import math
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

from reproduce import render_table

HERE = Path(__file__).resolve().parent


def calculate():
    iers.conf.auto_download = False
    rows = np.genfromtxt(HERE/'ps1_j0720.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
    source = SkyCoord(110.0124128*u.deg, 64.6631933*u.deg)
    site = EarthLocation.from_geodetic(-156.2572, 20.7074, 3055)
    times = Time(rows['obsTime'], format='mjd', scale='utc', location=site)
    bjd = (times.tdb + times.light_travel_time(source)).jd
    period = .0314116865 * 86400
    phase_seconds = ((bjd-2458796.922867)*86400+period/2) % period - period/2
    output = []
    for band_id, band in [(1,'g'), (2,'r')]:
        mask = rows['filterID']==band_id
        selected = mask & (abs(phase_seconds)<35)
        psf_base = np.median(rows['psfFlux'][mask])
        ap_base = np.median(rows['apFlux'][mask])
        for i in np.flatnonzero(selected):
            r = rows[i]
            rank = 1 + np.count_nonzero(rows['psfFlux'][mask] < r['psfFlux'])
            output.append(dict(date_utc=times[i].isot, detection_id=int(r['detectID']), band=band,
                               observations_in_band=int(mask.sum()), faintness_rank=rank,
                               offset_from_linear_eclipse_s=round(phase_seconds[i],2),
                               integration_s=float(r['expTime']),
                               psf_deficit_mag=round(-2.5*math.log10(r['psfFlux']/psf_base),4),
                               formal_psf_error_mag=round(2.5*r['psfFluxErr']/(math.log(10)*r['psfFlux']),4),
                               aperture_deficit_mag=round(-2.5*math.log10(r['apFlux']/ap_base),4)))
    return render_table(sorted(output, key=lambda row: row['date_utc']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = HERE/'ps1_recovery.csv'
    output = calculate()
    if args.check:
        if path.read_text()!=output:
            raise SystemExit('PS1 recovery differs from frozen output')
    else:
        path.write_text(output)

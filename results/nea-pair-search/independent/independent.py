"""Refit VE138 with ASSIST/REBOUND and separate topocentric astrometric code."""
import argparse
import csv
import json
from pathlib import Path

import assist
import numpy as np
import rebound
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.utils import iers
from scipy.linalg import lstsq

AU_KM = 149597870.7
C_AU_DAY = 299792.458 * 86400.0 / AU_KM
ARCSEC = 180 * 3600 / np.pi
EPOCH = 2461200.5
SCALE = np.array([1e-4] * 3 + [1e-6] * 3 + [1e-8])
OBLIQUITY = np.radians(84381.448 / 3600.0)
ROTATION = np.array([[1, 0, 0], [0, np.cos(OBLIQUITY), -np.sin(OBLIQUITY)],
                     [0, np.sin(OBLIQUITY), np.cos(OBLIQUITY)]])
iers.conf.auto_download = False


def positions(ephemeris, rows, stations_file):
    utc = Time([row['utc'].removesuffix('Z') for row in rows], format='isot', scale='utc')
    dates = utc.tdb.jd
    geo = np.array([ephemeris.get_particle('Earth', jd - ephemeris.jd_ref).xyz for jd in dates])
    sites = np.array([row['site'] for row in rows])
    for line in stations_file.read_text().splitlines():
        mask = sites == line[:3]
        if not np.any(mask):
            continue
        lon = np.radians(float(line[4:13]))
        rc = 6378137 * float(line[13:21])
        rs = 6378137 * float(line[21:30])
        location = EarthLocation.from_geocentric(rc*np.cos(lon), rc*np.sin(lon), rs, unit='m')
        geo[mask] += location.get_gcrs_posvel(utc[mask])[0].xyz.to_value('au').T
    ra, dec = np.radians([[float(row[key]) for row in rows] for key in ['ra_deg', 'dec_deg']])
    line_of_sight = np.column_stack((np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)))
    east = np.column_stack((-np.sin(ra), np.cos(ra), np.zeros(len(ra))))
    return dates, geo, east, np.cross(line_of_sight, east)


class Astrometry:
    def __init__(self, kernel, directory, bias_archive=None):
        self.eph = assist.Ephem(str(kernel))
        with (directory / 'astrometry.csv').open() as handle:
            rows = sorted((r for r in csv.DictReader(handle) if r['designation'] in
                           ('2007 VE138', '2026 VV15')), key=lambda r: r['utc'])
        self.corrections = []
        if bias_archive:
            from catalog_bias import corrected
            rows, self.corrections = corrected(rows, bias_archive)
        self.dates, self.observer, self.east, self.north = positions(self.eph, rows, directory / 'stations.txt')
        self.sites = np.array([r['site'] for r in rows])
        self.sun = self.eph.get_particle('Sun', EPOCH - self.eph.jd_ref)

    def residuals(self, state):
        """Independent ICRF orbit integration, observer computation, and first-order light time."""
        result = np.empty((len(self.dates), 2))
        early = np.flatnonzero(self.dates < EPOCH)[::-1]
        late = np.flatnonzero(self.dates >= EPOCH)
        for indices in (early, late):
            sim = rebound.Simulation()
            sim.t = EPOCH - self.eph.jd_ref
            xyz = ROTATION @ state[:3] + self.sun.xyz
            vel = ROTATION @ state[3:6] + self.sun.vxyz
            sim.add(x=xyz[0], y=xyz[1], z=xyz[2], vx=vel[0], vy=vel[1], vz=vel[2])
            extras = assist.Extras(sim, self.eph)
            extras.forces = ['SUN', 'PLANETS', 'NON_GRAVITATIONAL', 'EARTH_HARMONICS',
                             'SUN_HARMONICS', 'GR_EIH']
            # Normalized Marsden law g(r)=r^-2, positive A1 pointing away from the Sun.
            extras.alpha, extras.r0, extras.nm, extras.nk = 1., 1., 2., 0.
            acceleration = np.array([state[6], 0., 0.], dtype=np.float64)
            extras.particle_params = acceleration
            for k in indices:
                extras.integrate_or_interpolate(self.dates[k] - self.eph.jd_ref)
                particle = sim.particles[0]
                delta = np.array(particle.xyz) - self.observer[k]
                velocity = np.array(particle.vxyz)
                emitted = delta.copy()
                for _ in range(3):
                    emitted = delta - velocity * np.linalg.norm(emitted) / C_AU_DAY
                unit = emitted / np.linalg.norm(emitted)
                result[k] = [unit @ self.east[k], unit @ self.north[k]]
        return result * ARCSEC


def fit(model, seed, allow_radial):
    active = np.arange(7 if allow_radial else 6)
    step = SCALE[active] * 0.002
    state = np.array(seed)
    for _ in range(6):
        r = model.residuals(state).ravel()
        derivatives = []
        for j, h in zip(active, step):
            direction = np.zeros(7)
            direction[j] = h
            high = model.residuals(state + direction)
            low = model.residuals(state - direction)
            derivatives.append(((high - low) / .004).ravel())
        jacobian = np.column_stack(derivatives)
        change = lstsq(jacobian, -r)[0] * SCALE[active]
        state[active] += change
        if np.max(np.abs(change / SCALE[active])) < 1e-5:
            break
    residual = model.residuals(state)
    return {'heliocentric_ecliptic_j2000_state_and_A1': state.tolist(),
            'positions': len(residual),
            'radial_rms_arcsec': float(np.sqrt(np.sum(residual**2, axis=1).mean())),
            'radial_median_arcsec': float(np.median(np.linalg.norm(residual, axis=1)))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('de440', type=Path)
    parser.add_argument('orbits_dir', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--bias-archive', type=Path)
    args = parser.parse_args()
    seed = json.loads((args.orbits_dir / 'fits/VE138.json').read_text())['full']['state'] + [0.]
    model = Astrometry(args.de440, args.orbits_dir, args.bias_archive)
    out = {name: fit(model, seed, allow) for name, allow in [('gravity', False), ('radial', True)]}
    if args.bias_archive:
        out['catalog_corrections'] = model.corrections
    args.output.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({name: [v['radial_rms_arcsec'], v['heliocentric_ecliptic_j2000_state_and_A1'][-1]] for name, v in out.items() if name in ('gravity', 'radial')}))

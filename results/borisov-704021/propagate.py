"""Integrate the independent ESA asteroid orbit to the encounter with Borisov.

Usage: python propagate.py /path/to/de440s.bsp
Requires numpy, scipy and jplephem. The output uses TDB, ecliptic J2000, AU and days.
"""
import argparse
from pathlib import Path
import numpy as np
from jplephem.spk import SPK
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parent
AU = 149597870.7
DAY = 86400.0
GM = 132712440041.279419 / AU**3 * DAY**2
PLANET_GM = np.array([22031.86855, 324858.592, 403503.235502, 42828.375816,
                      126712764.10, 37940584.8418, 5794556.4, 6836527.10058]) / AU**3 * DAY**2


def initial_state():
    lines = (ROOT / 'data/704021.ke1').read_text().splitlines()
    elements = next(row.split()[1:] for row in lines if row.startswith(' KEP'))
    epoch = next(float(row.split()[1]) for row in lines if row.startswith(' MJD'))
    a, e, inc, node, peri, mean = map(float, elements)
    inc, node, peri, mean = np.deg2rad([inc, node, peri, mean])
    anomaly = mean
    for _ in range(12):
        anomaly -= (anomaly - e * np.sin(anomaly) - mean) / (1 - e * np.cos(anomaly))
    ci, si = np.cos(inc), np.sin(inc)
    co, so = np.cos(node), np.sin(node)
    cw, sw = np.cos(peri), np.sin(peri)
    p = np.array([cw * co - sw * ci * so, cw * so + sw * ci * co, sw * si])
    q = np.array([-sw * co - cw * ci * so, -sw * so + cw * ci * co, cw * si])
    cs, sn, beta = np.cos(anomaly), np.sin(anomaly), np.sqrt(1 - e**2)
    pos = a * ((cs - e) * p + beta * sn * q)
    vel = np.sqrt(GM / a) / (1 - e * cs) * (-sn * p + beta * cs * q)
    return epoch, np.r_[pos, vel]


def integrate(kernel):
    epoch, state = initial_state()
    times = np.genfromtxt(ROOT / 'data/jpl-vectors.csv', delimiter=',', names=True)['mjd_tdb']
    days = np.arange(np.floor(times[0]) - 1, epoch + 2, 0.5)
    jd = days + 2400000.5
    obl = np.deg2rad(23.4392911111)
    rot = np.array([[1, 0, 0], [0, np.cos(obl), -np.sin(obl)], [0, np.sin(obl), np.cos(obl)]])
    with SPK.open(str(kernel)) as eph:
        sun = eph[0, 10].compute(jd)
        planets = np.stack([(eph[0, j].compute(jd) - sun).T @ rot / AU for j in range(1, 9)], axis=1)
    positions = CubicSpline(days, planets, axis=0)

    def derivative(t, y):
        p = positions(t)
        x, v = y[:3], y[3:]
        r = np.linalg.norm(x)
        d = p - x
        planets = np.sum(PLANET_GM[:, None] * (d / np.linalg.norm(d, axis=1)[:, None]**3
                         - p / np.linalg.norm(p, axis=1)[:, None]**3), axis=0)
        relativistic = GM / (173.144632674240**2 * r**3) * ((4 * GM / r - v @ v) * x + 4 * (x @ v) * v)
        return np.r_[v, -GM * x / r**3 + planets + relativistic]

    fit = solve_ivp(derivative, (epoch, times[0] - 0.1), state, rtol=2e-13, atol=8e-15,
                    method='DOP853', dense_output=True, max_step=4)
    if not fit.success:
        raise RuntimeError(fit.message)
    header = 'mjd_tdb,x_au,y_au,z_au,vx_au_day,vy_au_day,vz_au_day'
    np.savetxt(ROOT / 'data/esa-propagation.csv', np.c_[times, fit.sol(times).T], delimiter=',',
               header=header, comments='', fmt='%.15g')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernel', type=Path)
    integrate(parser.parse_args().kernel)

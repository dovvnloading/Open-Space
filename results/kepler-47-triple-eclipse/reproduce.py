"""Reintegrate Kepler-47, validate Orosz d timings, and test the 2029 shadow."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parent
EPOCH = 2454965.0
CENTER = 7357.22  # days after EPOCH
G = 0.01720209895**2  # Orosz et al. 2019
C = 173.144632674  # AU/day


def orbit(end=CENTER + 1, rtol=3e-12):
    rows = np.loadtxt(ROOT / 'initial-state.csv', delimiter=',', skiprows=1, usecols=range(1, 9))
    masses, radii = rows[:, 0], rows[:, -1]
    initial = np.concatenate((rows[:, 1:4], rows[:, 4:7])).ravel()
    total = masses[0] + masses[1]
    eta = masses[0] * masses[1] / total**2

    def derivative(_, state):
        pos, vel = state.reshape(2, 5, 3)
        dx = pos[None, :, :] - pos[:, None, :]
        square = np.einsum('ijk,ijk->ij', dx, dx)
        np.fill_diagonal(square, np.inf)
        acc = G * np.einsum('ij,j,ijk->ik', square**-1.5, masses, dx)
        # Inner binary point-mass first post-Newtonian relative acceleration.
        vector, velocity = dx[0, 1], vel[1] - vel[0]
        r = np.sqrt(square[0, 1])
        rv = vector @ velocity
        radial = (4 + 2*eta)*G*total/r - (1 + 3*eta)*(velocity @ velocity) + 1.5*eta*(rv/r)**2
        correction = G*total/(C*C*r**3) * (radial*vector + (4 - 2*eta)*rv*velocity)
        acc[0] -= masses[1]/total * correction
        acc[1] += masses[0]/total * correction
        return np.concatenate((vel, acc)).ravel()

    grid = np.arange(0, end + .12, .12)
    sol = solve_ivp(derivative, (0, grid[-1]), initial, method='DOP853', rtol=rtol,
                    atol=rtol*1e-3, t_eval=grid)
    if not sol.success:
        raise RuntimeError(sol.message)
    pos, vel = sol.y.T.reshape(-1, 2, 5, 3).transpose(1, 0, 2, 3)
    spline = CubicHermiteSpline(grid, pos, vel)
    speed = spline.derivative()

    def apparent(time):
        xyz = spline(time)
        return xyz + speed(time) * xyz[..., 2:3]/C  # retarded emission to first order in v/c

    return apparent, radii


def separation(pos, first, second):
    return np.linalg.norm(pos[..., first, :2] - pos[..., second, :2], axis=-1)


def closest(apparent, first, second, low, high):
    return minimize_scalar(lambda t: separation(apparent(t), first, second),
                           method='bounded', bounds=(low, high),
                           options={'xatol': 1e-9}).x


def model_checks(apparent, radii):
    rows = []
    with (ROOT / 'published-d-transits.csv').open() as fh:
        for row in csv.DictReader(fh):
            reference = float(row['bjd_minus_2455000']) + 35
            t = closest(apparent, 0, 3, reference - .35, reference + .35)
            rows.append({'source': row['source'], 'reference_bjd': EPOCH + reference,
                         'calculated_minus_published_minutes': round((t-reference)*1440, 3),
                         'reported_sigma_minutes': round(float(row['time_sigma_d'])*1440, 2),
                         'calculated_impact': round(separation(apparent(t), 0, 3)/radii[0], 5)})
    return rows


def sensitivity(apparent, radii):
    times = CENTER + np.linspace(-.27, .27, 3101)
    stars = apparent(times)/radii[0]
    r = radii/radii[0]
    def margins(shift, height):
        planet = apparent(times + shift)[:, 3, :2]/radii[0]
        planet[:, 1] += height
        separation_a = np.linalg.norm(planet - stars[:, 0, :2], axis=1)
        separation_b = np.linalg.norm(planet - stars[:, 1, :2], axis=1)
        # When the planet is entirely inside A and overlaps B, a luminous
        # part of A reappears relative to the sum of two isolated transits.
        margin_inside_a = 1 - r[3] - separation_a
        partial = np.max(np.minimum(margin_inside_a, r[1] + r[3] - separation_b))
        full = np.max(np.minimum(margin_inside_a, r[1] - r[3] - separation_b))
        return partial, full
    results = []
    for height in [-.10, -.05, 0, .05, .10]:
        grid = np.linspace(-.5, .5, 61)
        values = np.array([margins(offset, height) for offset in grid])
        result = {'cross_track_primary_radii': height}
        for i, label in enumerate(['any_overlap', 'full_containment']):
            bounds = []
            for k in np.flatnonzero(np.diff(values[:, i] > 0)):
                root = brentq(lambda offset: margins(offset, height)[i], grid[k], grid[k+1])
                bounds.append(round(24*root, 3))
            result[label + '_phase_limits_hr'] = bounds
        results.append(result)
    return results


def forecast(apparent, radii):
    # Here B is in front of A and d is in front of both.
    t = closest(apparent, 1, 3, CENTER-.06, CENTER+.06)
    pos = apparent(t)
    assert pos[3, 2] > pos[1, 2] > pos[0, 2]
    def contacts(i, j, radius, left, middle, right):
        def residual(time):
            return separation(apparent(time), i, j) - radius
        return [brentq(residual, left, middle), brentq(residual, middle, right)]
    outer = contacts(1, 3, radii[1]+radii[3], CENTER-.18, t, CENTER+.18)
    inner = contacts(1, 3, radii[1]-radii[3], outer[0], t, outer[1])
    # Check if the entire planet is on the primary throughout B–d contact.
    sample = np.linspace(*outer, 1001)
    cap = (separation(apparent(sample), 0, 3) + radii[3]) / radii[0]
    def triple_margin(time):
        # The earliest intersection occurs just after outer B–d contact. Find
        # the point of their common lens closest to the center of A.
        xyz = apparent(time)[:, :2]/radii[0]
        b, d = xyz[1]-xyz[0], xyz[3]-xyz[0]
        rb, rd = radii[[1, 3]]/radii[0]
        sep = np.linalg.norm(d-b)
        if sep >= rb+rd:
            return max(sep-rb-rd, 1e-10)
        pts = []
        for center, radius, other, other_radius in [(b, rb, d, rd), (d, rd, b, rb)]:
            near = center * max(0, 1-radius/np.linalg.norm(center))
            if np.linalg.norm(near-other) <= other_radius:
                pts.append(near)
        if abs(rb-rd) <= sep:
            ux = (d-b)/sep
            chord = (rb*rb-rd*rd+sep*sep)/(2*sep)
            offset = np.sqrt(max(0, rb*rb-chord*chord))*np.array([-ux[1], ux[0]])
            pts.extend([b+chord*ux+offset, b+chord*ux-offset])
        return min(np.linalg.norm(v) for v in pts) - 1

    triple_start = brentq(triple_margin, outer[0], inner[0])
    assert max((separation(apparent(np.linspace(*inner, 201)), 0, 3)+radii[3])/radii[0]) < 1
    # Integrate primary surface brightness over the planet disk at closest B–d approach.
    rho = np.sqrt((np.arange(80) + .5)/80)
    theta = 2*np.pi*(np.arange(360)+.5)/360
    x,y = (pos[3, :2]-pos[0, :2])/radii[0]
    px = x + radii[3]/radii[0] * rho[:, None] * np.cos(theta)
    py = y + radii[3]/radii[0] * rho[:, None] * np.sin(theta)
    mu = np.sqrt(np.maximum(1-px*px-py*py, 0))
    q1, q2 = .502, .188
    u1 = 2*np.sqrt(q1)*q2
    u2 = np.sqrt(q1)*(1-2*q2)
    intensity = 1-u1*(1-mu)-u2*(1-mu)**2
    intensity[px*px+py*py>=1] = 0.
    recovered = (radii[3]/radii[0])**2 * np.mean(intensity)/(1-u1/3-u2/6)
    # The visible foreground star B loses its own light due to the planet.
    secondary_fraction = .005
    delta = (recovered - secondary_fraction*(radii[3]/radii[1])**2)/(1 + secondary_fraction)
    return {'outer_contact_bjd_tdb': [round(EPOCH+v, 6) for v in outer],
            'inner_contact_bjd_tdb': [round(EPOCH+v, 6) for v in inner],
            'closest_bjd_tdb': round(EPOCH+t, 6),
            'triple_contact_bjd_tdb': [round(EPOCH+triple_start, 6), round(EPOCH+outer[1], 6)],
            'planet_B_distance_over_radius_sum': round(separation(pos, 1, 3)/(radii[1]+radii[3]), 5),
            'planet_A_distance_over_A_radius': round(separation(pos, 0, 3)/radii[0], 5),
            'max_planet_outer_edge_on_A_during_B_contact': round(max(cap), 5),
            'mutual_contact_minutes': round((outer[1]-outer[0])*1440, 2),
            'full_containment_minutes': round((inner[1]-inner[0])*1440, 2),
            'recovery_vs_independent_A_transits_ppm': round(recovered/(1+secondary_fraction)*1e6),
            'net_recovery_including_uniform_B_occultation_ppm': round(delta*1e6)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rtol', type=float, default=3e-12)
    args = parser.parse_args()
    position, radii = orbit(rtol=args.rtol)
    print(json.dumps({'forecast': forecast(position, radii), 'validation': model_checks(position, radii),
                      'stress_tests': sensitivity(position, radii)}, indent=2))

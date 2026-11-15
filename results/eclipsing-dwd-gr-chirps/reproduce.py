"""Recompute the two orbital timing forecasts from the frozen literature table."""
import csv
import io
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
C = 299792458.0  # m/s, exact
GM_SUN = 1.3271244e20  # IAU 2015 nominal solar mass parameter, m^3/s^2
DAY = 86400.0
YEAR = 365.25 * DAY


def chirp_rate(period_s, primary, secondary):
    """Signed circular-orbit quadrupole dP/dt, in seconds per second."""
    mass_factor = primary * secondary / (primary + secondary) ** (1 / 3)
    return -(192 * math.pi / 5) * (math.tau * GM_SUN / (C**3 * period_s)) ** (5 / 3) * mass_factor


def energy_balance_rate(period_s, primary, secondary):
    """Different derivation from Kepler's law, orbital energy and quadrupole power."""
    mu1, mu2 = GM_SUN * primary, GM_SUN * secondary
    a = ((mu1 + mu2) * (period_s / math.tau) ** 2) ** (1 / 3)
    # Both energy and power include a common factor of G, which cancels.
    energy_slope = mu1 * mu2 / (2 * a**2)
    power = 32 * (mu1 * mu2) ** 2 * (mu1 + mu2) / (5 * C**5 * a**5)
    return -(3 / 2) * period_s * power / (a * energy_slope)


def advance(period_s, rate):
    # Magnitude of O-C relative to tangent ephemeris, baseline of ten Julian years.
    return -0.5 * rate * (10 * YEAR) ** 2 / period_s


def semi_amplitude_rate(k1_kms, k2_kms, inclination_deg):
    """Double-lined circular orbit. Substitution removes period and mass-distance scale."""
    k1, k2 = 1000 * k1_kms, 1000 * k2_kms
    return -(192 * math.pi / 5) * k1 * k2 * (k1 + k2) ** 3 / (C * math.sin(math.radians(inclination_deg))) ** 5


def apsidal_rate(period_s, total_mass):
    n = math.tau / period_s
    rate = 3 * n * (n * GM_SUN * total_mass / C**3) ** (2 / 3)
    return math.degrees(rate) * YEAR


def render_table(rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, list(rows[0]), lineterminator='\n')
    writer.writeheader()
    for row in rows:
        writer.writerow({k: (format(v, '.8g') if isinstance(v, float) else v) for k, v in row.items()})
    return stream.getvalue()


def calculate():
    rows, apsides, velocities = [], [], []
    with (HERE / 'parameters.csv').open() as source:
        for raw in csv.DictReader(source):
            name = raw['target']
            x = {k: float(v) for k, v in raw.items() if k != 'target' and v}
            p = x['period_d'] * DAY
            m1, m2 = x['m1_msun'], x['m2_msun']
            pair_low = (m1 - x['m1_minus'], m2 - x['m2_minus'])
            pair_high = (m1 + x['m1_plus'], m2 + x['m2_plus'])
            rate = chirp_rate(p, m1, m2)
            if not math.isclose(rate, energy_balance_rate(p, m1, m2), rel_tol=2e-14):
                raise ArithmeticError('Independent energy-balance calculation disagrees')
            forecast = advance(p, rate)
            rows.append(dict(target=name, signed_pdot=rate, printed_abs_pdot=x['printed_abs_pdot_e15']*1e-15,
                             ratio=-rate/(x['printed_abs_pdot_e15']*1e-15), advance_10yr_s=forecast,
                             mass_corner_min_s=advance(p, chirp_rate(p, *pair_low)),
                             mass_corner_max_s=advance(p, chirp_rate(p, *pair_high)),
                             midpoint_residual_in_10yr_chord_s=forecast/4))
            if 'k2_kms' in x:
                def rv_advance(offset):
                    rv = semi_amplitude_rate(x['k1_kms'] + offset*x['k1_err'],
                                             x['k2_kms'] + offset*x['k2_err'], x['incl_deg'])
                    return advance(p, rv)
                velocities.append(dict(target=name, rv_advance_10yr_s=rv_advance(0),
                                       rv_corner_min_s=rv_advance(-1), rv_corner_max_s=rv_advance(1),
                                       mass_fit_advance_10yr_s=forecast))
            if 'secondary_offset_s' in x:
                rate_deg = apsidal_rate(p, m1 + m2)
                printed_eq4 = 5.447e-4 * ((m1 + m2) / x['period_d'])**(2/3) * 365.25 / x['period_d']
                if not math.isclose(printed_eq4, rate_deg, rel_tol=3e-5):
                    raise ArithmeticError('Equation 4 cross-check disagrees')
                offset = x['secondary_offset_s'] - x['romer_s']
                theta = math.radians(rate_deg * 10)
                # For times 0,T,2T with known first offset and unknown sine coefficient,
                # minimizing max(|D(T)-D(0)|, |D(2T)-D(0)|) is attained at equal,
                # opposite changes when cos(theta)>0. See README for the recurrence.
                cosine = math.cos(theta)
                forward_bound = abs(offset) * 2 * (1 - cosine) / (1 + 2 * cosine)
                apsides.append(dict(target=name, deg_per_year=rate_deg, printed_eq4_deg_per_year=printed_eq4,
                                    mass_corner_min_deg_per_year=apsidal_rate(p, sum(pair_low)),
                                    mass_corner_max_deg_per_year=apsidal_rate(p, sum(pair_high)),
                                    apsidal_period_years=360/rate_deg, eccentric_offset_s=offset,
                                    e_cos_omega=math.pi*offset/(2*p),
                                    change_10yr_if_omega45_s=offset*(math.cos(theta)-math.sin(theta)-1),
                                    min_largest_change_at_10_or_20yr_s=forward_bound,
                                    centered_second_difference_at_10yr_spacing_s=2*offset*(math.cos(theta)-1)))
    control = chirp_rate(765.206543, .26, .50)
    if abs(control + 9.8e-12) > 2.8e-12:
        raise ArithmeticError('J0651 lies outside the measured one-sigma rate')
    j0221 = rows[0]['advance_10yr_s']
    if not math.isclose(j0221, 3.6, abs_tol=.1):
        raise ArithmeticError('J0221 differs from published ten-year forecast')
    checks = dict(missing_normalization=math.tau**(5/3), j0651_calculated_pdot=control,
                  j0651_measured_pdot=-9.8e-12, j0651_measured_1sigma=2.8e-12,
                  j0221_prior_advance_10yr_s=3.6, j0221_prior_1sigma_s=.1,
                  j0221_recalculated_advance_10yr_s=j0221)
    return {'chirps.csv':render_table(rows), 'precession.csv':render_table(apsides),
            'velocity_check.csv':render_table(velocities), 'controls.json':json.dumps(checks, indent=2)+'\n'}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='verify frozen outputs byte for byte')
    args = parser.parse_args()
    for filename, output in calculate().items():
        path = HERE / filename
        if args.check:
            if path.read_text() != output:
                raise SystemExit(f'Output mismatch in {filename}')
        else:
            path.write_text(output)

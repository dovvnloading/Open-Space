"""Propagate the target's Gaia joint covariance locally. Run with Python 3 and NumPy."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PARAMETERS = (
    'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'a_thiele_innes', 'b_thiele_innes',
    'f_thiele_innes', 'g_thiele_innes', 'c_thiele_innes', 'h_thiele_innes',
    'center_of_mass_velocity', 'eccentricity', 'period', 't_periastron',
)
# Tangent-plane position errors are in mas. The means of positions and proper
# motions do not enter the orbit transformation; sample standardized residuals.
AU_KM = 149597870.7
DAYS_PER_YEAR = 365.25
SEED = 96100077
DRAWS = 150000


def read_one(name):
    with (ROOT / name).open() as file:
        row, = csv.DictReader(file)
    return row


def interval(x):
    return dict(zip(('p16', 'median', 'p84'), np.quantile(x, [0.16, 0.50, 0.84]).tolist()))


def analyze(scales):
    row = read_one('data/nss-orbit.csv')
    assert row['bit_index'] == '65535' and row['nss_solution_type'] == 'AstroSpectroSB1'
    coefficients = np.fromstring(row['corr_vec'].strip('()'), sep=',')
    n = len(PARAMETERS)
    assert len(coefficients) == n * (n - 1) // 2
    lower = np.zeros((n, n))
    lower[np.tril_indices(n, -1)] = coefficients
    correlation = np.eye(n) + lower + lower.T
    rng = np.random.default_rng(SEED)
    deviations = rng.standard_normal((DRAWS, n)) @ np.linalg.cholesky(correlation).T
    # Include the exact catalogue point as the first row for reference values.
    deviations = np.vstack([np.zeros(n), deviations])
    samples = {key: float(row[key]) + float(row[key + '_error']) * deviations[:, i]
               for i, key in enumerate(PARAMETERS) if i >= 2}
    A, B, F, G, C, H = (samples[axis + '_thiele_innes'] for axis in 'abfgch')
    # The two singular values of the projected orbital ellipse are a0 and a0*abs(cos i).
    sum_axis = np.hypot(A + G, B - F)
    diff_axis = np.hypot(A - G, B + F)
    a0 = 0.5 * (sum_axis + diff_axis)
    inclination = np.arccos((sum_axis - diff_axis) / (2 * a0))
    a_photo = a0 / samples['parallax']
    projected = np.hypot(C, H)
    a_rv = projected / np.sin(inclination)
    ratio = a_photo / a_rv
    period, ecc, systemic = (samples[k] for k in ('period', 'eccentricity', 'center_of_mass_velocity'))
    k_rv = (2 * np.pi * AU_KM / 86400) * projected / (period * np.sqrt(1 - ecc**2))
    w_rv = np.arctan2(C, H)
    w_ast = (np.arctan2(B - F, A + G) - np.arctan2(B + F, A - G)) / 2
    delta = (np.rad2deg(w_ast - w_rv) + 90) % 180 - 90
    wrapped_ast = (np.rad2deg(w_ast) + 90) % 180 - 90
    values = {
        'inclination_deg': np.rad2deg(inclination), 'a0_mas': a0, 'a_photo_AU': a_photo,
        'a_rv_AU': a_rv, 'a_rv_sin_i_AU': projected, 'a_photo_over_a_rv': ratio,
        'omega_ast_deg_mod_180': wrapped_ast, 'omega_rv_deg': np.rad2deg(w_rv),
        'omega_ast_minus_rv_deg_mod_180': delta, 'K_RV_km_s': k_rv,
        'rv_max_km_s': systemic + k_rv * (1 + ecc * np.cos(w_rv)),
        'rv_min_km_s': systemic + k_rv * (-1 + ecc * np.cos(w_rv)),
    }
    def transform(v):
        return {'catalogue_value': float(v[0]), **interval(v[1:])}

    eb = read_one('data/eclipse-model.csv')
    assert eb['model_type'] == 'TWOGAUSSIANS_WITH_ELLIPSOIDAL_ON_ECLIPSE1'
    phase = np.arange(20000) / 20000
    mag = (float(eb['geom_model_cosine_half_period_amplitude']) *
           np.cos(4 * np.pi * (phase - float(eb['geom_model_cosine_half_period_phase']))))
    for i in (1, 2):
        prefix = f'geom_model_gaussian{i}_'
        distance = (phase - float(eb[prefix + 'phase']) + 0.5) % 1 - 0.5
        mag += float(eb[prefix + 'depth']) * np.exp(-0.5 * (distance / float(eb[prefix + 'sigma']))**2)
    flux_at_phase = 10**(-0.4 * (mag - mag.min()))
    mean_inverse_flux = float(np.mean(1 / flux_at_phase))
    # For a constant luminous tertiary, beta(phase) = 1 - F_tertiary / F_total(phase).
    # Uniform sampling of inner phase is an illustrative limit for a mean photocentre.
    def eb_fraction_at_max(mean_fraction):
        return 1 - (1 - mean_fraction) / mean_inverse_flux
    conditional = {}
    for scale in scales:
        a_rel = a_rv * (1 + scale)
        mass = a_rel**3 / (period / DAYS_PER_YEAR)**2
        conditional[str(scale)] = {label: transform(val) for label, val in {
            'outer_relative_AU': a_rel, 'effective_total_mass_solar': mass,
            'RV_tracer_mass_solar': mass * scale / (1 + scale),
            'eclipsing_pair_mass_solar': mass / (1 + scale),
            'EB_light_fraction_same_side': (1 - ratio) / (1 + scale),
            'EB_light_fraction_opposite_side': (1 + ratio) / (1 + scale),
            'same_side_EB_light_fraction_at_photometric_max_uniform_phase': eb_fraction_at_max((1 - ratio) / (1 + scale)),
            'opposite_side_EB_light_fraction_at_photometric_max_uniform_phase': eb_fraction_at_max((1 + ratio) / (1 + scale)),
            'EB_barycentre_K_km_s': scale * k_rv,
        }.items()}
    losses = []
    for num, kind in enumerate(('primary', 'secondary'), 1):
        full_mag = float(eb[f'derived_{kind}_ecl_depth'])
        gaussian_mag = float(eb[f'geom_model_gaussian{num}_depth'])
        # The second derived eclipse center deviates slightly from Gaussian center.
        delta_phase = float(eb[f'derived_{kind}_ecl_phase']) - float(eb[f'geom_model_gaussian{num}_phase'])
        width = float(eb[f'geom_model_gaussian{num}_sigma'])
        gaussian_mag *= np.exp(-0.5 * (delta_phase / width)**2)
        min_flux = 10**(-0.4 * full_mag)
        losses.append(dict(full_fraction_of_maximum=1 - min_flux,
                           gaussian_fraction_of_maximum=min_flux * np.expm1(0.4 * np.log(10) * gaussian_mag)))
    massrow = read_one('data/binary-masses.csv')
    fractions = {key: float(massrow[key]) / (1 + float(massrow[key]))
                 for key in ('fluxratio', 'fluxratio_lower', 'fluxratio_upper')}
    result = {
        'source_id': row['source_id'], 'input_flag_bits': [i for i in range(64) if int(row['flags']) & (1 << i)],
        'normal_draws': DRAWS, 'seed': SEED, 'catalogue_efficiency': float(row['efficiency']),
        'correlation_min_eigenvalue': float(np.linalg.eigvalsh(correlation).min()),
        'correlation_det_geometric_root': float(np.linalg.det(correlation)**(1/n)),
        'orbital_quantities': {key: transform(v) for key,v in values.items()},
        'fixed_abs_clock_scales': conditional, 'eclipse_losses': losses,
        'uniform_inner_phase_mean_inverse_normalized_flux': mean_inverse_flux,
        'uniform_inner_phase_mean_normalized_flux': float(np.mean(flux_at_phase)),
        'minimum_phase_mean_EB_fraction_nonnegative_at_deepest_eclipse': 1 - float(flux_at_phase.min()) * mean_inverse_flux,
        'sum_two_gaussian_losses': sum(v['gaussian_fraction_of_maximum'] for v in losses),
        'Gaia_conditional_secondary_fraction': fractions,
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scales', nargs='+', type=float, default=[0.80, 0.965])
    args = parser.parse_args()
    if any(s <= 0 for s in args.scales):
        parser.error('Absolute clock scales must be positive.')
    print(json.dumps(analyze(args.scales), indent=2, allow_nan=False))

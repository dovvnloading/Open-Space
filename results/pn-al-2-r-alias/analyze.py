"""Reproduce the PN Al 2-R cadence experiment from the frozen OGLE and Gaia data."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
LONG = 501.418
SHORT = 0.9952767
SIDEREAL_FREQUENCY = 1.00273790935  # cycles per mean solar day
SEED = 91473


def harmonic(times, period):
    phase = 2 * np.pi * (times - times[0]) / period
    return np.column_stack((np.cos(phase), np.sin(phase)))


def sinusoid(times, values, periods):
    columns = [np.ones((len(times), 1))] + [harmonic(times, p) for p in periods]
    design = np.column_stack(columns)
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    prediction = design @ coefficients
    return coefficients, prediction


def rms(values):
    return float(np.sqrt(np.mean(values**2)))


def group_stats(groups, design, values):
    matrices, vectors, counts, labels = [], [], [], []
    for group in np.unique(groups):
        mask = groups == group
        if mask.sum() < 2:
            continue
        x = design[mask] - design[mask].mean(axis=0)
        y = values[mask] - values[mask].mean()
        matrices.append(x.T @ x)
        vectors.append(x.T @ y)
        counts.append(int(mask.sum()))
        labels.append(int(group))
    return np.array(matrices), np.array(vectors), counts, labels


def block_fit(times, values, period=SHORT):
    # Each observing visit is a sequence of transits separated by <= 12 hours.
    # A separate offset absorbs brightness differences between visits.
    groups = np.cumsum(np.r_[True, np.diff(times) > 0.5])
    design = harmonic(times, period)
    matrices, vectors, counts, labels = group_stats(groups, design, values)
    total = matrices.sum(axis=0)
    rhs = vectors.sum(axis=0)
    beta = np.linalg.solve(total, rhs)
    norm = float(np.linalg.norm(beta))
    draws = np.random.default_rng(SEED).integers(len(counts), size=(5000, len(counts)))
    sampled_matrices = matrices[draws].sum(axis=1)
    sampled_vectors = vectors[draws].sum(axis=1)
    valid = np.linalg.cond(sampled_matrices) < 1e8
    estimates = np.linalg.solve(sampled_matrices[valid], sampled_vectors[valid, :, None])[:, :, 0]
    radius = float(np.quantile(np.linalg.norm(estimates - beta, axis=1), 0.95))
    leverages = np.einsum('ij,gji->g', np.linalg.inv(total), matrices) / 2
    drop = np.linalg.solve(total - matrices, (rhs - vectors)[:, :, None])[:, :, 0]
    group_spans = [float(np.ptp(times[groups == k])) for k in labels]
    # Min/max phase response: inject a true 0.10 mag semi-amplitude signal before projection.
    injections = []
    for phase in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        injected = values + design @ (0.1 * np.array([np.cos(phase), np.sin(phase)]))
        _, alt, _, _ = group_stats(groups, design, injected)
        injections.append(np.linalg.norm(np.linalg.solve(total, alt.sum(axis=0))))
    return {
        'observations': len(times), 'repeat_visits': len(counts), 'repeat_observations': sum(counts),
        'max_visit_span_days': max(group_spans), 'largest_visit_n': max(counts),
        'largest_single_visit_information_fraction': float(max(leverages)),
        'two_coefficient_condition_number': float(np.linalg.cond(total)),
        'semiamplitude_mag': norm, 'omit_one_visit_range_mag': [float(min(np.linalg.norm(drop,axis=1))), float(max(np.linalg.norm(drop,axis=1)))],
        'drop_most_populated_visit_mag': float(np.linalg.norm(drop[np.argmax(counts)])),
        'bootstrap_draws_used': int(valid.sum()),
        'cluster_bootstrap_95pct_radius_plus_amplitude_mag': norm + radius,
        'injected_0p10_mag_recovery_range': [float(min(injections)), float(max(injections))],
    }


def analyze():
    expected = json.loads((ROOT / 'sources.json').read_text())
    for name, meta in expected.items():
        if hashlib.sha256((ROOT / 'inputs' / name).read_bytes()).hexdigest() != meta['sha256']:
            raise ValueError(f'Input hash differs for {name}')
    t, y, _ = np.loadtxt(ROOT / 'inputs/ogle-I.dat').T
    with (ROOT / 'inputs/gaia.ecsv').open() as stream:
        gaia = list(csv.DictReader(line for line in stream if not line.startswith('#')))
    long_b, long_pred = sinusoid(t, y, [LONG])
    short_b, short_pred = sinusoid(t, y, [SHORT])
    synthetic_b, _ = sinusoid(t, long_pred, [SHORT])
    joint_b, joint_pred = sinusoid(t, y, [LONG, SHORT])
    b_data, b_fake = short_b[1:], synthetic_b[1:]
    phase = np.arctan2(b_data[1], b_data[0]) - np.arctan2(b_fake[1], b_fake[0])
    difference = 1 / SHORT - 1 / LONG
    standalone = {
        'long_semiamplitude_mag': float(np.linalg.norm(long_b[1:])),
        'short_semiamplitude_mag': float(np.linalg.norm(b_data)),
        'long_only_replayed_short_semiamplitude_mag': float(np.linalg.norm(b_fake)),
        'data_minus_replay_phase_abs_cycles': abs(float(np.angle(np.exp(1j * phase)) / (2 * np.pi))),
        'joint_short_semiamplitude_mag': float(np.linalg.norm(joint_b[-2:])),
        'rms_constant_long_short_joint_mag': [rms(y-y.mean()),rms(y-long_pred),rms(y-short_pred),rms(y-joint_pred)],
        'frequency_separation_per_day': difference, 'difference_from_sidereal_per_day': difference - SIDEREAL_FREQUENCY,
        'ogle_span_days': float(np.ptp(t)), 'rayleigh_spacing_per_day': 1 / float(np.ptp(t)),
        'window_modulus_at_frequency_separation': float(abs(np.exp(2j*np.pi*(t-t[0])*difference).mean())),
    }
    blocks = {'OGLE_I': block_fit(t, y)}
    for band, timecol, magcol in [('G','g_transit_time','g_transit_mag'), ('RP','rp_obs_time','rp_mag')]:
        usable = [r for r in gaia if r['rejected_by_photometry'].lower() == 'false'
                  and r[f'variability_flag_{band.lower()}_reject'].lower() == 'false' and r[magcol]
                  and np.isfinite(float(r[magcol]))]
        series = sorted((float(r[timecol]), float(r[magcol])) for r in usable)
        gt, gy = np.transpose(series)
        blocks[f'Gaia_{band}'] = block_fit(gt, gy)
    # Check potential rounding discrepancies from the prose and sensitivity to the tabulated frequency.
    grid = np.linspace(0.990, 1.005, 301)
    groups = np.cumsum(np.r_[True, np.diff(t) > 0.5])
    estimates = []
    for period in grid:
        mat, vec, _, _ = group_stats(groups, harmonic(t, period), y)
        estimates.append(float(np.linalg.norm(np.linalg.solve(mat.sum(axis=0),vec.sum(axis=0)))))
    return {'published_periods_days': [LONG, SHORT], 'ogle_full_series': standalone,
            'within_visit': blocks, 'max_ogle_within_visit_amplitude_grid_0p990_to_1p005_step_0p00005_mag': max(estimates)}


if __name__ == '__main__':
    print(json.dumps(analyze(), indent=2))

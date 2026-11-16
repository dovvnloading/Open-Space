"""Exposure membership of the frozen December training and January comparison fits."""
from collections import Counter

FROZEN_EXPOSURES = (
    '2025W52_1A_0135_2', '2025W52_1A_0135_3', '2025W52_1A_0135_4',
    '2025W52_2A_0077_1', '2025W52_2A_0077_2', '2025W52_2A_0077_3',
    '2026W01_2A_0193_1', '2026W01_2A_0193_2', '2026W01_2A_0193_3',
    '2026W01_2A_0193_4',
)


def verify_fit_sample(rows):
    counts = Counter(row['observationid'] for row in rows)
    expected = set(FROZEN_EXPOSURES)
    missing = sorted(expected - counts.keys())
    extra = sorted(counts.keys() - expected)
    duplicates = sorted(key for key, count in counts.items() if count > 1)
    if missing or extra or duplicates:
        raise ValueError(
            f'Frozen SWAN fit exposure mismatch: missing={missing}, '
            f'unexpected={extra}, repeated={duplicates}'
        )


def residual_dof(fits, shared_parameters=0):
    """Residual dimensions left after fitting the epoch models and shared shape."""
    return max(0, sum(fit['n'] - fit['rank'] for fit in fits) - shared_parameters)

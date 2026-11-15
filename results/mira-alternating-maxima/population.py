"""Full-survey upper-tail rank test with leave-one-survey calibration across thresholds.

NumPy is required. Catalogue values, cycle counts, ties and star-specific
amplitudes are retained. With blocks, pick a random block-grid origin in
each simulated star, randomly reorder chunks of its chronologically ordered
observations, and place them at the original observed cycle numbers.
"""

import argparse
import numpy as np

from replay import load_series, longest, segments

BENCHMARKS = {"chi Cyg", "R Cyg", "RU Her", "R Hor"}
THRESHOLDS = np.arange(15, 34)


def draws(n, width, trials, rng):
    # Randomize the grid without wrapping an end-of-series measurement to
    # the first peak as if they were naturally adjacent.
    if width == 1:
        return np.argsort(rng.random((trials, n)), axis=1)
    orders = np.empty((trials, n), dtype=np.int16)
    offsets = rng.integers(0, min(width, n), trials)
    for offset in np.unique(offsets):
        chosen = np.flatnonzero(offsets == offset)
        start = np.arange(n)
        chunks = np.split(start, np.arange(offset or width, n, width))
        perm = rng.random((len(chosen), len(chunks))).argsort(axis=1)
        for j, p in zip(chosen, perm):
            orders[j] = np.concatenate([chunks[k] for k in p])
    return orders


def run_lengths(curves, cycles):
    trials = len(curves)
    step = np.sign(np.diff(curves, axis=1))
    run = np.ones(trials, dtype=np.int16)
    high = run.copy()
    prior = np.zeros(trials, dtype=np.int16)
    for i, jump in enumerate(np.diff(cycles)):
        if jump == 1:
            cur = step[:, i]
            run = np.where(cur == 0, 1, np.where(cur * prior == -1, run + 1, 2))
            prior = cur
        else:
            prior = np.zeros(trials, dtype=np.int16)
            run.fill(1)
        np.maximum(run, high, out=high)
    return high


def ar1_rank_curves(values, cycles, trials, rng):
    from scipy.special import ndtri
    from scipy.stats import rankdata

    n = len(values)
    scores = ndtri((rankdata(values) - 0.5) / n)
    contiguous = np.diff(cycles) == 1
    rho = 0.0
    if contiguous.sum() >= 5:
        before, after = scores[:-1][contiguous], scores[1:][contiguous]
        if before.std() > 0 and after.std() > 0:
            rho = np.clip(np.corrcoef(before, after)[0, 1], -0.98, 0.98)
    z = rng.standard_normal((trials, n))
    for j, gap in enumerate(np.diff(cycles), 1):
        decay = rho ** int(gap)
        z[:, j] = decay * z[:, j - 1] + np.sqrt(1 - decay * decay) * z[:, j]
    # Ordering the observed values by simulated ranks preserves every tie.
    result = np.empty((trials, n), dtype=values.dtype)
    np.put_along_axis(result, np.argsort(z, axis=1), np.sort(values)[None, :], axis=1)
    return result


def maxima(seq, trials, block, rng):
    cycles = np.array(sorted(seq))
    values = np.array([seq[e] for e in cycles], dtype=np.int16)
    if block == 0:
        return run_lengths(ar1_rank_curves(values, cycles, trials, rng), cycles)
    return run_lengths(values[draws(len(cycles), block, trials, rng)], cycles)


def calibrated(counts):
    # Row zero is observed. Rank every row as if it were the observed one.
    # This incorporates integer ties conservatively, including the row itself.
    n = len(counts)
    by_k = np.empty(counts.shape)
    for j, x in enumerate(counts.T):
        tails = np.bincount(x)[::-1].cumsum()[::-1]
        by_k[:, j] = tails[x] / n
    min_rank = by_k.min(axis=1)
    return (
        np.flatnonzero(by_k[0] == min_rank[0]),
        min_rank[0],
        int(np.count_nonzero(min_rank <= min_rank[0])),
        n,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument(
        "--blocks",
        nargs="+",
        type=int,
        default=(1, 4, 8, 12, 0),
        help="0 uses fitted per-star Gaussian-rank AR(1)",
    )
    opt = parser.parse_args()
    data = load_series()
    for block in opt.blocks:
        rng = np.random.default_rng(29112026)
        groups = {
            k: np.zeros((opt.trials + 1, len(THRESHOLDS)), dtype=np.int16)
            for k in ("all", "excluding benchmarks")
        }
        for name, seq in sorted(data.items()):
            high = max(longest(chain)[0] for chain in segments(seq))
            pseudo = maxima(seq, opt.trials, block, rng)
            counts = np.concatenate(([high], pseudo))[:, None] >= THRESHOLDS
            groups["all"] += counts
            if name not in BENCHMARKS:
                groups["excluding benchmarks"] += counts
        for cohort, counts in groups.items():
            where, p, ext, n = calibrated(counts)
            print(
                block,
                cohort,
                "observed counts",
                list(map(int, counts[0])),
                "strongest threshold(s)",
                list(map(int, THRESHOLDS[where])),
                "pointwise minimum",
                p,
                "calibrated rank",
                ext,
                "/",
                n,
                "p",
                round(ext / n, 5),
                flush=True,
            )


if __name__ == "__main__":
    main()

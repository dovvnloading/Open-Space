"""Joint forward/backward trapping fractions at fixed baselines for the two MPC fits."""

from concurrent.futures import ProcessPoolExecutor
import csv
import numpy as np
from current import TARGETS, setup
from dynamics import WORK, integrate, angle


def propagate(task):
    obj, sign = task
    n = 64
    sim = setup(obj, n, sign * 0.5)
    _, data = integrate(sim, sign * 10000000, 250)
    baselines = np.array(
        [
            50000,
            100000,
            200000,
            400000,
            1000000,
            2000000,
            3000000,
            5000000,
            7500000,
            10000000,
        ]
    )
    extrema = np.zeros((len(baselines), n + 1, 2))
    for k in range(n + 1):
        phi = np.rad2deg(np.unwrap(angle(data[:, k], *obj[2:])))
        phi -= phi[0]
        extrema[:, k, 0] = np.minimum.accumulate(phi)[baselines // 250]
        extrema[:, k, 1] = np.maximum.accumulate(phi)[baselines // 250]
    return baselines, extrema


if __name__ == "__main__":
    tasks = [(obj, s) for obj in TARGETS for s in (-1, 1)]
    with ProcessPoolExecutor(4) as pool:
        result = list(pool.map(propagate, tasks))
    with (WORK / "duration.csv").open("w") as f:
        wr = csv.writer(f)
        wr.writerow(
            [
                "object",
                "half_span_yr",
                "nominal_joint_range_deg",
                "both_directions_bounded",
                "past_only_test_bounded",
                "future_only_test_bounded",
                "draws",
            ]
        )
        for k, obj in enumerate(TARGETS):
            baselines, back = result[2 * k]
            _, fwd = result[2 * k + 1]
            for i, baseline in enumerate(baselines):
                b, future = back[i], fwd[i]
                span = np.maximum(b[:, 1], future[:, 1]) - np.minimum(
                    b[:, 0], future[:, 0]
                )
                wr.writerow(
                    [
                        obj[0],
                        baseline,
                        round(span[0], 4),
                        sum(span[1:] < 355),
                        sum(np.diff(b[1:], axis=1)[:, 0] < 355),
                        sum(np.diff(future[1:], axis=1)[:, 0] < 355),
                        64,
                    ]
                )

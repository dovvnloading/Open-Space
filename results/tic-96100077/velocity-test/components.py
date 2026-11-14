"""Test a labelled, three-line velocity series against the pre-existing outer orbit."""

import argparse
import json

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from orbit import CLOCK_SCALE, NOMINAL, predict_jd

NAMES = ["exposure_id", "jd_tcb", "rv_a", "err_a", "rv_b", "err_b", "rv_c", "err_c"]


def score(table):
    """A/B are matched inner-binary velocities, C is the putative outer tracer."""
    if not set(NAMES).issubset(table):
        raise ValueError(f"Required columns: {', '.join(NAMES)}")
    if table.exposure_id.isna().any() or table.exposure_id.duplicated().any():
        raise ValueError("Each exposure must have one unique identifier")
    if len(table) < 4:
        raise ValueError("At least four matched exposures are required")
    jd, a, ea, b, eb, c, ec = table[NAMES[1:]].to_numpy(dtype=float).T
    if (
        not np.isfinite([jd, a, ea, b, eb, c, ec]).all()
        or min(ea.min(), eb.min(), ec.min()) <= 0
    ):
        raise ValueError("Use finite velocities and strictly positive uncertainties")
    if np.any((jd < 2.4e6) | (jd > 2.6e6)):
        raise ValueError(
            "Use at least four matched exposures with full barycentric JD(TCB)"
        )
    x = predict_jd(jd) - NOMINAL[2]
    if np.ptp(x) < 1e-6:
        raise ValueError("A single predicted outer velocity cannot constrain its scale")

    def regression(y, sigma, scale):
        if scale is None:
            design = np.column_stack([np.ones(len(x)), x])
            coef, _, rank, _ = np.linalg.lstsq(
                design / sigma[:, None], y / sigma, rcond=None
            )
            if rank < 2:
                raise ValueError("Degenerate outer velocity sampling")
            intercept, scale = coef
        else:
            intercept = np.average(y - scale * x, weights=1 / sigma**2)
        residual = (y - intercept - scale * x) / sigma
        return dict(
            offset_kms=float(intercept),
            outer_scale=float(scale),
            chi2=float(residual @ residual),
        )

    def profile(scale):
        # w=q/(1+q) gives the secondary mass fraction and includes either single-line limit.
        def fit_at(w):
            # Profiling out a free inner orbital velocity per exposure leaves this residual.
            y = (1 - w) * a + w * b
            sigma = np.hypot((1 - w) * ea, w * eb)
            return regression(y, sigma, scale)

        grid = np.linspace(0, 1, 101)
        costs = np.array([fit_at(w)["chi2"] for w in grid])
        candidates = [0.0, 1.0]
        left = np.r_[np.inf, costs[:-1]]
        right = np.r_[costs[1:], np.inf]
        for i in np.flatnonzero((costs <= left) & (costs <= right)):
            optimum = minimize_scalar(
                lambda w: fit_at(w)["chi2"],
                method="bounded",
                bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]),
                options={"xatol": 1e-11},
            )
            candidates.append(float(optimum.x))
        w = min(candidates, key=lambda w: fit_at(w)["chi2"])
        return dict(
            **fit_at(w),
            weight_b=w,
            mass_ratio_b_over_a=w / (1 - w) if w < 1 else None,
            single_line_boundary=w in (0.0, 1.0),
        )

    inner_design = np.column_stack([np.ones(len(x)), x, b - a])
    return dict(
        free_pair_model_identifiable=bool(np.linalg.matrix_rank(inner_design) == 3),
        n_exposures=len(jd),
        outer_template_span_kms=float(np.ptp(x)),
        inner_velocity_difference_span_kms=float(np.ptp(a - b)),
        inner_pair={
            name: profile(s)
            for name, s in (
                ("constant", 0.0),
                ("predicted", CLOCK_SCALE),
                ("reversed", -CLOCK_SCALE),
                ("free", None),
            )
        },
        third_line={
            name: regression(c, ec, s)
            for name, s in (
                ("constant", 0.0),
                ("predicted", 1.0),
                ("reversed", -1.0),
                ("free", None),
            )
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("matched_csv")
    args = p.parse_args()
    print(
        json.dumps(
            score(pd.read_csv(args.matched_csv, dtype={"exposure_id": str})),
            indent=2,
            allow_nan=False,
        )
    )

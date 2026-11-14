"""Resolve independent timing coefficients for the two orbital halves and null-inject."""

import numpy as np
import pandas as pd

from ground_clock import (
    PERIOD_SECONDS,
    ROOT,
    SCALE,
    analyze,
    basis,
    gaia_delay,
    photometry,
    templates,
)


def separate_clocks(groups):
    rows = []
    for g in groups:
        phase, w = g["p"], g["w"]
        X = basis(phase, 4) * w[:, None]
        Q, _ = np.linalg.qr(X, mode="reduced")
        y = g["y"] * w
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        # The phase-0 and phase-0.5 halves are named by Gaia, irrespective of relative depth.
        grid = np.linspace(-0.12, 0.12, 241)
        center = grid[np.argmax(basis(grid, 4) @ b)]
        use = abs((phase - center + 0.5) % 1.0 - 0.5) < 0.25
        clock = g["rate"] * gaia_delay(g["t"]) * w
        drift = (g["t"] - g["t"].mean()) * g["rate"] * w
        Y = np.column_stack((y, clock * use, clock * ~use, drift))
        rows.append(Y - Q @ (Q.T @ Y))
    Y = np.vstack(rows)
    for survey in {g["survey"] for g in groups}:
        mask = np.concatenate(
            [np.full(len(g["t"]), g["survey"] == survey) for g in groups]
        )
        drift = Y[mask, -1]
        Y[mask, :3] -= np.outer(drift, drift @ Y[mask, :3] / (drift @ drift))
    y, X = Y[:, 0], Y[:, 1:3]
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ b
    _, cluster = np.unique(
        np.concatenate([g["block"] for g in groups]), return_inverse=True
    )
    scores = np.column_stack(
        [np.bincount(cluster, weights=col * residual) for col in X.T]
    )
    bread = np.linalg.inv(X.T @ X)
    cov = bread @ (scores.T @ scores) @ bread
    delta = np.array([1.0, -1.0])
    return [*b, *np.sqrt(np.diag(cov)), delta @ b, np.sqrt(delta @ cov @ delta)]


def inject(groups, scale):
    # Noise-free forward injection on the real irregular calendar with the measured
    # per-season templates. This checks sign and attenuation with method-fixed weights.
    injected = []
    for g in groups:
        X = basis(g["p"], 4)
        b = np.linalg.lstsq(X * g["w"][:, None], g["y"] * g["w"], rcond=None)[0]
        phase = g["p"] - scale * gaia_delay(g["t"]) / PERIOD_SECONDS
        injected.append(dict(g, y=basis(phase, 4) @ b))
    return analyze(injected, 4, controls=False)[0]["scale"]


def main():
    all_groups = templates(photometry()[0], 4)
    rows = []
    for survey in ("ASAS-3", "ASAS-SN", "combined"):
        groups = [
            g for g in all_groups if survey == "combined" or g["survey"] == survey
        ]
        rows.append(
            [
                survey,
                *separate_clocks(groups),
                inject(groups, 0.0),
                inject(groups, SCALE),
                inject(groups, -SCALE),
            ]
        )
    pd.DataFrame(
        rows,
        columns=[
            "survey",
            "phase0_scale",
            "phase0_5_scale",
            "phase0_se",
            "phase0_5_se",
            "difference",
            "difference_se",
            "recover_zero",
            "recover_negative",
            "recover_positive",
        ],
    ).to_csv(ROOT / "eclipse-controls.csv", index=False, float_format="%.9g")


if __name__ == "__main__":
    main()

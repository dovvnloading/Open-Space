"""Matched-template injection check for the illustrative +60 degree orientation."""

import numpy as np
import pandas as pd
from replay import HERE, GEOMETRIES, brightness, build_events, circle_area


def main():
    rng = np.random.default_rng(17018)
    geom = GEOMETRIES["ShortCadence"]
    theta = np.radians(60)
    (bc, rc, _), (bd, rd, _), _ = geom
    ntrials = 2000
    noise_cc = np.zeros(ntrials)
    measured = power = 0.0
    for _, (xc, xd, _), basis, resid, err in build_events(geom, "pdc", "best"):
        x = xd * np.cos(theta) - bd * np.sin(theta)
        y = xd * np.sin(theta) + bd * np.cos(theta)
        good = (np.hypot(xc, bc) < 1 - rc) & (np.hypot(x, y) < 1 - rd)
        overlap = circle_area(np.hypot(x - xc, y - bc), rc, rd)
        added = (1e6 * overlap * brightness(np.hypot((x + xc) / 2, (y + bc) / 2)) * good).mean(
            axis=1
        ) / err
        template = added - basis @ (basis.T @ added)
        power += template @ template
        measured += template @ resid
        null = np.array([template @ np.roll(resid, offset) for offset in range(len(resid))])
        noise_cc += rng.choice(null, ntrials)
    sig = np.sqrt(power)
    amplitude = measured / power
    rows = []
    for inject in (0.0, 0.5, 1.0):
        scores = (noise_cc + inject * power) / sig
        rows.append(
            [
                inject,
                ntrials,
                amplitude,
                sig,
                np.percentile(scores, 5),
                np.median(scores),
                np.percentile(scores, 95),
                np.mean(scores >= 4.0),
            ]
        )
    pd.DataFrame(
        rows,
        columns=[
            "injected_fraction",
            "trials",
            "observed_amplitude",
            "nominal_snr_full",
            "q05_score",
            "median_score",
            "q95_score",
            "fraction_score_ge_4",
        ],
    ).to_csv(HERE / "injections.csv", index=False, float_format="%.6f")


if __name__ == "__main__":
    main()

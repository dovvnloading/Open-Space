"""Stress-test the catalog ephemeris and the assumed residual correlation timescale."""

import numpy as np
import pandas as pd
from astropy.time import Time
from ground_clock import ROOT, SCALE, gaia_delay, photometry, profile, templates
from measure import ORBIT


def main():
    groups = templates(photometry()[0], 4)
    t = np.concatenate([g["t"] for g in groups])
    epoch = (
        Time(
            Time(2016.0, format="jyear", scale="tcb").jd + ORBIT.t_periastron,
            format="jd",
            scale="tcb",
        ).tdb.jd
        - 2457000
    )
    cuts = np.cumsum([len(g["t"]) for g in groups])[:-1]
    delays = np.split(gaia_delay(t), cuts)
    baseline = profile(groups, delays, 0.0, 4)
    rows = []
    # A sensitivity grid around the reference epoch and period, not a Gaia posterior.
    for dp in (-1, 0, 1):
        for de in (-1, 0, 1):
            p = ORBIT.period + dp * ORBIT.period_error
            t_alt = (
                epoch + (t - epoch - de * ORBIT.t_periastron_error) * ORBIT.period / p
            )
            alt = np.split(gaia_delay(t_alt), cuts)
            rows.append([dp, de, baseline - profile(groups, alt, SCALE, 4)])
    pd.DataFrame(
        rows, columns=["period_sigma", "epoch_sigma", "frozen_improvement"]
    ).to_csv(ROOT / "orbit-uncertainty-grid.csv", index=False, float_format="%.9g")
    rng = np.random.default_rng(96100077)
    res, signal = profile(groups, delays, 0.0, 4, sensitivity=True)
    scale = (res @ signal) / (signal @ signal)
    rows = []
    for label, block in [
        ("30_days", np.floor((t + 2457000) / 30)),
        ("90_days", np.floor((t + 2457000) / 90)),
        (
            "calendar_year",
            np.concatenate([np.repeat(g["year"], len(g["t"])) for g in groups]),
        ),
    ]:
        _, idx = np.unique(block, return_inverse=True)
        s = np.bincount(idx, weights=res * signal)
        sc = np.bincount(idx, weights=(res - scale * signal) * signal)
        se = np.sqrt(sc @ sc) / (signal @ signal)
        null = rng.choice([-1, 1], (4999, len(s))) @ s
        p = (1 + np.sum(null <= s.sum())) / (1 + len(null))
        rows.append([label, len(s), scale, se, p])
    pd.DataFrame(
        rows,
        columns=[
            "cluster",
            "groups",
            "linearized_scale",
            "cluster_se",
            "one_sided_null_p",
        ],
    ).to_csv(ROOT / "cluster-controls.csv", index=False, float_format="%.9g")


if __name__ == "__main__":
    main()

"""Compare the fixed signed template with calendar aliases at the observed times."""

import numpy as np
import pandas as pd
from astropy.time import Time
from ground_clock import ROOT, SCALE, gaia_delay, photometry, profile, templates
from measure import ORBIT


def main():
    d, _ = photometry()
    pooled = templates(d, 4)
    t0 = Time(2016.0, format="jyear", scale="tcb").tdb.jd - 2457000
    rows = []
    for survey in ("ASAS-3", "ASAS-SN", "combined"):
        groups = [g for g in pooled if survey == "combined" or g["survey"] == survey]
        t = np.concatenate([g["t"] for g in groups])
        cuts = np.cumsum([len(g["t"]) for g in groups])[:-1]
        zero = profile(groups, np.split(gaia_delay(t), cuts), 0.0, 4)
        for name, period in (
            ("Gaia", ORBIT.period),
            ("year", 365.25),
            ("half_year", 365.25 / 2),
            ("synodic_month", 29.53059),
        ):
            # Retain the Gaia waveform and amplitude, changing only angular rate.
            t_alt = t0 + (t - t0) * ORBIT.period / period
            score = []
            for shift in np.arange(100) / 100:
                delays = np.split(gaia_delay(t_alt, shift), cuts)
                score.append(zero - profile(groups, delays, SCALE, 4))
            best = np.argmax(score)
            rows.append([survey, name, period, score[0], best / 100, score[best]])
    pd.DataFrame(
        rows,
        columns=[
            "survey",
            "waveform_period",
            "days",
            "nominal_improvement",
            "best_phase",
            "best_improvement",
        ],
    ).to_csv(ROOT / "calendar-controls.csv", index=False, float_format="%.9g")


if __name__ == "__main__":
    main()

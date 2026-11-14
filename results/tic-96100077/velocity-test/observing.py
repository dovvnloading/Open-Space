"""Compute geometrical observing opportunities and the fixed Gaia/EB predictions."""

import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.utils import iers
from orbit import CLOCK_SCALE, DOSSIER, NOMINAL, predict_jd, sample_orbits, velocity

sys.path.insert(0, str(DOSSIER))
from measure import DEC, PERIOD, RA, REF, fit, gaia_delay

ROOT = Path(__file__).resolve().parent
TARGET = SkyCoord(RA * u.deg, DEC * u.deg)
SITE = EarthLocation.from_geodetic(
    -(70 + 44 / 60 + 4.543 / 3600) * u.deg,
    -(29 + 15 / 60 + 15.433 / 3600) * u.deg,
    2335 * u.m,
)
START = "2026-11-14T12:00:00"
NIGHTS = 365
SECONDS_PER_DAY = 86400.0


def timing_phase(bjd_tdb, timings):
    """TESS quadratic-plus-Gaia model for the common eclipse clock."""
    t = bjd_tdb - 2457000
    beta, _, _ = fit(timings, quadratic=True)
    x = t - timings.dropna(subset=["both"]).time.mean()
    delay = beta[0] + beta[1] * x + beta[2] * (x / 1000) ** 2 + beta[3] * gaia_delay(t)
    return ((t - REF - delay / SECONDS_PER_DAY) / PERIOD) % 1


def run():
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    t0 = Time(START, scale="utc", location=SITE)
    # Each row spans UTC noon to noon, keeping one Chilean night together.
    t = t0 + np.arange(NIGHTS * 144) * (10 * u.min)
    sky = AltAz(obstime=t, location=SITE, pressure=0 * u.hPa)
    alt = TARGET.transform_to(sky).alt.deg
    sun = get_sun(t).transform_to(sky).alt.deg
    visible = (alt >= 35) & (sun <= -18)
    bjd = t.tdb + t.light_travel_time(TARGET, kind="barycentric", ephemeris="builtin")
    timings = pd.read_csv(DOSSIER / "outputs/timings.csv")
    ph = timing_phase(bjd.jd, timings)
    which = np.where(ph < 0.5, 0.25, 0.75)
    quads = visible & (abs(ph - which) <= 0.06)
    elements = sample_orbits()
    gamma = elements[:, 2]
    rows = []
    for night in range(NIGHTS):
        ii = np.arange(night * 144, (night + 1) * 144)
        eligible = ii[visible[ii]]
        if not len(eligible):
            continue
        mid = eligible[len(eligible) // 2]
        velocities = velocity(bjd.tcb.jd[mid], elements)
        delta = velocities - gamma
        low, med, high = np.quantile(delta, [0.1, 0.5, 0.9])
        rows.append(
            dict(
                night=night,
                evening_utc_date=t[ii[0]].utc.isot[:10],
                representative_utc=t[mid].utc.isot[:16],
                visible_hours=len(eligible) / 6,
                quadrature_025=bool(np.any(quads[ii] & (which[ii] == 0.25))),
                quadrature_075=bool(np.any(quads[ii] & (which[ii] == 0.75))),
                delta_rv_p10=low,
                delta_rv_p50=med,
                delta_rv_p90=high,
                nominal_tracer_rv=float(predict_jd(bjd.tcb.jd[mid])),
                nominal_pair_gamma=NOMINAL[2]
                + CLOCK_SCALE * (predict_jd(bjd.tcb.jd[mid]) - NOMINAL[2]),
            )
        )
    nights = pd.DataFrame(rows).set_index("night")
    nights.to_csv(ROOT / "nightly-visibility.csv", float_format="%.4f")
    weeks = []
    for n in range(NIGHTS - 6):
        w = nights.reindex(range(n, n + 7))
        if (
            w.visible_hours.isna().any()
            or min(w.visible_hours) < 2
            or min(w.quadrature_025.sum(), w.quadrature_075.sum()) < 2
        ):
            continue
        weeks.append(
            dict(
                start=n,
                end=n + 6,
                low_score=w.delta_rv_p90.max(),
                high_score=w.delta_rv_p10.min(),
            )
        )
    if not weeks:
        raise ValueError("No qualifying week")
    weeks = pd.DataFrame(weeks)
    # Score the weakest night in a complete week at the conservative tail for each sign.
    high = weeks.loc[weeks.high_score.idxmax()]
    low = weeks.loc[weeks.low_score.idxmin()]
    chosen = {"low": low, "high": high}
    # One intermediate epoch after the high week tests the predicted fall back towards gamma.
    later = weeks[weeks.start >= high.end + 14]
    if len(later):
        chosen["return"] = later.loc[later.low_score.idxmin()]
    visits = []
    representative = {}
    for role, week in chosen.items():
        indexes = np.arange(int(week.start) * 144, (int(week.end) + 1) * 144)
        in_window = indexes[quads[indexes]]
        # Choose the visible grid sample closest to the quadrature of each passage.
        for group in np.split(in_window, np.flatnonzero(np.diff(in_window) > 1) + 1):
            if not len(group):
                continue
            j = group[np.argmin(abs(ph[group] - which[group]))]
            if j // 144 == int(week.start) + 3 and role not in representative:
                representative[role] = j
            v = velocity(bjd.tcb.jd[j], elements)
            eb = gamma + CLOCK_SCALE * (v - gamma)
            vq = np.quantile(v, [0.1, 0.5, 0.9])
            ebq = np.quantile(eb, [0.1, 0.5, 0.9])
            visits.append(
                dict(
                    role=role,
                    week_start=nights.at[int(week.start), "evening_utc_date"],
                    utc_midpoint=t[j].utc.isot[:16],
                    phase=ph[j],
                    nominal_quadrature=which[j],
                    altitude_deg=alt[j],
                    sun_altitude_deg=sun[j],
                    midpoint_window_start_utc=t[group[0]].utc.isot[:16],
                    midpoint_window_end_utc=t[group[-1]].utc.isot[:16],
                    tracer_p10=vq[0],
                    tracer_p50=vq[1],
                    tracer_p90=vq[2],
                    pair_gamma_p10=ebq[0],
                    pair_gamma_p50=ebq[1],
                    pair_gamma_p90=ebq[2],
                    nominal_tracer_rv=predict_jd(bjd.tcb.jd[j]),
                )
            )
    pd.DataFrame(visits).to_csv(
        ROOT / "quadrature-windows.csv", index=False, float_format="%.4f"
    )
    # Difference the same Gaia draw at both epochs to preserve the covariance of predictions.
    comparisons = []
    for a, b in (("low", "high"), ("high", "return")):
        if b not in representative:
            continue
        i, j = representative[a], representative[b]
        diff = velocity(bjd.tcb.jd[j], elements) - velocity(bjd.tcb.jd[i], elements)
        row = dict(
            first_utc=t[i].utc.isot[:16],
            second_utc=t[j].utc.isot[:16],
            first_role=a,
            second_role=b,
        )
        for label, values in (("tracer", diff), ("pair_gamma", CLOCK_SCALE * diff)):
            for label_percentile, val in zip(
                ("p025", "p50", "p975"), np.quantile(values, [0.025, 0.50, 0.975])
            ):
                row[f"{label}_{label_percentile}"] = val
        comparisons.append(row)
    pd.DataFrame(comparisons).to_csv(
        ROOT / "paired-predictions.csv", index=False, float_format="%.4f"
    )


if __name__ == "__main__":
    run()

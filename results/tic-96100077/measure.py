"""Measure an eclipsing-binary clock and compare it to the Gaia spectroscopic orbit."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "outputs"
MODEL = pd.read_csv(DATA / "gaia-model.csv").iloc[0]
ORBIT = pd.read_csv(DATA / "gaia-orbit.csv").iloc[0]
RA, DEC = ORBIT.ra, ORBIT.dec
# Gaia variability epochs use JD(TCB)-2455197.5. Convert both epoch and rate to TDB.
REF = (
    Time(
        2455197.5
        + MODEL.reference_time
        + MODEL.derived_primary_ecl_phase / MODEL.frequency,
        format="jd",
        scale="tcb",
    ).tdb.jd
    - 2457000
)
PERIOD = (1 - 1.550519768e-8) / MODEL.frequency


def clock(t, flux):
    """Common second-harmonic phase plus independent local fits to the two minima."""
    phase = ((t - REF) / PERIOD) % 1
    design = np.column_stack(
        [
            np.ones(len(t)),
            *(
                func(k * 2 * np.pi * phase)
                for k in range(1, 5)
                for func in (np.cos, np.sin)
            ),
        ]
    )
    coef = np.linalg.lstsq(design, flux, rcond=None)[0]
    common = np.arctan2(-coef[4], -coef[3]) / (4 * np.pi)
    offsets = [common]
    for center in [common, common + 0.5]:
        x = (phase - center + 0.5) % 1 - 0.5
        mask = abs(x) < 0.08
        if mask.sum() < 6 or np.ptp(x[mask]) < 0.10:
            offsets.append(np.nan)
            continue
        a, b, _ = np.polyfit(x[mask], flux[mask], 2)
        correction = -b / (2 * a) if a > 0 else np.nan
        offsets.append(common + correction if abs(correction) < 0.04 else np.nan)
    return np.array(offsets) * 86400 * PERIOD


def gaia_delay(t, shift=0):
    """Line-of-sight displacement of the Gaia RV primary, in light-seconds."""
    epoch = Time(2016.0, format="jyear", scale="tcb").jd + ORBIT.t_periastron
    tt = Time(np.asarray(t) + 2457000, format="jd", scale="tdb").tcb.jd
    m = 2 * np.pi * ((tt - epoch) / ORBIT.period + shift) % (2 * np.pi)
    e = ORBIT.eccentricity
    anomaly = m.copy()
    for _ in range(10):
        anomaly -= (anomaly - e * np.sin(anomaly) - m) / (1 - e * np.cos(anomaly))
    return 499.00478384 * (
        ORBIT.c_thiele_innes * (np.cos(anomaly) - e)
        + ORBIT.h_thiele_innes * np.sqrt(1 - e * e) * np.sin(anomaly)
    )


def windows(curves):
    result = []
    for sector, lc in curves.groupby("sector"):
        t = lc.time.to_numpy()
        for group in np.unique(np.floor((t - t.min()) / 3).astype(int)):
            idx = np.floor((t - t.min()) / 3).astype(int) == group
            tt = t[idx]
            phase_counts = np.bincount(
                (((tt - REF) / PERIOD) % 1 * 12).astype(int), minlength=12
            )
            row = {
                "sector": int(sector),
                "time": tt.mean(),
                "n": len(tt),
                "phase_bins": int((phase_counts > 0).sum()),
            }
            if np.ptp(tt) >= 1.5 and len(tt) >= 24 and np.all(phase_counts >= 1):
                row.update(
                    zip(
                        ("both", "primary", "secondary"),
                        clock(tt, lc.flux.to_numpy()[idx]),
                    )
                )
            result.append(row)
    return pd.DataFrame(result)


def fit(times, method="both", shift=0, use_orbit=True, local=False, quadratic=False):
    d = times.dropna(subset=[method])
    t = d.time.to_numpy()
    terms = [np.ones(len(t)), t - t.mean()]
    if quadratic:
        terms.append(((t - t.mean()) / 1000) ** 2)
    if local:
        terms.extend((d.sector == s).to_numpy() for s in sorted(d.sector.unique())[1:])
    if use_orbit:
        terms.append(gaia_delay(t, shift))
    x = np.column_stack(terms)
    y = d[method].to_numpy()
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    res = y - x @ beta
    return beta, float(np.sqrt(np.mean(res**2))), res


def evaluate(timings):
    results = {}
    for method in ["both", "primary", "secondary"]:
        for local in [False, True]:
            label = method + ("_within_sector" if local else "_global")
            b, rms, _ = fit(timings, method, local=local)
            _, linear, _ = fit(timings, method, local=local, use_orbit=False)
            results[label] = {
                "scale": b[-1],
                "rms_seconds": rms,
                "linear_rms_seconds": linear,
            }
    configurations = {
        "linear": {},
        "quadratic": {"quadratic": True},
        "within_sector": {"local": True},
    }
    sectors = sorted(timings.sector.unique())
    controls = []
    for name, settings in configurations.items():
        sub = {}
        b, rms, _ = fit(timings, **settings)
        _, base, _ = fit(timings, use_orbit=False, **settings)
        sub.update(scale=b[-1], rms_seconds=rms, baseline_rms_seconds=base)
        sub["leave_one_sector_out_scale"] = {
            str(s): fit(timings[timings.sector != s], **settings)[0][-1]
            for s in sectors
        }
        for i in range(100):
            b, score, _ = fit(timings, shift=i / 100, **settings)
            controls.append(
                {
                    "model": name,
                    "outer_phase_shift": i / 100,
                    "free_scale": b[-1],
                    "rms_seconds": score,
                }
            )
        sub["shifted_as_good_of_99"] = sum(
            row["rms_seconds"] <= rms
            for row in controls
            if row["model"] == name and row["outer_phase_shift"] > 0
        )
        results[name] = sub
    pd.DataFrame(controls).to_csv(
        OUT / "outer-phase-controls.csv", index=False, float_format="%.9f"
    )
    (OUT / "timing-fits.json").write_text(json.dumps(results, indent=2) + "\n")
    return results


if __name__ == "__main__":
    inputs = ROOT / "raw" / "curves.csv"
    if inputs.exists():
        timings = windows(pd.read_csv(inputs))
        timings.to_csv(OUT / "timings.csv", index=False, float_format="%.12f")
    evaluate(pd.read_csv(OUT / "timings.csv"))

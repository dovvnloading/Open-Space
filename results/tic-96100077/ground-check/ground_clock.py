"""Test the frozen Gaia/TESS timing prediction against independent ground photometry.

Place this folder next to the dossier's measure.py (or add its directory to PYTHONPATH).
The default run writes the tabulated outputs in this folder.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from measure import DEC, PERIOD, RA, REF, gaia_delay

SCALE = -0.8016409589619811
TARGET = SkyCoord(RA * u.deg, DEC * u.deg).cartesian.xyz.value
SECONDS_PER_DAY = 86400.0
PERIOD_SECONDS = PERIOD * SECONDS_PER_DAY


def projection(body, time):
    """Projected barycentric displacement in seconds using Astropy's builtin ephemeris."""
    return (TARGET @ get_body_barycentric(body, time).xyz.to_value(u.au)) * 499.00478384


def photometry():
    d = pd.read_csv(ROOT / "ground-photometry.csv.gz")
    # Validate photometric timestamp convention on overlapping old V observations.
    bg = d[d.channel == "bg"].jd.to_numpy()
    legacy = pd.read_csv(ROOT / "legacy-clock.csv").hjd.to_numpy()
    nearest = np.argmin(abs(bg[:, None] - legacy), axis=0)
    dt = (legacy - bg[nearest]) * SECONDS_PER_DAY
    matched = abs(dt) < 1200
    matched_t = Time(bg[nearest[matched]], format="jd", scale="utc")
    correction = projection("earth", matched_t) - projection("sun", matched_t)
    audit = {
        "legacy_total": len(legacy),
        "matched": int(matched.sum()),
        "median_abs_uncorrected_sec": np.median(abs(dt[matched])),
        "median_abs_corrected_sec": np.median(abs(dt[matched] - correction)),
        "correlation": np.corrcoef(dt[matched], correction)[0, 1],
    }
    sn = d.survey == "ASAS-SN"
    quality = (sn & (d.flag == "G")) | (~sn & d.flag.isin(["A", "B"]))
    d = d[
        quality & d.err.between(0, 0.2, inclusive="neither") & d.mag.between(5, 20)
    ].copy()
    for name, part in d.groupby("survey"):
        t = Time(part.jd.to_numpy(), format="jd", scale="utc")
        body = "sun" if name == "ASAS-3" else "earth"
        d.loc[part.index, "t"] = (
            t.tdb.jd - 2457000 + projection(body, t) / SECONDS_PER_DAY
        )
        d.loc[part.index, "year"] = (
            t.datetime64.astype("datetime64[Y]").astype(int) + 1970
        )
    return d, audit


def basis(phase, order):
    a = (2 * np.pi) * phase[:, None] * np.arange(1, order + 1)
    return np.column_stack((np.ones(len(phase)), np.cos(a), np.sin(a)))


def derivative(phase, coeff):
    order = len(coeff) // 2
    k = np.arange(1, order + 1)
    a = (2 * np.pi) * phase[:, None] * k
    return (
        np.sin(a) @ (k * coeff[1 : order + 1]) - np.cos(a) @ (k * coeff[order + 1 :])
    ) * (2 * np.pi / PERIOD_SECONDS)


def templates(d, order):
    result = []
    for (survey, channel, year), g in d.groupby(["survey", "channel", "year"]):
        t = g.t.to_numpy()
        phase = ((t - REF) / PERIOD) % 1
        if len(t) < 30 or len(np.unique((12 * phase).astype(int))) < 9:
            continue
        X = basis(phase, order)
        y = g.mag.to_numpy()
        weights = np.ones(len(y))
        for _ in range(5):
            b = np.linalg.lstsq(X * weights[:, None], y * weights, rcond=None)[0]
            r = y - X @ b
            sigma = max(0.005, 1.4826 * np.median(abs(r - np.median(r))))
            weights = (
                np.sqrt(np.minimum(1.0, 1.5 * sigma / np.maximum(abs(r), 1e-9))) / sigma
            )
        result.append(
            dict(
                survey=survey,
                channel=channel,
                year=int(year),
                t=t,
                y=y,
                p=phase,
                w=weights,
                rate=derivative(phase, b),
                block=np.floor((t + 2457000) / 30).astype(int),
            )
        )
    return result


def profile(groups, delays, scale, order, local_rate=False, sensitivity=False):
    """Profile camera-season Fourier shapes, plus a shared clock drift per survey."""
    rows = []
    for g, delay in zip(groups, delays):
        phase = g["p"] - scale * delay / PERIOD_SECONDS
        X = basis(phase, order) * g["w"][:, None]
        rate = (g["t"] - g["t"].mean()) * g["rate"] * g["w"]
        if local_rate:
            X = np.column_stack([X, rate])
        Q, _ = np.linalg.qr(X, mode="reduced")
        col = np.column_stack((g["y"] * g["w"], rate, g["rate"] * delay * g["w"]))
        rows.append(col - Q @ (Q.T @ col))
    stacked = np.vstack(rows)
    y, drift, signal = stacked.T
    for survey in {g["survey"] for g in groups}:
        use = np.concatenate(
            [np.full(len(g["t"]), g["survey"] == survey) for g in groups]
        )
        if not local_rate:
            v = drift[use]
            for arr in (y, signal):
                arr[use] -= v * (v @ arr[use] / (v @ v))
    if sensitivity:
        return y, signal
    return y @ y


def analyze(groups, order, local_rate=False, controls=True):
    t = np.concatenate([g["t"] for g in groups])
    cuts = np.cumsum([len(g["t"]) for g in groups])[:-1]
    orbit = lambda shift: np.split(gaia_delay(t, shift), cuts)
    delays = orbit(0.0)
    f = lambda scale: profile(groups, delays, scale, order, local_rate)
    opt = minimize_scalar(f, bounds=(-5, 5), method="bounded")
    zero, fixed = f(0.0), f(SCALE)
    # Efficient score at the zero-delay null, retaining all original observing dates.
    res, signal = profile(groups, delays, 0.0, order, local_rate, sensitivity=True)
    blocks = np.concatenate([g["block"] for g in groups])
    _, index = np.unique(blocks, return_inverse=True)
    score = np.bincount(index, weights=res * signal)
    lin_scale = (res @ signal) / (signal @ signal)
    # Independent 30-day sign flips under the null; a conditional sensitivity check.
    rng = np.random.default_rng(96100077)
    null = rng.choice([-1, 1], (4999, len(score))) @ score
    p = (1 + np.sum(null <= score.sum())) / (1 + len(null))
    cluster = np.bincount(index, weights=(res - lin_scale * signal) * signal)
    se = np.sqrt(cluster @ cluster) / (signal @ signal)
    out = dict(
        points=len(t),
        groups=len(groups),
        blocks=len(score),
        scale=opt.x,
        linearized_scale=lin_scale,
        cluster_se=se,
        wild_one_sided_p=p,
        frozen_improvement=zero - fixed,
        free_improvement=zero - opt.fun,
    )
    control_table = []
    if controls:
        for i in range(100):
            cost = profile(groups, orbit(i / 100), SCALE, order, local_rate)
            control_table.append([i / 100, zero - cost])
        out["wrong_phases_as_good"] = sum(
            v >= zero - fixed for _, v in control_table[1:]
        )
    return out, control_table


def main():
    d, audit = photometry()
    import json

    print(json.dumps(audit, indent=2))
    (ROOT / "timestamp-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    summaries, phase_rows, station_rows = [], [], []
    for order in (4, 3, 5):
        all_groups = templates(d, order)
        for survey in ("ASAS-3", "ASAS-SN", "combined"):
            groups = [
                g for g in all_groups if g["survey"] == survey or survey == "combined"
            ]
            for local in (False, True) if order == 4 else (False,):
                summary, controls = analyze(
                    groups, order, local, controls=order == 4 and not local
                )
                row = dict(survey=survey, harmonics=order, local_rate=local, **summary)
                summaries.append(row)
                print(row, flush=True)
                phase_rows.extend([survey, *c] for c in controls)
            if order == 4:
                for year in sorted({g["year"] for g in groups}):
                    other = [g for g in groups if g["year"] != year]
                    row, _ = analyze(other, order, controls=False)
                    station_rows.append(
                        [survey, year, row["scale"], row["frozen_improvement"]]
                    )
    pd.DataFrame(summaries).to_csv(
        ROOT / "ground-fits.csv", index=False, float_format="%.8g"
    )
    pd.DataFrame(
        phase_rows, columns=["survey", "phase_shift", "frozen_improvement"]
    ).to_csv(ROOT / "ground-phase-controls.csv", index=False, float_format="%.8g")
    pd.DataFrame(
        station_rows, columns=["survey", "omitted_year", "scale", "frozen_improvement"]
    ).to_csv(ROOT / "ground-year-controls.csv", index=False, float_format="%.8g")


if __name__ == "__main__":
    main()

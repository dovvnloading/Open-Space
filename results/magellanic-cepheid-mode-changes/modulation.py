"""Measure the slower envelope of both radial overtones using short time windows."""

import csv
import json
import numpy as np
from scipy.optimize import minimize_scalar
from screen import catalog, read_curves, estimate, ROOT, INPUT
from follow_fetch import TARGETS


def eros(target):
    rows = [
        s.split("\t")
        for s in (INPUT / f"{target}_eros_xcat.tsv").read_text().splitlines()
        if s and s[0].isdigit() and "\tLC\t" in s
    ]
    if not rows:
        return
    row = min(rows, key=lambda r: float(r[0]))
    source = row[2].strip()
    if float(row[0]) > 1:
        return
    for k, band in [(2, "R"), (4, "B")]:
        points = []
        for l in (INPUT / f"{target}_eros_xtime.tsv").read_text().splitlines():
            r = l.split("\t")
            if (
                len(r) > k + 1
                and r[0].strip() == source
                and r[k].strip()
                and r[k + 1].strip()
            ):
                t, y, e = map(float, [r[1], r[k], r[k + 1]])
                if 0 < e < 0.2 and np.isfinite(y):
                    points.append([t, y, e])
        if points:
            yield band, source, float(row[0]), np.array(points)


def segments(data, p, q, width, shift):
    t = data[:, 0]
    keys = np.floor((t - shift) / width)
    result = []
    for key in np.unique(keys):
        chunk = data[keys == key]
        if len(chunk) < 24 or np.ptp(chunk[:, 0]) < 28:
            continue
        amplitudes = estimate(chunk, p, q)
        if amplitudes is not None:
            result.append((np.mean(chunk[:, 0]), *amplitudes, len(chunk)))
    return np.array(result).reshape(-1, 4)


def design(t, period):
    x = 2 * np.pi * (t - 3500) / period
    return np.column_stack(
        [np.ones(len(x)), np.sin(x), np.cos(x), np.sin(2 * x), np.cos(2 * x)]
    )


def model(arr):
    t = arr[:, 0]
    y = arr[:, 2]

    def fit(P):
        X = design(t, P)
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        return np.mean((y - X @ b) ** 2), b

    grid = np.arange(600, 2602, 4)
    key = min(grid, key=lambda p: fit(p)[0])
    P = minimize_scalar(
        lambda p: fit(p)[0], bounds=(key - 4, key + 4), method="bounded"
    ).x
    err, b = fit(P)
    return P, b


def score(seq, period, beta):
    if seq.shape[0] < 4:
        return {}
    pred = design(seq[:, 0], period) @ beta
    obs = seq[:, 2]
    return dict(
        n=len(obs),
        r=float(np.corrcoef(obs, pred)[0, 1]),
        mean=float(np.mean(obs)),
        rms=float(np.sqrt(np.mean((pred - obs) ** 2))),
        constant_rms=float(np.sqrt(np.mean((obs - beta[0]) ** 2))),
    )


def main():
    names = ["OGLE-SMC-CEP-" + n for n in TARGETS]
    curves = dict(read_curves(names))
    stars = catalog()
    allrows = []
    results = []
    for name in names:
        k = name[-4:]
        pair, p, q = stars[name]
        er = list(eros(k))
        for width, shift in [(190, 0), (190, 95), (160, 0), (220, 0), (270, 0)]:
            III = segments(curves[name, 3, "I"], p, q, width, shift)
            IV = segments(curves[name, 4, "I"], p, q, width, shift)
            if len(III) < 8 or len(IV) < 4:
                continue
            P, bet = model(III)
            # EROS supplies independent measurements in two bands. Compare their envelopes at the OGLE phase.
            both = np.vstack([III, IV])
            allP, allbet = model(both)
            result = {
                "star": name,
                "width": width,
                "shift": shift,
                "old_period": P,
                "old_n": len(III),
                "ogle_iv": score(IV, P, bet),
                "all_period": allP,
                "all_n": len(both),
                "r_A1_A2": float(np.corrcoef(both[:, 1], both[:, 2])[0, 1]),
                "eros": {},
            }
            if width == 190 and shift == 0:
                for phase, ar in [("OGLE II/III I", III), ("OGLE IV I", IV)]:
                    allrows.extend((name, phase, "", *r) for r in ar)
            for band, id_, sep, data in er:
                seq = segments(data, p, q, width, shift)
                result["eros"][band] = dict(
                    id=id_, sep_arcsec=sep, **score(seq, allP, allbet)
                )
                if width == 190 and shift == 0:
                    allrows.extend((name, "EROS " + band, id_, *r) for r in seq)
            results.append(result)
            if width == 190 and shift == 0:
                print(
                    k,
                    round(P),
                    result["ogle_iv"],
                    round(allP),
                    [(x, v.get("r")) for x, v in result["eros"].items()],
                    flush=True,
                )
    with (ROOT / "segments.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "star",
                "survey_passband",
                "eros_id",
                "mean_HJD_minus_2450000",
                "semiamp_1O_mag",
                "semiamp_2O_mag",
                "observations",
            ]
        )
        w.writerows(allrows)
    (ROOT / "modulation.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()

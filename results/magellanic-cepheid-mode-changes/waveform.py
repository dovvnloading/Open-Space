"""Fit a double-overtone modulation ephemeris in OGLE-II/III and score other instruments."""

import csv
import json
import numpy as np
from scipy.optimize import minimize
from screen import catalog, read_curves, ROOT, INPUT
from follow_fetch import TARGETS


def components(time, frequencies, mod_phase=0.0):
    f, g, h = frequencies
    t = time - 3500.0
    x, y, z = 2 * np.pi * f * t, 2 * np.pi * g * t, 2 * np.pi * h * t + mod_phase
    # Columns grouped by which overtone they belong to for external colour calibration.
    args = [x, 2 * x, 3 * x, y, 2 * y, x + y, x - y]
    groups = [0, 0, 0, 1, 1, 2, 2]
    if h:
        for k in [1, 2]:
            for phase, group in [(x, 0), (y, 1)]:
                args.extend([phase - k * z, phase + k * z])
                groups.extend([group, group])
    X = np.column_stack(
        [np.ones(len(t))] + [tr(a) for a in args for tr in [np.sin, np.cos]]
    )
    return X, np.array([-1] + [g for g in groups for _ in range(2)])


def solve(X, y, w):
    coeff = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]
    return coeff, np.sum(((y - X @ coeff) * w) ** 2)


def fit(data, p, q):
    time, y, e = data.T
    w = 1 / np.hypot(e, 0.009)
    use = np.full(len(y), True)
    f0, g0 = 1 / p, 1 / q

    def loss(pars, ret=False):
        X, _ = components(time, pars)
        coeff, rss = solve(X[use], y[use], w[use])
        return (coeff, rss) if ret else rss

    def optimize(start, mod):
        # Use frequency offsets in cycles per thousand days for numerical conditioning.
        def freq(v):
            return (
                f0 + v[0] * 0.001,
                g0 + v[1] * 0.001,
                (v[2] * 0.001 if mod else 0.0),
            )

        s = np.array(
            [(start[0] - f0) * 1000, (start[1] - g0) * 1000]
            + ([start[2] * 1000] if mod else [])
        )
        bounds = [(-0.6, 0.6), (-0.8, 0.8)] + (
            [(1000 / 3600, 1000 / 520)] if mod else []
        )
        out = minimize(
            lambda v: loss(freq(v)),
            s,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 2e-12, "maxiter": 100},
        )
        return freq(out.x)

    single = optimize((f0, g0, 0), False)
    X, _ = components(time, single)
    b, _ = loss(single, True)
    r = (y - X @ b) * w
    use = abs(r - np.median(r)) < max(8, 7 * 1.4826 * np.median(abs(r - np.median(r))))
    single = optimize(single, False)
    grid = sorted(
        [
            (loss((single[0], single[1], h)), h)
            for h in np.linspace(1 / 3500, 1 / 530, 110)
        ]
    )
    candidates = []
    for _, h in grid:
        if any(abs(h - k) < 0.00009 for k in candidates):
            continue
        candidates.append(h)
        if len(candidates) == 6:
            break
    options = [optimize((single[0], single[1], h), True) for h in candidates]
    best = min(options, key=loss)
    return [(single, *loss(single, True)), (best, *loss(best, True))], int(sum(use))


def prediction(t, model, phi=0):
    pars, b, _ = model
    X, groups = components(t, pars, phi)
    return np.column_stack([X[:, groups == g] @ b[groups == g] for g in range(3)])


def time_hjd(mjd, ra, dec):
    jd = mjd + 2400000.5
    days = jd - 2451545
    lon = np.deg2rad((280.460 + 0.9856474 * days) % 360)
    anomaly = np.deg2rad((357.528 + 0.9856003 * days) % 360)
    lon += np.deg2rad(1.915 * np.sin(anomaly) + 0.020 * np.sin(2 * anomaly))
    obl = np.deg2rad(23.439 - 0.0000004 * days)
    radius = 1.00014 - 0.01671 * np.cos(anomaly) - 0.00014 * np.cos(2 * anomaly)
    a, d = np.deg2rad([ra, dec])
    projection = np.cos(lon) * np.cos(a) * np.cos(d) + np.sin(lon) * (
        np.cos(obl) * np.sin(a) * np.cos(d) + np.sin(obl) * np.sin(d)
    )
    return jd - 2450000.0 - radius * projection * 499.005 / 86400


def delve(target, cut=None):
    rows = list(csv.DictReader((INPUT / f"{target}_delve_phot.csv").open()))
    best = {}
    for r in rows:
        e = float(r["err"])
        if (
            not 0 < e < 0.08
            or not 12 < float(r["mag"]) < 23
            or r["filter"] not in ("g", "r", "i", "z")
        ):
            continue
        key = r["filter"], r["exposure"].split("_")[0]
        if key not in best or e < float(best[key]["err"]):
            best[key] = r
    rows = list(best.values())
    ra = np.median([float(r["ra"]) for r in rows])
    dec = np.median([float(r["dec"]) for r in rows])
    t = time_hjd(np.array([float(r["mjd"]) for r in rows]), ra, dec)
    data = np.column_stack(
        [t, [float(r["mag"]) for r in rows], [float(r["err"]) for r in rows]]
    )
    bands = np.array([r["filter"] for r in rows])
    if cut is None:
        return data, bands
    return data[t > cut], bands[t > cut]


def scores(data, bands, models):
    t, y, e = data.T
    # Equal total leverage per night per filter. Same calibration degrees of freedom for each model.
    nights = np.floor(t).astype(int)
    valid = np.ones(len(t), bool)
    for band in np.unique(bands):
        if len(np.unique(nights[bands == band])) < 6:
            valid[bands == band] = False
    t, y, e, bands, nights = t[valid], y[valid], e[valid], bands[valid], nights[valid]
    if len(t) < 10:
        return None
    errors = []
    for model in models:
        pred = prediction(t, model)
        residual = np.zeros(len(t))
        for band in np.unique(bands):
            mask = bands == band
            nn = np.unique(nights[mask], return_inverse=True, return_counts=True)
            w = 1 / np.sqrt(nn[2][nn[1]])
            # One passband offset and two colour scales, integrating combination terms with first mode.
            X = np.column_stack(
                [np.ones(sum(mask)), pred[mask, 0] + pred[mask, 2], pred[mask, 1]]
            )
            coeff, _ = solve(X, y[mask], w)
            residual[mask] = y[mask] - X @ coeff
        per_night = [np.mean(residual[nights == n] ** 2) for n in np.unique(nights)]
        errors.append(per_night)
    a, b = np.array(errors)
    return dict(
        n=len(t),
        nights=len(a),
        range=t[[np.argmin(t), np.argmax(t)]].tolist(),
        rms_const=float(np.sqrt(a.mean())),
        rms_mod=float(np.sqrt(b.mean())),
        win=int(sum(b < a)),
        loss=int(sum(a < b)),
    )


if __name__ == "__main__":
    names = ["OGLE-SMC-CEP-" + n for n in TARGETS]
    curves = dict(read_curves(names))
    st = catalog()
    results = {}
    for name in names:
        k = name[-4:]
        pair, p, q = st[name]
        train = curves[name, 3, "I"]
        test = curves[name, 4, "I"]
        models, n = fit(train, p, q)
        out = dict(
            train_used=n,
            models=[
                {"freq": m[0], "beta": m[1].tolist(), "chi2": m[2]} for m in models
            ],
            ogle=scores(test, np.repeat("I", len(test)), models),
            decam=scores(*delve(k), models),
            decam_after_ogle=scores(*delve(k, cut=max(test[:, 0])), models),
        )
        results[k] = out
        print(
            k,
            round(1 / models[1][0][2], 1),
            [
                (
                    key,
                    (
                        None
                        if v is None
                        else (
                            v["nights"],
                            round(v["rms_const"], 4),
                            round(v["rms_mod"], 4),
                        )
                    ),
                )
                for key, v in out.items()
                if key in ["ogle", "decam", "decam_after_ogle"]
            ],
            flush=True,
        )
        (ROOT / "clock_results.json").write_text(json.dumps(results, indent=2) + "\n")

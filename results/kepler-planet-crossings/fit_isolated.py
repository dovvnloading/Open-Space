from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import least_squares, minimize_scalar
from batman import _quadratic_ld


def model(x, z, lds):
    b, r, dur = z
    v = 2 * np.sqrt((1 + r) ** 2 - b * b) / (dur / 24)
    # three-minute bins; planet contact profiles change on a timescale of minutes.
    t = x[:, None] + np.array([-0.0007, 0, 0.0007])
    sep = np.ascontiguousarray(np.hypot(t * v, b).ravel())
    return (_quadratic_ld._quadratic_ld(sep, r, *lds, 1).reshape(len(x), 3).mean(axis=1) - 1) * 1e6


def solve(d, lds):
    q = d.copy()
    q["offset"] = 0.0
    z = np.array([0.55, 0.05, 3.6])
    for loop in range(4):
        for ep, g in q.groupby("epoch"):
            x, y, e = g[["dt", "flux", "error"]].values.T

            def f(off):
                resid = (y - model(x - off, z, lds)) / e
                return np.sum(resid**2)

            of = minimize_scalar(f, bounds=(-0.015, 0.015), method="bounded").x
            q.loc[g.index, "offset"] = of
        x = q.dt.values - q.offset.values
        y = q.flux.values
        e = q.error.values
        fit = least_squares(
            lambda a: (model(x, a, lds) - y) / e,
            z,
            bounds=([0, 0.02, 2.7], [0.94, 0.12, 4.4]),
            diff_step=1e-4,
            x_scale=[0.3, 0.01, 0.1],
        )
        z = fit.x
    return z, q


if __name__ == "__main__":
    d = pd.read_csv("isolated.csv.gz")
    sol = []
    for L in [(0.5256, 0.186), (0.43, 0.23), (0.62, 0.12)]:
        for name, g in d.groupby("planet"):
            z, q = solve(g, L)
            print(L, name, np.round(z, 5), flush=True)
            sol.append([*L, name, *z, len(g), g.epoch.nunique()])
            if L == (0.5256, 0.186):
                e = np.sort(g.epoch.unique())
                for i, eps in enumerate(np.array_split(e, 4)):
                    zz, qq = solve(g[g.epoch.isin(eps)], L)
                    print("block", name, i, np.round(zz, 4), flush=True)
                    sol.append([*L, name + "_" + str(i), *zz, len(qq), len(eps)])
                q["minute_bin"] = np.round((q.dt - q.offset) * 1440 / 2).astype(int) * 2
                b = q.groupby("minute_bin").agg(
                    time=("dt", "size"), flux=("flux", "mean"), flux_std=("flux", "std")
                )
                b["error"] = b.flux_std / np.sqrt(b.time)
                b.to_csv("fold_" + name + ".csv", float_format="%.5f")
    pd.DataFrame(
        sol, columns="u1 u2 planet impact radius duration_hours bins events".split()
    ).to_csv("isolated_fits.csv", index=False)

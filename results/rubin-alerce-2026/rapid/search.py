import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.linalg import lstsq


def scan(g):
    g = g.sort_values("mjd")
    t = (g.mjd.values - g.mjd.min()) * 1440
    y = g.scienceflux.values
    e = np.hypot(g.sciencefluxerr.values, 0.005 * abs(np.median(y)))
    w = 1 / e**2
    if len(t) < 16 or t[-1] < 14 or np.max(np.diff(t)) > 0.6 * t[-1]:
        return
    # Single bad fluxes can supply otherwise convincing periodogram power.
    m = np.median(y)
    mad = 1.4826 * np.median(abs(y - m))
    good = abs(y - m) < 8 * np.maximum(mad, e)
    if good.sum() < 16:
        return
    t, y, w = t[good], y[good], w[good]
    if t[-1] <= 12:
        return
    u = t - np.average(t, weights=w)
    trend = np.stack([np.ones_like(t), u], axis=1)
    null = lstsq(trend * np.sqrt(w[:, None]), y * np.sqrt(w), rcond=None)[0]
    res = y - trend @ null
    f = np.arange(max(1 / 60, 2 / t[-1]), 1 / 4, 1 / (6 * t[-1]))
    if len(f) < 4:
        return
    basis = np.stack(
        [np.sin(2 * np.pi * f[:, None] * t), np.cos(2 * np.pi * f[:, None] * t)], axis=2
    )
    fit = np.linalg.solve(
        trend.T @ (trend * w[:, None]), np.einsum("ni,n,fnk->fik", trend, w, basis)
    )
    basis = basis - np.einsum("ni,fik->fnk", trend, fit)
    mat = np.einsum("fni,n,fnj->fij", basis, w, basis)
    rhs = np.einsum("fni,n->fi", basis, res * w)
    coef = np.linalg.solve(mat, rhs[..., None])[..., 0]
    delta = np.einsum("ij,ij->i", rhs, coef)
    j = delta.argmax()
    frac = delta[j] / np.dot(res * w, res)
    amp = np.linalg.norm(coef[j])
    return dict(
        n=len(t),
        df=1 / t[-1],
        minutes=t[-1],
        mjd=g.mjd.min(),
        p=1 / f[j],
        amp=amp,
        median=m,
        frac=frac,
        delta=delta[j],
        scatter=mad,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit single-band nightly sinusoids to Rubin forced science fluxes."
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    p = pd.concat(pd.read_csv(path) for path in args.inputs)
    p = p[np.isfinite(p.scienceflux) & (p.sciencefluxerr > 0)]
    p["night"] = np.floor(p.mjd - 0.5).astype(int)
    rows = []
    for key, g in p.groupby(["oid", "night", "band"]):
        r = scan(g)
        if r:
            rows.append(dict(zip(["oid", "night", "band"], key)) | r)
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    print(f"{len(p)} input measurements; {len(out)} eligible object/night/band series.")

"""Exploratory Kepler long-cadence double-transit residual screen."""

import pandas as pd, numpy as np, itertools, glob, sys
from pathlib import Path
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.optimize import least_squares

kois = pd.read_csv("koi.csv")
pairs = pd.read_csv("pairs.csv")


def readtimes():
    a = pd.read_csv("times.tsv", sep="\t", comment="#", skipinitialspace=True, low_memory=False)
    a = a[["KOI", "tn", "O-C", "Out", "Over"]].apply(pd.to_numeric, errors="coerce").dropna()
    a["tn"] += 67
    a["tt"] = a.tn + a["O-C"] / 1440
    return a


hs = readtimes()


def near(ts, q):
    i = np.searchsorted(ts, q).clip(1, len(ts) - 1)
    return np.where(abs(ts[i] - q) < abs(ts[i - 1] - q), ts[i], ts[i - 1])


def readstar(k):
    files = sorted(glob.glob(f"cache/kplr{k:09}-*_llc.npz"))
    data = []
    for f in files:
        d = np.load(f)
        t = d["TIME"]
        y = d["PDCSAP_FLUX"]
        valid = np.isfinite(t + y) & (d["SAP_QUALITY"] == 0)
        if np.sum(valid) < 20:
            continue
        data.append((t[valid], y[valid] / np.nanmedian(y[valid]) * 1e6))
    return [np.concatenate(x) for x in zip(*data)] if data else ([], [])


def prior(r):
    ep = r.koi_time0bk
    per = r.koi_period
    dur = r.koi_duration / 24
    # Archive KOI ephemerides sometimes drift significantly from Holczer's. Use their ephemeris when present.
    h = hs[(hs.KOI == float(r.kepoi_name[1:]))].sort_values("tn").drop_duplicates("tn")
    if len(h) > 3:
        N = np.round((h.tn - h.tn.iloc[0]) / per)
        coef = np.polyfit(N, h.tn, 1)
        if abs(coef[0] / per - 1) < 0.001:
            per, ep = coef
    ts = ep + np.arange(np.floor((110 - ep) / per), np.ceil((1610 - ep) / per) + 1) * per
    if len(h) > 3:
        clean = h[(h.Over == 0) & (h.Out == 0)]
        if len(clean) > 2:
            ts += np.interp(ts, clean.tn, clean["O-C"]) / 1440
        # do not allow overlapping event anomalies to bias anticipated centers.
        for o in h[(h.Over == 0) & (h.Out == 0)].itertuples():
            j = np.argmin(abs(ts - o.tt))
            if abs(ts[j] - o.tt) < 0.3 * per:
                ts[j] = o.tt
    return ts, dur


def runstar(k, write=False):
    t, y = readstar(k)
    if len(t) < 1000:
        return []
    group = kois[(kois.kepid == k) & (kois.koi_disposition != "FALSE POSITIVE")].copy()
    items = []
    for r in group.itertuples():
        ts, d = prior(r)
        if len(ts) < 2 or not np.isfinite(d):
            continue
        phase = (t - near(ts, t)) / d
        items.append(
            dict(
                r=r,
                ts=ts,
                d=d,
                phase=phase,
                name=r.kepler_name if isinstance(r.kepler_name, str) else r.kepoi_name,
            )
        )
    if len(items) < 2:
        return []
    # low-frequency trend estimated only from cadences away from every predicted transit, to avoid detrending a bump
    out = np.all([abs(p["phase"]) > 0.72 for p in items], axis=0)
    if out.sum() < 40:
        return []
    for rng in np.split(np.arange(len(t)), np.where(np.diff(t) > 0.5)[0] + 1):
        if len(rng) < 15:
            continue
        j = rng[out[rng]]
        if len(j) < 5:
            continue
        trend = median_filter(y[j], size=min(37, len(j) // 2 * 2 + 1), mode="nearest")
        y[rng] -= np.interp(t[rng], t[j], trend)
    noise = np.nanmedian(abs(np.diff(y[out]))) / 0.954
    if noise < 1:
        return []
    for p in items:
        good = (np.abs(p["phase"]) < 0.85) & np.all(
            [abs(a["phase"]) > 0.70 for a in items if a is not p], axis=0
        )
        dt = p["phase"][good]
        yy = y[good]
        dx = max(0.025, min(0.13, 2.0 / max(30, len(dt)) ** 0.55))
        xx = np.arange(-0.7, 0.71, dx)
        med = np.array(
            [np.median(yy[abs(dt - z) < dx]) if sum(abs(dt - z) < dx) >= 4 else np.nan for z in xx]
        )
        finite = np.isfinite(med)
        if finite.sum() < 4:
            p["model"] = lambda q: np.zeros_like(q)
            continue
        sh = gaussian_filter1d(np.interp(xx, xx[finite], med[finite]), 0.65)
        # fixed empirical mean with photocenter asymmetry suppressed. stop extrapolation beyond physical contact.
        sh = (sh + np.interp(-xx, xx, sh)) / 2
        sh = np.minimum(sh, 0)
        sh[abs(xx) > 0.62] = 0
        p["model"] = lambda q, xx=xx, sh=sh: np.interp(q, xx, sh, left=0, right=0)
    out = []
    for a, b in itertools.combinations(items, 2):
        for ta in a["ts"]:
            tb = near(b["ts"], ta)
            dmax = max(a["d"], b["d"])
            if abs(ta - tb) > (a["d"] + b["d"]) * 0.56:
                continue
            cen = (ta + tb) / 2
            mask = abs(t - cen) < max(0.46, dmax * 2.15)
            if sum(mask) < 28:
                continue
            x = t[mask]
            v = y[mask].copy()
            guard = (
                np.all(
                    [abs(p["phase"][mask]) > 0.64 for p in items if p is not a and p is not b],
                    axis=0,
                )
                if len(items) > 2
                else np.ones(len(x), bool)
            )
            x = x[guard]
            v = v[guard]
            if (
                len(x) < 23
                or sum((abs(x - ta) < a["d"] * 0.48) & (abs(x - tb) < b["d"] * 0.48)) < 1
            ):
                continue

            # fit both centers and local baseline (marginalize actual offset).
            def resid(z):
                model = (
                    a["model"]((x - ta - z[0]) / a["d"])
                    + b["model"]((x - tb - z[1]) / b["d"])
                    + z[2]
                    + z[3] * (x - cen)
                )
                return (v - model) / noise

            da = min(a["d"] * 0.30, 0.10)
            db = min(b["d"] * 0.30, 0.10)
            fit = least_squares(
                resid,
                [0, 0, 0, 0],
                bounds=([-da, -db, -np.inf, -np.inf], [da, db, np.inf, np.inf]),
                loss="soft_l1",
                x_scale=[0.012, 0.012, 100, 100],
            )
            res = resid(fit.x) * noise
            off = (abs(x - ta) > a["d"] * 0.66) & (abs(x - tb) > b["d"] * 0.66)
            if sum(off) < 8:
                continue
            s = max(noise, np.median(abs(res[off] - np.median(res[off]))) * 1.4826)
            ok = (abs(x - ta - fit.x[0]) < a["d"] * 0.54) & (abs(x - tb - fit.x[1]) < b["d"] * 0.54)
            best = (-np.inf, 0, 0, 0)
            for i in np.where(ok)[0]:
                for j in range(i + 1, min(i + 9, len(x) + 1)):
                    if not np.all(ok[i:j]) or x[j - 1] - x[i] > (j - i) * 0.023:
                        break
                    sig = res[i:j].sum() / (s * np.sqrt(j - i))
                    best = max(best, (sig, j - i, np.mean(x[i:j]), np.mean(res[i:j])))
            if not np.isfinite(best[0]):
                continue
            score, n, tm, amp = best
            row = dict(
                k=k,
                a=a["name"],
                b=b["name"],
                ta=ta,
                tb=tb,
                dt1=fit.x[0],
                dt2=fit.x[1],
                snr=score,
                n=n,
                tm=tm,
                amp=amp,
                s=s,
            )
            out.append(row)
            if write and score > 3.5:
                Path("events").mkdir(exist_ok=True)
                pd.DataFrame(dict(x=x, v=v, m=v - res, res=res, ok=ok)).to_csv(
                    "events/%s_%.2f.csv" % (k, cen), index=False
                )
    return out


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor

    ks = pairs[pairs.snr > 6].k.drop_duplicates()
    with ProcessPoolExecutor(5) as pool:
        data = [q for rows in pool.map(runstar, ks) for q in rows]
    df = pd.DataFrame(data).sort_values("snr", ascending=False)
    df.to_csv("screen.csv", index=False, float_format="%.5f")
    print(df.head(15).round(3).to_string(index=False))

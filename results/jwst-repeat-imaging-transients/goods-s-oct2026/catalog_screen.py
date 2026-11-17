"""Compare public JWST stage-3 aperture catalogs for the October 2025/2026 GOODS-S visits."""

from pathlib import Path
from glob import glob
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

bands = ["115", "150", "200", "277", "356", "444"]
RA = 53.16
DEC = -27.79


def project(df):
    return (
        np.c_[
            (df["sky_centroid.ra"] - RA) * np.cos(np.deg2rad(DEC)),
            df["sky_centroid.dec"] - DEC,
        ]
        * 3600
    )


def readcat(pid, obs, band):
    path = next(Path(str(pid)).glob(f"*-o{obs:03}_*f{band}w_cat.ecsv*"))
    cat = pd.read_csv(
        path, sep=r"\s+", comment="#", na_values=["nan", "NaN", "--", '""']
    )
    return cat.dropna(
        subset=[
            "sky_centroid.ra",
            "sky_centroid.dec",
            "aper_total_flux",
            "aper_total_flux_err",
        ]
    ).reset_index(drop=True)


def match(x, y, r):
    tree = cKDTree(y)
    sep, j = tree.query(x)
    _, ix = cKDTree(x).query(y)
    return (sep < r) & (ix[j] == np.arange(len(x))), j


if __name__ == "__main__":
    allrows = []
    stats = []
    for obs in [1, 2]:
        for band in bands:
            a = readcat(8060, obs, band)
            b = readcat(12577, obs, band)
            x = project(b)
            y = project(a)
            use, j = match(x, y, 1.0)
            consistent = use & (
                np.abs(
                    np.log(
                        np.maximum(1e-16, b.aper_total_flux)
                        / np.maximum(1e-16, a.aper_total_flux.values[j])
                    )
                )
                < 0.3
            )
            off = np.median((x - y[j])[consistent], axis=0)
            x -= off
            use, j = match(x, y, 0.12)
            new = b[use].copy().reset_index(drop=True)
            old = a.iloc[j[use]].reset_index(drop=True)
            fn, fo = [c.aper_total_flux.values * 1e9 for c in [new, old]]
            en, eo = [c.aper_total_flux_err.values * 1e9 for c in [new, old]]
            dist = np.linalg.norm((x - y[j])[use], axis=1)
            comb = np.hypot(en, eo)
            flux = (fn + fo) / 2
            mag = np.full(len(fn), np.nan)
            np.log10(np.maximum(1e-20, flux) / 3631e9, out=mag)
            mag *= -2.5
            bins = np.arange(16, 34)
            noise = np.full_like(flux, np.nan)
            bias = noise.copy()
            for lo, hi in zip(bins, bins[1:]):
                g = (mag >= lo) & (mag < hi)
                h = (abs(mag - (lo + 0.5)) < 1) & (comb > 0)
                if h.sum() < 25:
                    h = (abs(mag - (lo + 0.5)) < 2.5) & (comb > 0)
                if h.sum() < 12:
                    continue
                z = (fn[h] - fo[h]) / comb[h]
                med = np.median(z)
                bias[g] = med
                noise[g] = max(1.0, 1.4826 * np.median(np.abs(z - med)))
            z = ((fn - fo) / comb - bias) / noise
            for i in np.flatnonzero(np.abs(z) > 3):
                allrows.append(
                    dict(
                        obs=obs,
                        filter="F" + band + "W",
                        ra=new.loc[i, "sky_centroid.ra"],
                        dec=new.loc[i, "sky_centroid.dec"],
                        old_njy=fo[i],
                        new_njy=fn[i],
                        old_err_njy=eo[i],
                        new_err_njy=en[i],
                        empirical_z=z[i],
                        separation_arcsec=dist[i],
                        old_extended=old.loc[i, "is_extended"],
                        new_extended=new.loc[i, "is_extended"],
                        new_label=new.loc[i, "label"],
                        old_label=old.loc[i, "label"],
                    )
                )
            stats.append(
                dict(
                    obs=obs,
                    filter="F" + band + "W",
                    old=len(a),
                    new=len(b),
                    mutual=int(use.sum()),
                    offset_east_arcsec=off[0],
                    offset_north_arcsec=off[1],
                    median_resid_arcsec=np.median(dist),
                    median_noise_multiplier=np.nanmedian(noise),
                )
            )
    rows = pd.DataFrame(allrows)
    rows.to_csv("catalog_outliers.csv", index=False)
    pd.DataFrame(stats).to_csv("catalog_stats.csv", index=False)
    xy = np.c_[(rows.ra - RA) * np.cos(np.deg2rad(DEC)), rows.dec - DEC] * 3600
    tree = cKDTree(xy)
    used = set()
    idx = []
    for i in np.argsort(-abs(rows.empirical_z.values)):
        if i in used:
            continue
        ix = tree.query_ball_point(xy[i], 0.16)
        group = rows.iloc[ix]
        good = group[
            (group.empirical_z * np.sign(rows.empirical_z[i]) > 5)
            & (
                abs(group.new_njy - group.old_njy)
                > 0.16 * np.maximum(group.new_njy, group.old_njy)
            )
        ]
        if good["filter"].nunique() >= 2 and max(abs(good.empirical_z)) >= 7:
            used.update(ix)
            idx.extend((i, j) for j in good.index)
    peaks = rows.loc[[q[1] for q in idx]].copy()
    peaks["group"] = [q[0] for q in idx]
    peaks.to_csv("multiband.csv", index=False)
    print(len(stats), len(rows), peaks.group.nunique(), flush=True)
    print(
        peaks.groupby("group")
        .agg({"ra": "mean", "dec": "mean", "empirical_z": "median", "filter": "count"})
        .sort_values("empirical_z", key=abs, ascending=False)
        .head(15)
    )

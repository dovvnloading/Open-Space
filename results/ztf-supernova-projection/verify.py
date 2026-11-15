"""Replay the all-sky BTS screen, source localization and spectral template check."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from pixels import measure

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "inputs"
IDS = ("ZTF19aakrnqh", "ZTF22aamtomc")


def decimal_angle(s):
    h, m, sec = map(float, s.split(":"))
    return (-1 if s.startswith("-") else 1) * (abs(h) + m / 60 + sec / 3600)


def unit_vectors(ra, dec):
    a, b = np.radians(ra), np.radians(dec)
    return np.column_stack((np.cos(a) * np.cos(b), np.sin(a) * np.cos(b), np.sin(b)))


def arcsec(first, second):
    v = unit_vectors([first[0], second[0]], [first[1], second[1]])
    return np.degrees(2 * np.arcsin(np.linalg.norm(v[1] - v[0]) / 2)) * 3600


def catalog_screen():
    catalog = pd.read_csv(DATA / "bts.csv")
    catalog["redshift"] = pd.to_numeric(catalog.redshift, errors="coerce")
    sn = catalog[catalog.type.str.startswith("SN") & catalog.IAUID.ne("-")].dropna(subset=["redshift"])
    # If BTS split one explosion across multiple IDs, keep the shortest separation
    # between rows that represent the same pair of distinct names.
    vectors = unit_vectors(sn.RA.map(decimal_angle) * 15, sn.Dec.map(decimal_angle))
    tree = cKDTree(vectors)
    pairs = {}
    for i, j in sorted(tree.query_pairs(2 * np.sin(np.radians(100 / 3600) / 2))):
        left, right = sn.iloc[i], sn.iloc[j]
        if left.IAUID == right.IAUID:
            continue
        key = tuple(sorted((left.IAUID, right.IAUID)))
        separation = np.degrees(2 * np.arcsin(np.linalg.norm(vectors[i] - vectors[j]) / 2)) * 3600
        row = (separation, abs(left.redshift - right.redshift))
        if key not in pairs or row[0] < pairs[key][0]:
            pairs[key] = row
    print(f"{len(catalog)} BTS rows, {len(sn)} eligible rows, {sn.IAUID.nunique()} names")
    target = ("SN2019bns", "SN2022lxq")
    for radius in (4, 10, 30, 60, 100):
        subset = sorted(((v[1], k) for k, v in pairs.items() if v[0] <= radius), reverse=True)
        assert subset[0][1] == target
        print(f"{radius:3d} arcsec: {len(subset)} distinct pairs, highest Δz {subset[0][0]:.3f}, next {subset[1][0]:.3f}")
    print(f"Target separation {pairs[target][0]:.4f} arcsec")
    matches = catalog.set_index("ZTFID").loc[list(IDS)]
    print(f"BTS peak-epoch separation {np.ptp(matches.peakt):.2f} days")
    ia = sn[sn.type.eq("SN Ia")].sort_values("peakmag").drop_duplicates("IAUID")
    far = ia[ia.redshift.between(0.09, 0.13)]
    print(f"{len(far)} normal Ia names at z=0.09–0.13, brightest {far.peakmag.min():.4f} mag")


def lightcurves():
    alerts = []
    for obj in IDS:
        raw = json.loads((DATA / f"{obj}.json").read_text())
        alerts.append(pd.DataFrame(raw["detections"]).query("isdiffpos > 0 and fid in [1, 2] and sigmapsf < 0.17")
                      .sort_values("sigmapsf").drop_duplicates(["pid", "fid"]))
    positions = [d[["ra", "dec"]].median().tolist() for d in alerts]
    host = pd.read_csv(DATA / "globalhost_data.csv").iloc[0]
    print("Distinct usable g/r alerts:", *map(len, alerts))
    print("ZTF centroid separation: %.3f arcsec" % arcsec(*positions))
    print("DR2 host offsets in arcsec:", [round(arcsec(pos, (host.ra_host, host.dec_host)), 3) for pos in positions])
    print("Alert r-band brightest points:", [d[d.fid.eq(2)].magpsf.min() for d in alerts])
    forced = pd.read_csv(DATA / "dr2-lightcurve.txt", sep=r"\s+", comment="#")
    good = forced[forced.flag.eq(0)].copy()
    good["snr"] = good.flux / good.flux_err
    for lo, hi in [(59660, 59725), (59743, 59765)]:
        print(f"DR2 MJD {lo}–{hi}")
        print(good[good.mjd.between(lo, hi)].groupby("filter").agg(n=("flux", "size"), median=("flux", "median"),
                                                                    min_snr=("snr", "min"), max_snr=("snr", "max")).round(2))
    for fid, band in [(1, "ztfg"), (2, "ztfr")]:
        aperture = good[good.mjd.between(59733, 59785) & good["filter"].eq(band)].sort_values("mjd")
        # Add any detections below the centroid quality cut for the exposure match.
        raw = pd.DataFrame(json.loads((DATA / f"{IDS[1]}.json").read_text())["detections"])
        source = raw[raw.fid.eq(fid) & raw.isdiffpos.gt(0)].sort_values("sigmapsf").drop_duplicates("pid").sort_values("mjd")
        matched = pd.merge_asof(aperture, source[["mjd", "magpsf"]], on="mjd", direction="nearest",
                                tolerance=0.00005).dropna(subset=["magpsf"])
        target_flux = 10 ** (0.4 * (matched.ZP - matched.magpsf))
        print(f"{band}: {len(matched)} matched frames, median flux fraction {np.median(matched.flux / target_flux):.3f},"
              f" Spearman rho {spearmanr(matched.flux, target_flux).statistic:.3f}")
    gaia = pd.read_csv(DATA / "Gaia22cvp.csv", skiprows=1)
    numeric = pd.to_numeric(gaia.averagemag, errors="coerce").dropna()
    print(f"Gaia G: {len(numeric)} numeric magnitudes, range {numeric.min():.2f}–{numeric.max():.2f}")


def spectral_test():
    x, y, variance = np.loadtxt(next(DATA.glob("*.ascii"))).T
    templates = pd.read_csv(DATA / "hsiao-window.csv")
    clean = np.isfinite(y) & (y > 0) & (variance > 0)
    clean &= ~((x > 7550) & (x < 7700)) & ~((x > 6820) & (x < 6940))
    for low, high in [(4300, 8200), (4300, 5800), (5800, 8000)]:
        use = clean & (x >= low) & (x <= high)
        obs = y[use] / np.median(y[use])
        wave = x[use]
        weights = np.sqrt(np.median(variance[use]) / variance[use])
        scaled_wave = (wave - 6000) / 1000
        scores = []
        for phase in range(3, 15):
            template = templates[f"p{phase}"].to_numpy()
            for z in np.arange(181) / 1000:
                shape = np.interp(wave / (1 + z), templates.wavelength, template)
                shape /= np.median(shape)
                basis = np.column_stack((shape, shape * scaled_wave, shape * scaled_wave**2,
                                         np.ones_like(wave), scaled_wave))
                coeff = np.linalg.lstsq(basis * weights[:, None], obs * weights, rcond=None)[0]
                residual = (obs - basis @ coeff) * weights
                scores.append((residual @ residual, z, phase))
        rss, z, phase = min(scores)
        low_rss = min(s[0] for s in scores if 0.006 <= s[1] <= 0.036)
        print(f"{low}–{high} Angstrom: z={z:.3f}, phase +{phase} days, low-z / minimum RSS {low_rss/rss:.3f}")


if __name__ == "__main__":
    for name, entry in json.loads((ROOT / "provenance.json").read_text()).items():
        if hashlib.sha256((DATA / name).read_bytes()).hexdigest() != entry["file_sha256"]:
            raise RuntimeError(f"Snapshot mismatch: {name}")
    catalog_screen()
    lightcurves()
    measure()
    spectral_test()

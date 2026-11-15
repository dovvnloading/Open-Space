"""Screen DESI DR1 repeat DA exposures, then fit individual Balmer lines.

Run from any directory with Python 3.10+, numpy, pandas, scipy and astropy.
Products are written to --output. Raw inputs are downloaded once into --cache.
"""
import argparse
import hashlib
import io
import itertools
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
C_KMS = 299792.458
LINES = {"alpha": 6564.61, "beta": 4862.68, "gamma": 4341.68, "delta": 4102.89}
IDENTITY = ["WDJname", "NIGHT", "EXPOSURE", "FIBER"]
CONTROL = "WDJ101606.87-011917.14"
TARGETS = [CONTROL, "WDJ042155.26-132531.35", "WDJ023809.44+273329.79"]
ARMS = ["B", "R"]


def retrieve(cache):
    sources = json.loads((ROOT / "sources.json").read_text())
    cache.mkdir(parents=True, exist_ok=True)
    for name, source in sources.items():
        path = cache / name
        if not path.exists():
            partial = path.with_suffix(path.suffix + ".part")
            urllib.request.urlretrieve(source["url"], partial)
            partial.replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Unexpected bytes for {name}")


def screen(catalogue, rv, exposures):
    stars = catalogue[(catalogue.specType == "DA") &
                      (catalogue.fluxContamination_DESI < .05) & (catalogue.nExposures >= 2)]
    measurements = rv[(rv.fitType == "joint") & (rv.model == "H") &
                      rv.EXPOSURE.notna() & rv.WDJname.isin(stars.WDJname)]
    # Exclude any ambiguous exposure keys instead of multiplying catalogue rows.
    measurements = measurements.drop_duplicates(IDENTITY, keep=False)
    exposures = exposures.drop_duplicates(IDENTITY, keep=False)
    measurements = measurements.merge(exposures, on=IDENTITY, validate="one_to_one")
    keep = (measurements.SN_B > 6) & (measurements.SN_R > 4)
    for arm in ARMS:
        keep &= (measurements[f"RV_{arm}_error"].between(.5, 70) &
                 (measurements[f"RV_{arm}"].abs() < 1500) &
                 (measurements[f"fractionMasked_{arm}"] < .05))
    good = measurements.loc[keep]
    results, pairs_searched, paired_stars = [], 0, 0
    for name, series in good.groupby("WDJname"):
        if len(series) < 2:
            continue
        paired_stars += 1
        best = None
        for i, j in itertools.combinations(series.itertuples(), 2):
            if i.EXPOSURE == j.EXPOSURE:
                continue
            pairs_searched += 1
            db, dr = i.RV_B - j.RV_B, i.RV_R - j.RV_R
            eb = np.sqrt(i.RV_B_error**2 + j.RV_B_error**2 + 2 * 8**2)
            er = np.sqrt(i.RV_R_error**2 + j.RV_R_error**2 + 2 * 8**2)
            if db * dr <= 0 or abs(db - dr) >= 3 * np.hypot(eb, er):
                continue
            score = min(abs(db) / eb, abs(dr) / er)
            if best is None or score > best["score"]:
                best = dict(name=name, n=len(series), first=int(i.EXPOSURE), second=int(j.EXPOSURE),
                            score=score, blue_delta=db, red_delta=dr, blue_error=eb, red_error=er)
        if best is not None:
            results.append(best)
    ranks = pd.DataFrame(results).sort_values("score", ascending=False).reset_index(drop=True)
    coverage = dict(catalogue_rows=len(catalogue), initially_selected=len(stars),
                    quality_exposures=len(good), stars_with_two_exposures=paired_stars,
                    pairs_searched=pairs_searched, stars_with_concordant_pair=len(ranks))
    return ranks, good, coverage


def fit_line(data, centre, window):
    wavelength, flux, _, ivar, mask = data.T
    use = (abs(wavelength - centre) < window) & (ivar > 0) & (mask == 0) & np.isfinite(flux)
    x = wavelength[use] - centre
    scale = np.median(flux[use][abs(x) > .65 * window])
    y = flux[use] / scale
    error = ivar[use] ** -.5 / scale

    def residual(parameters):
        level, slope, shift, depth1, width1, depth2, width2 = parameters
        model = level + slope * x / window
        model -= depth1 / (1 + ((x - shift) / width1)**2)
        model -= depth2 / (1 + ((x - shift) / width2)**2)
        return (model - y) / error

    fit = least_squares(residual, [1, 0, 0, .15, 1.7, .3, 13],
                        bounds=([.3, -.6, -18, 0, .4, 0, 5], [2.9, .6, 18, 2.2, 8, 2.2, 100]),
                        xtol=1e-10)
    if not fit.success:
        raise ValueError(f"Line fit failed at {centre}")
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
    return (fit.x[2] * C_KMS / centre, np.sqrt(covariance[2, 2]) * C_KMS / centre,
            np.sqrt(np.mean(fit.fun**2)))


def measure(cache):
    rows = []
    for name in TARGETS:
        with tarfile.open(cache / (name + ".tar.gz")) as archive:
            for member in archive:
                fields = member.name.rsplit("_", 3)
                if len(fields) != 4 or not fields[1].isdigit() or not fields[-1].endswith(("-b.dat", "-r.dat")):
                    continue
                data = np.loadtxt(io.BytesIO(archive.extractfile(member).read()))
                for line, centre in LINES.items():
                    if (line == "alpha") != member.name.endswith("-r.dat"):
                        continue
                    for window in (35, 50, 65):
                        velocity, error, rms = fit_line(data, centre, window)
                        rows.append((name, int(fields[1]), int(fields[2]), line, window, velocity, error, rms))
    return pd.DataFrame(rows, columns=["name","night","exposure","line","window_A","v_kms","error_kms","resid_rms"])


def field_controls(catalogue, rv, selected):
    da = catalogue.loc[catalogue.specType == "DA", "WDJname"]
    q = rv[(rv.model == "H") & (rv.fitType == "joint") & rv.WDJname.isin(da)]
    coadd = q[q.EXPOSURE.isna()]
    rows = []
    for row in selected.itertuples():
        others = q[(q.EXPOSURE == row.EXPOSURE) & (q.WDJname != row.WDJname)]
        others = others.merge(coadd, on="WDJname", suffixes=("", "_coadd"))
        for arm in ARMS:
            v = "RV_" + arm
            # Single-epoch coadds share all photons and cannot control a shift.
            multi = catalogue.loc[catalogue.nExposures > 1, "WDJname"]
            ref = others[others.WDJname.isin(multi) & others[v+"_error"].between(.1, 55) &
                         others[v+"_error_coadd"].between(.1, 40)]
            differences = ref[v] - ref[v+"_coadd"]
            rows.append((row.WDJname, int(row.EXPOSURE), arm, len(ref), differences.median()))
    return pd.DataFrame(rows, columns=["target", "exposure", "arm", "other_stars", "median_epoch_minus_coadd_kms"])


def catalogue_checks(cache, catalogue):
    sample = catalogue[catalogue.WDJname.isin(TARGETS)].copy()
    munday = pd.read_csv(cache / "CloseDoubleWhiteDwarfs_20261109.csv")
    catalogue_coordinates, names = [], []
    for row in munday.to_dict("records"):
        try:
            point = SkyCoord(str(row["RA"]), str(row["Dec"]).replace("−", "-"), unit=(u.hourangle, u.deg))
        except ValueError:
            continue
        catalogue_coordinates.append(point)
        names.append(row["Unique ID"].strip())
    sky = SkyCoord(catalogue_coordinates)
    result = []
    with fits.open(cache / "sdssv.fits") as hdus:
        for row in sample.itertuples():
            s = row.WDJname[3:]
            coord = SkyCoord(s[:2]+"h"+s[2:4]+"m"+s[4:9]+"s "+s[9:12]+"d"+s[12:14]+"m"+s[14:]+"s")
            separations = coord.separation(sky).arcsec
            nearest = int(np.argmin(separations))
            gaia = int(row.designation.split()[-1])
            nmatches = np.count_nonzero(hdus[1].data["gaia_dr3_id"] == gaia)
            result.append((row.WDJname, gaia, len(sky), names[nearest], separations[nearest], nmatches))
    return pd.DataFrame(result, columns=["name","gaia_dr3","munday_positions","munday_nearest",
                                         "munday_sep_arcsec","sdssv_vac_matches"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("cache"))
    parser.add_argument("--output", type=Path, default=Path("replayed"))
    args = parser.parse_args()
    retrieve(args.cache)
    def load(table):
        return pd.read_csv(args.cache / f"DESI_DR1_{table}_1.0.csv.gz")
    catalogue, rv, exposures = (load(k) for k in ("catalogue", "RVs", "exposures"))
    ranks, good, coverage = screen(catalogue, rv, exposures)
    selected = good[good.WDJname.isin(TARGETS)].copy()
    cols = IDENTITY + ["TARGETID","MJD","BMJD","EXPTIME","SN_B","SN_R","RV_B","RV_B_error","RV_R","RV_R_error"]
    products = {"screen-top10": ranks.head(10), "selected-exposures": selected[cols],
                "balmer-lines": measure(args.cache), "field-controls": field_controls(catalogue, rv, selected),
                "prior-catalogues": catalogue_checks(args.cache, catalogue)}
    args.output.mkdir(parents=True, exist_ok=True)
    for name, table in products.items():
        table.to_csv(args.output / (name + ".csv"), index=False, float_format="%.6f")
    (args.output / "coverage.json").write_text(json.dumps(coverage, indent=2) + "\n")


if __name__ == "__main__":
    main()

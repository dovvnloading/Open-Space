"""Targeted NEOWISE positive control; run in this directory with Python, astropy, scipy, numpy and pandas."""
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from scipy.cluster.hierarchy import fclusterdata

HERE = Path(__file__).resolve().parent
RA, DEC = 133.773, -7.243
COSDEC = np.cos(np.deg2rad(DEC))
COLUMNS = "ra,dec,mjd,w2mpro,w2sigmpro,w2snr,w2rchi2,w1mpro,w1snr,sso_flg,cc_flags,scan_id,frame_num,sigra,sigdec,sigradec,qual_frame,qi_fact"
QUERY = (f"SELECT {COLUMNS} FROM neowiser_p1bs_psd WHERE "
         "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',133.773,-7.243,0.019))")
URL = "https://irsa.ipac.caltech.edu/TAP/sync?" + urlencode(
    dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=QUERY))


def get_epochs(data):
    clean = data.loc[(data.w2snr >= 5) & (data.w1snr < 3) & (data.sso_flg == 0)
                     & (data.qual_frame >= 5) & (data.w2rchi2 < 3)
                     & (data.cc_flags.str[1] == "0")].sort_values("mjd").copy()
    clean["visit"] = (clean.mjd.diff() > 20).cumsum()
    clean["x"] = (clean.ra - RA) * 3600 * COSDEC
    clean["y"] = (clean.dec - DEC) * 3600
    epochs = []
    for visit, detections in clean.groupby("visit"):
        labels = fclusterdata(detections[["x", "y"]], 2.5, criterion="distance", method="single")
        groups = [g for _, g in detections.groupby(labels) if len(g) >= 3]
        if len(groups) != 1:
            raise ValueError(f"Ambiguous control field in visit {visit}: {len(groups)} clusters")
        g = groups[0]
        entry = dict(visit=visit, mjd=g.mjd.median(), n=len(g), w2=g.w2mpro.median())
        for axis in ["x", "y"]:
            values = g[axis]
            med = values.median()
            # Gaussian median standard error via MAD, plus an explicit per-visit error floor.
            sem = 1.2533 * 1.4826 * (values - med).abs().median() / np.sqrt(len(values))
            entry[axis] = med
            entry["s" + axis] = np.hypot(sem, 0.12)
        epochs.append(entry)
    return pd.DataFrame(epochs)


def solve(epochs, scan=False):
    mjd = epochs.mjd.to_numpy()
    dt = (mjd - 58500.) / 365.25
    a, d = np.deg2rad([RA, DEC])
    X, Y, Z = get_body_barycentric("earth", Time(mjd, format="mjd", scale="utc")).xyz.to_value("AU")
    px = X*np.sin(a) - Y*np.cos(a)
    py = X*np.cos(a)*np.sin(d) + Y*np.sin(a)*np.sin(d) - Z*np.cos(d)
    n = len(mjd)
    design = np.zeros((2*n, 7 if scan else 5))
    design[:n, 0] = design[n:, 1] = 1.
    design[:n, 2] = design[n:, 3] = dt
    design[:, 4] = np.r_[px, py]
    if scan:
        # WISE alternates scan orientation at this low ecliptic latitude.
        parity = np.sign(px)
        design[:n, 5] = design[n:, 6] = parity
    sigma = np.r_[epochs.sx, epochs.sy]
    values = np.r_[epochs.x, epochs.y]
    A = design / sigma[:, None]
    params = np.linalg.lstsq(A, values/sigma, rcond=None)[0]
    residual = (values - design @ params) / sigma
    chi2 = float(residual @ residual)
    dof = len(values) - design.shape[1]
    cov = np.linalg.inv(A.T @ A)*max(1., chi2/dof)
    return dict(parameters=params.tolist(), stderr=np.sqrt(np.diag(cov)).tolist(),
                covariance=cov.tolist(), chi2=chi2, dof=dof, condition=float(np.linalg.cond(A)))


def main():
    path = HERE / "single-exposures.csv.gz"
    if not path.exists():
        raw = urlopen(URL, timeout=240).read()
        if not raw.startswith(b"ra,dec,mjd,"):
            raise ValueError("IRSA did not return the requested table")
        path.write_bytes(gzip.compress(raw, mtime=0))
    raw = gzip.decompress(path.read_bytes())
    data = pd.read_csv(path, dtype={"cc_flags": str, "scan_id": str})
    epochs = get_epochs(data)
    epochs.to_csv(HERE/"epochs.csv", index=False, float_format="%.9f")
    base = solve(epochs)
    robustness = [solve(epochs.drop(index=[i]))["parameters"][4] for i in epochs.index]
    result = dict(source=URL, query=QUERY, sha256_raw=hashlib.sha256(raw).hexdigest(),
                  n_field=len(data), n_target=int(epochs.n.sum()), n_visits=len(epochs),
                  mjd_range=[float(epochs.mjd.min()), float(epochs.mjd.max())],
                  reference_mjd=58500, tangent_point_deg=[RA, DEC],
                  parameter_order=["x_arcsec","y_arcsec","pmra_cosdec_arcsec_yr",
                                   "pmdec_arcsec_yr","relative_parallax_arcsec"],
                  fit=base, fit_with_independent_scan_offsets=solve(epochs, scan=True),
                  leave_one_visit_out_parallax_range_arcsec=[min(robustness), max(robustness)])
    (HERE/"fit.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fit", base["parameters"], "formal stderr", base["stderr"], "chi2/dof", base["chi2"],base["dof"])
    nuisance=result["fit_with_independent_scan_offsets"]
    print("Scan nuisance",nuisance["parameters"][4:],nuisance["stderr"][4:],robustness[0])


if __name__ == "__main__":
    main()

"""Retrieve EROS-2 cone photometry and test eclipse contrast at three OGLE-IV SMC coordinates.

From an Open-Space checkout with numpy installed, pass a cache path, the OGLE-IV
catalog file (next to phot/I as released by OGLE) and an output JSON path.
"""

import argparse
from collections import defaultdict
import numpy as np
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
import hashlib
import json

TARGETS = {
    "07730": ("01:07:06.79", "-72:45:31.9"),
    "06703": ("00:48:53.97", "-73:13:41.8"),
    "09564": ("00:53:09.24", "-71:59:43.0"),
}


def deg(s):
    a, b, c = map(float, s.split(":"))
    return (abs(a) + b / 60 + c / 3600) * (-1 if s.startswith("-") else 1)


def url_for(target, table):
    r, d = TARGETS[target]
    q = {
        "-source": "II/390/" + table,
        "-c": f"{15 * deg(r):.7f} {deg(d):.7f}",
        "-c.rs": "8",
        "-out.max": "30000",
    }
    if table == "xcat":
        q["-out.all"] = "1"
    else:
        q["-out"] = "EROS,HJD,Rmag,e_Rmag,Bmag,e_Bmag"
    return "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(q)


def fetch(directory, target, table):
    url = url_for(target, table)
    path = directory / f"{target}_{table}.tsv"
    if not path.exists():
        data = urlopen(url, timeout=90).read()
        path.write_bytes(data)
    raw = path.read_bytes()
    return dict(
        target=target,
        table=table,
        url=url,
        bytes=len(raw),
        sha256=hashlib.sha256(raw).hexdigest(),
    )


def parse(path):
    raw = path.read_bytes()
    if b"#ERROR" in raw or b"#Name: II/390/" not in raw:
        raise ValueError(path)
    lines = [
        l for l in raw.decode().splitlines() if l.strip() and not l.startswith("#")
    ]
    if not lines[2].startswith("---"):
        raise ValueError(path)
    rows = [
        dict(zip(lines[0].split("\t"), (v.strip() for v in l.split("\t")), strict=True))
        for l in lines[3:]
    ]
    if len(rows) >= 30000:
        raise ValueError("truncated")
    return rows, hashlib.sha256(
        ("\n".join(sorted(lines[3:])) + "\n").encode()
    ).hexdigest()


def estimate(
    curve, period, t0, shift=0, error_max=0.5, width=0.05, time_zero=2450000.0
):
    t, m, e = curve.T
    ok = np.isfinite(t + m + e) & (e > 0) & (e <= error_max)
    t, m, e = t[ok] + time_zero, m[ok], np.hypot(e[ok], 0.03)
    order = np.argsort(t, kind="stable")
    t, m, e = t[order], m[order], e[order]
    phase = ((t - t0) / period - shift) % 1
    distance = abs((phase + 0.25) % 0.5 - 0.25)
    use = (distance < width) | (distance >= 0.13)
    inside = distance[use] < width
    t, m, e = t[use], m[use], e[use]
    out = {"good_in": int(inside.sum()), "good_out": int((~inside).sum())}
    if inside.sum() < 12 or (~inside).sum() < 30:
        return out
    season = (t // 365.25).astype(int)
    X = np.column_stack([inside] + [season == i for i in np.unique(season)])
    # Huber iterative fit, downweighting gross photometric outliers.
    base = 1 / e**2
    w = base.copy()
    for i in range(8):
        beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], m * np.sqrt(w), rcond=None)[0]
        resid = m - X @ beta
        scale = max(1, 1.4826 * np.median(abs(resid / e - np.median(resid / e))))
        w = base * np.minimum(1.0, 2.5 * scale / np.maximum(abs(resid / e), 1e-12))
    # Seven-day cluster-sandwich SE with Huber scores and derivative of their estimating equations.
    cutoff = 2.5 * scale
    deriv = base * (abs(resid / e) <= cutoff)
    bread = np.linalg.pinv(X.T @ (deriv[:, None] * X))
    groups = (t // 7).astype(int)
    scores = X * (w * resid)[:, None]
    block = np.stack([scores[groups == i].sum(axis=0) for i in np.unique(groups)])
    n, k = X.shape
    g = len(block)
    variance = bread @ (block.T @ block) @ bread * g / (g - 1) * (n - 1) / (n - k)
    out.update(
        contrast=float(beta[0]),
        se=float(np.sqrt(max(0, variance[0, 0]))),
        blocks=g,
        first_hjd=float(t.min()),
        last_hjd=float(t.max()),
        median_error=float(np.median(e)),
    )
    return out


def main(raw, catalog):
    cat = {
        l.split()[0][-5:]: l.split()
        for l in catalog.read_text().splitlines()
        if l.split()[0][-5:] in TARGETS
    }
    out = {
        "method": "Locked OGLE orbital ephemeris. Fit a common primary/secondary core contrast (phase within 0.05 of 0 or 0.5) versus baseline (at least 0.13 from either core), with year-bin intercepts, Huber scores, positive errors <=0.5 mag and 0.03 mag floor. Seven-day block sandwich standard errors. Sensitivities and time-held-out phase scans are diagnostic only.",
        "time_reference": "EROS HJD = table date + 2450000.0. OGLE-IV ecl times are full HJD.",
        "catalog_source": "https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/ecl.dat",
        "catalog_sha256": hashlib.sha256(catalog.read_bytes()).hexdigest(),
        "targets": {},
    }
    for name in TARGETS:
        period, t0 = map(float, cat[name][3:5])
        roster, rhash = parse(raw / f"{name}_xcat.tsv")
        ts, thash = parse(raw / f"{name}_xtime.tsv")
        grouped = defaultdict(list)
        for r in ts:
            grouped[r["EROS"]].append(
                [
                    float(r[a] or "nan")
                    for a in ["HJD", "Rmag", "e_Rmag", "Bmag", "e_Bmag"]
                ]
            )
        assert {c["EROS"] for c in roster} == set(grouped)
        target = out["targets"][name] = {
            "period": period,
            "epoch": t0,
            "total_records": len(roster),
            "total_exposures": len(ts),
            "catalog_data_sha256": rhash,
            "photometry_data_sha256": thash,
            "sources": [],
        }
        ogle_path = catalog.parent / "phot" / "I" / ("OGLE-SMC-ECL-" + name + ".dat")
        ogle = np.loadtxt(ogle_path)
        target["ogle_i_lightcurve_sha256"] = hashlib.sha256(
            ogle_path.read_bytes()
        ).hexdigest()
        edges = [2455000, 2457000, 2458300, 2459700, 2500000]
        target["ogle_i_same_statistic"] = [
            {
                "start_hjd": lo,
                "end_hjd": hi,
                **estimate(
                    ogle[(ogle[:, 0] >= lo) & (ogle[:, 0] < hi)],
                    period,
                    t0,
                    time_zero=0,
                    error_max=0.1,
                ),
            }
            for lo, hi in zip(edges, edges[1:])
        ]
        for c in sorted(roster, key=lambda r: float(r["_r"])):
            arr = np.array(grouped[c["EROS"]])
            assert len(arr) == int(c["Nlc"])
            s = {
                "id": c["EROS"],
                "separation_arcsec": float(c["_r"]),
                "ra_deg": float(c["RAJ2000"]),
                "dec_deg": float(c["DEJ2000"]),
                "raw": len(arr),
            }
            target["sources"].append(s)
            for i, band in enumerate("RB"):
                curve = arr[:, [0, 1 + 2 * i, 2 + 2 * i]]
                s[band] = estimate(curve, period, t0)
                if s["separation_arcsec"] < 2:
                    s[band]["sensitivities"] = {
                        k: estimate(curve, period, t0, **opts)
                        for k, opts in [
                            ("phase_m0.04", {"shift": -0.04}),
                            ("phase_p0.04", {"shift": 0.04}),
                            ("core_0.035", {"width": 0.035}),
                            ("core_0.07", {"width": 0.07}),
                            ("error_max_1", {"error_max": 1}),
                        ]
                    }
            if s["separation_arcsec"] < 2:
                # Select the offset using early simultaneous red/blue EROS data only.
                # Use summed standardized depths just as a ranking statistic, never as a z score.
                split = np.median(arr[:, 0])
                tr = arr[arr[:, 0] <= split]
                te = arr[arr[:, 0] > split]
                candidates = []
                for shift in np.arange(-0.25, 0.25, 0.005):
                    values = [
                        estimate(
                            tr[:, [0, 1 + 2 * i, 2 + 2 * i]], period, t0, shift=shift
                        )
                        for i in range(2)
                    ]
                    if all("se" in v and v["se"] > 0 for v in values):
                        candidates.append(
                            (sum(v["contrast"] / v["se"] for v in values), shift)
                        )
                if candidates:
                    _, frozen = max(candidates)
                    s["time_holdout"] = {
                        "train_end_hjd": float(split + 2450000.0),
                        "train_selected_shift": float(frozen),
                        "shift_trials_measured": len(candidates),
                        "late": {
                            b: estimate(
                                te[:, [0, 1 + 2 * i, 2 + 2 * i]],
                                period,
                                t0,
                                shift=frozen,
                            )
                            for i, b in enumerate("RB")
                        },
                    }
    return out


def check_clock(cache):
    key = "1316.9098..1316.9100"
    params = {
        "-source": "II/390/xtinfo",
        "HJD": key,
        "-out.max": 50,
        "-out": "FileName,MJD,HJD,Date",
    }
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(params)
    path = cache / "epoch_check_0.tsv"
    if not path.exists():
        path.write_bytes(urlopen(url, timeout=60).read())
    rows, digest = parse(path)
    row = next(r for r in rows if r["FileName"] == "sm00300trr9e17119.fits")
    offset = float(row["HJD"]) + 2450000.0 - (float(row["MJD"]) + 2400000.5)
    if abs(offset) > 0.006:
        raise ValueError(
            "EROS and image observation clocks disagree by more than 0.006 d"
        )
    return {
        "filename": row["FileName"],
        "mjd": float(row["MJD"]),
        "eros_hjd": float(row["HJD"]),
        "heliocentric_minus_geocentric_days": offset,
        "url": url,
        "data_sha256": digest,
        "source": "https://doi.org/10.1051/0004-6361/202557321 (Appendix A.4)",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("ogle_catalog", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    receipts = [
        fetch(args.cache, target, table)
        for target in TARGETS
        for table in ("xcat", "xtime")
    ]
    result = main(args.cache, args.ogle_catalog)
    result["receipts"] = receipts
    result["clock_check"] = check_clock(args.cache)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

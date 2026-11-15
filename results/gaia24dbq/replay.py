"""Replay the eRO-ExTra / Gaia alerts match and the photometry of Gaia24dbq.

Requires Python 3.10+ with NumPy and SciPy. Run without options on the
bundled evidence or use --live to repeat the archive queries and sky search.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import statistics
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ERO = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2bA/693/A62/catalog&-out.all&-out.max=400"
ALERTS = "https://gsaweb.ast.cam.ac.uk/alerts/alerts.csv"
LIGHT = "https://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia24dbq/lightcurve.csv"
TAP = "https://datalab.noirlab.edu/tap/sync?"
TARGET = "1eRASS J111110.6-434214"


def fetch(url, hashes, tag):
    with urlopen(url, timeout=120) as response:
        raw = response.read()
    hashes[tag] = hashlib.sha256(raw).hexdigest()
    return raw.decode("utf-8-sig")


def tap(query, hashes, tag):
    url = TAP + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    return list(csv.DictReader(fetch(url, hashes, tag).splitlines()))


def coordinates(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.column_stack((np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)))


def screen(ero, alerts):
    ra = np.array([float(row["RALdeg"]) for row in ero])
    dec = np.array([float(row["DELdeg"]) for row in ero])
    alert_vec = coordinates([float(r["RaDeg"]) for r in alerts], [float(r["DecDeg"]) for r in alerts])
    tree = cKDTree(alert_vec)
    scale = 180*3600/math.pi
    radius = 2*math.sin(8/(2*scale))
    matches, controls = [], {}
    # Equal angular east-west shifts along declination circles to first order.
    for offset in [0, -600, -300, -120, -60, 60, 120, 300, 600]:
        vec = coordinates(ra + offset/(3600*np.cos(np.radians(dec))), dec)
        lists = tree.query_ball_point(vec, radius)
        if offset:
            controls[str(offset)] = sum(map(len, lists))
            continue
        for i, indices in enumerate(lists):
            for j in indices:
                r, a = ero[i], alerts[j]
                sep = 2*math.asin(np.linalg.norm(vec[i]-alert_vec[j])/2)*scale
                matches.append(dict(xray=r["ERO"], catalog_tns=r["TNS"], alert=a["#Name"],
                                    alert_ra=float(a["RaDeg"]), alert_dec=float(a["DecDeg"]),
                                    sep_arcsec=round(sep, 4)))
    return matches, controls


def live():
    hashes = {}
    table = fetch(ERO, hashes, "ero_tsv_with_dynamic_header")
    lines = [x for x in table.splitlines() if x and not x.startswith("#")]
    ero = [{k:v.strip() for k,v in r.items()} for r in csv.DictReader(lines, delimiter="\t")
           if r["recno"].strip().isdigit()]
    raw = fetch(ALERTS, hashes, "gaia_alerts")
    alerts = [{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(raw.splitlines())]
    matches, controls = screen(ero, alerts)
    r = next(r for r in ero if r["ERO"] == TARGET)
    raw = fetch(LIGHT, hashes, "gaia_photometry")
    g = [[row[0],float(row[1]),row[2]] for row in csv.reader(raw.splitlines()[2:]) if len(row) == 3]
    locs = tap("select objid,ra,dec,flux_w1,flux_w2 from ls_dr11.tractor "
               "where brickid=102468 and ra between 167.770 and 167.819 "
               "and dec between -43.722 and -43.684 order by objid", hashes, "tractor_field")
    # Fix proximity and baseline brightness before looking at light curves.
    target = next(x for x in locs if x["objid"] == "5969")
    def separation(row):
        a, d = float(row["ra"]),float(row["dec"])
        return math.hypot((a-float(target["ra"]))*math.cos(math.radians(d)),d-float(target["dec"]))*3600
    neighbors = sorted((x for x in locs if x["objid"]!="5969" and separation(x)<60
                        and float(x["flux_w1"])>30 and float(x["flux_w2"])>20), key=separation)
    ids = [target["objid"]] + [x["objid"] for x in neighbors]
    sql = "select * from ls_dr11.wise where brickid=102468 and objid in (" + ",".join(ids) + ") order by objid"
    w = tap(sql, hashes, "unwise_field")
    epochs = []
    for source in w:
        for band in ("w1","w2"):
            for i in range(1, 1 + sum(key.startswith("lc_mjd_w1_") for key in source)):
                def val(key):
                    return float(source[f"lc_{key}_{band}_{i}"])
                if val("nobs") <= 0 or val("flux_ivar") <= 0:
                    continue
                epochs.append([int(source["objid"]), band, val("mjd"), val("flux"),
                               val("flux_ivar"),val("nobs"),val("fracflux"),val("rchisq")])
    fields = ("ERO RA_ICRS DE_ICRS ePos RALdeg DELdeg FeRASS1 FeRASS2 FeRASS3 "
              "FeRASS4 FeRASS5 e_FeRASS1 e_FeRASS2 e_FeRASS3 e_FeRASS4 e_FeRASS5 "
              "ERODate Amp Signi Gamma e_Gamma LCClass z zType pany pi TNS MType").split()
    return dict(n_ero=len(ero), n_gaia=len(alerts), matches=matches, control_counts=controls,
                ero={k:r[k] for k in fields}, gaia=g, wise_epochs=epochs,
                comparison_sources=[{**x,"offset_arcsec":round(separation(x),3)} for x in neighbors],
                digests=hashes)


def summarize(data):
    curve = [(d,float(m)) for d,t,m in data["gaia"] if m not in ("null","untrusted","")]
    # First give each UT date equal weight to avoid weighting Gaia's scanning clusters.
    nights = {}
    for d, m in curve:
        nights.setdefault(d[:10], []).append(m)
    nights = {d:statistics.median(m) for d,m in nights.items()}
    windows = {}
    for first, last in [("2015","2018"),("2019-01","2019-06"),("2020-01","2021-01"),
                        ("2023-01","2024-01"),("2024-06","2024-11"),("2024-10","2024-11")]:
        values = [m for d,m in nights.items() if first <= d < last]
        windows[first+"/"+last] = [len(values), round(statistics.median(values),3)] if values else None
    ir = {}
    seasonal_ir = {}
    for objid in sorted({r[0] for r in data["wise_epochs"]}):
        for band in ("w1","w2"):
            rows = [r for r in data["wise_epochs"] if r[0]==objid and r[1]==band and r[5]>=5]
            # Compare sampled visit medians without mixing forced fluxes across bands.
            a = [r[3] for r in rows if 57000<=r[2]<58000]
            b = [r[3] for r in rows if 58600<=r[2]<59100]
            if len(a)>=3 and len(b)>=2 and statistics.median(a)>0:
                ir[f"{objid}_{band}"] = [len(a), len(b), round(statistics.median(a),3),
                                         round(statistics.median(b),3), round(statistics.median(b)/statistics.median(a),3)]
            # The two scan directions alternate across the ~six-month visits.
            season_ratios = []
            for anchor in (57026,57198):
                season = [r for r in rows if abs((r[2]-anchor)/365.25 - round((r[2]-anchor)/365.25))<0.15]
                early = [r[3] for r in season if 57000<=r[2]<58000]
                late = [r[3] for r in season if 58600<=r[2]<59100]
                if early and late and statistics.median(early)>0:
                    season_ratios.append(round(statistics.median(late)/statistics.median(early),3))
            seasonal_ir[f"{objid}_{band}"] = season_ratios
    return dict(n_ero=data["n_ero"],n_alerts=data["n_gaia"],matches=data["matches"],
                shifts=data["control_counts"],gaia_usable=len(curve),gaia_numeric_nights=len(nights),
                gaia_nightly_windows=windows,wise_median_windows=ir, wise_same_scan_ratios=seasonal_ir)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    args=parser.parse_args()
    data=live() if args.live else json.loads((HERE/"evidence.json").read_text())
    print(json.dumps(summarize(data),indent=2))


if __name__ == "__main__":
    main()

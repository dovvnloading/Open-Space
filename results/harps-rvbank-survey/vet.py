"""External signal checks and approximate joint single-companion mass fits."""
import argparse
import csv
from datetime import datetime, timezone
from collections import defaultdict
import hashlib
import io
import os
from pathlib import Path
import re
import time
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import requests
from astropy.constants import G, M_sun, M_jup, au
from astropy.io import fits
from astropy.table import Table
from scipy.linalg import solve_triangular
from scipy.optimize import brentq, least_squares
from scipy.stats import chi2

NASA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
EU = "https://voparis-tap-planeto.obspm.fr/tap/sync"
VIZIER = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
HGCA = "https://cdsarc.cds.unistra.fr/ftp/J/ApJS/254/42/HGCA_vEDR3.fits"
TWOPI = 2 * np.pi


def chunks(values, size):
    values = sorted(set(values))
    for start in range(0, len(values), size):
        yield values[start:start + size]


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def norm(value):
    return re.sub(r"[\s_*]", "", str(value)).casefold()


def names(star):
    raw = str(star).strip().replace("_", " ")
    spaced = re.sub(r"^(HD|HIP|GJ|HR|MCC|LHS|LTT|LSPM|TYC|2MASS|CD|BD)(?=[\d+-])",
                    r"\1 ", raw, flags=re.I)
    return {raw, spaced}


def sql_list(values):
    return ",".join("'" + str(v).replace("'", "''") + "'" for v in values)


def possep(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = np.radians([ra, dec, ra0, dec0])
    a = np.sin((dec - dec0) / 2)**2 + np.cos(dec) * np.cos(dec0) * np.sin((ra - ra0) / 2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))) * 3600


class Sources:
    def __init__(self, work, out):
        self.audit = out/"catalog-inputs.csv"
        self.fields = ["source","endpoint","query","status","format","rows","bytes","sha256",
                       "retrieved_utc","accessed_utc","cache_key","detail"]
        with self.audit.open("w",newline="") as handle:
            csv.DictWriter(handle,fieldnames=self.fields).writeheader()
        self.cache = work / "catalog-cache"
        self.cache.mkdir(parents=True, exist_ok=True)
        self.failures = defaultdict(int)
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "Open-Space-HARPS-survey/1.0 (research catalog queries)"

    def get(self, label, url, query="", fmt="votable"):
        params = {"QUERY": query, "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt,
                  "MAXREC": 1000000} if query else {}
        key = hashlib.sha256((url + repr(params)).encode()).hexdigest()
        path = self.cache / key
        entry = dict(source=label, endpoint=url, query=query, format=fmt, cache_key=key)
        entry["accessed_utc"] = datetime.now(timezone.utc).isoformat()
        def record():
            with self.audit.open("a",newline="") as handle:
                csv.DictWriter(handle,fieldnames=self.fields).writerow(entry)
        if self.failures[label] >= 3 and not path.exists():
            entry.update(status="unavailable after three consecutive failed requests")
            record()
            return None
        for attempt in range(2):
            try:
                if path.exists():
                    data = path.read_bytes()
                else:
                    if url == NASA or not query:
                        response = self.session.get(url, params=params, timeout=(12, 120))
                    else:
                        response = self.session.post(url, data=params, timeout=(12, 120))
                    response.raise_for_status()
                    data = response.content
                entry.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
                if fmt == "fits":
                    with fits.open(io.BytesIO(data), memmap=False) as hdus:
                        result = hdus[1].data.copy()
                elif fmt == "csv":
                    result = pd.read_csv(io.BytesIO(data), dtype=str, keep_default_na=False)
                    if len(result.columns) < 2 or len(result) >= params["MAXREC"]:
                        raise ValueError("unexpected or truncated TAP CSV")
                else:
                    root = ET.fromstring(data)
                    for item in root.iter():
                        if item.get("name") == "QUERY_STATUS" and item.get("value") != "OK":
                            raise ValueError(str(item.get("value")) + " " + str(item.text)[:250])
                    tab = Table.read(io.BytesIO(data), format="votable", use_names_over_ids=True)
                    # Convert identifiers before pandas can upcast masked INT64 to float64.
                    result = pd.DataFrame({k: ["" if np.ma.is_masked(x) else
                                              (x.decode() if isinstance(x, bytes) else str(x))
                                              for x in tab[k]] for k in tab.colnames})
                    result.columns = [c.lower() for c in result.columns]
                if not path.exists():
                    path.write_bytes(data)
                self.failures[label] = 0
                entry.update(status="ok", rows=len(result),
                             retrieved_utc=datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat())
                record()
                return result
            except (ValueError, OSError, ET.ParseError, requests.RequestException, IndexError) as exc:
                entry.update(status="failed", detail=str(exc)[:500])
                if attempt == 0:
                    time.sleep(1)
        self.failures[label] += 1
        record()
        return None


def batch_catalog(db, label, url, template, ids, key, size=100):
    found = defaultdict(list)
    checked = set()
    for group in chunks(ids, size):
        table = db.get(label, url, template.format(ids=sql_list(group) if key == "id" else ",".join(group)))
        if table is None or key not in table:
            continue
        checked.update(group)
        for row in table.to_dict("records"):
            found[row[key]].append(row)
    return found, checked


def identify(db, targets):
    choices = {star: names(star) for star in targets.index}
    all_ids = set().union(*choices.values()) if choices else set()
    table, _ = batch_catalog(db, "SIMBAD identity", SIMBAD,
        "SELECT i.id, b.oid, b.main_id,b.ra,b.dec,b.otype FROM ident AS i "
        "JOIN basic AS b ON b.oid=i.oidref WHERE i.id IN ({ids})", all_ids, "id")
    by_name = defaultdict(list)
    for key, rows in table.items():
        by_name[norm(key)].extend(rows)
    best = {}
    for star, alternatives in choices.items():
        hits = {r["oid"]: r for n in alternatives for r in by_name[norm(n)]}
        if len(hits) == 1:
            row = next(iter(hits.values()))
            # Both positions are nominally J2000. Reject gross disagreements.
            t = targets.loc[star]
            sep = possep(num(row["ra"]), num(row["dec"]), num(t.ra_deg), num(t.dec_deg))
            if np.isnan(sep) or sep < 30:
                best[star] = row
    # Coordinate rescue is limited to one SIMBAD source in 3 arcsec.
    missing = targets.loc[~targets.index.isin(best)].copy()
    missing = missing[np.isfinite(pd.to_numeric(missing.ra_deg, errors="coerce")) &
                      np.isfinite(pd.to_numeric(missing.dec_deg, errors="coerce"))]
    for group in chunks(missing.index, 15):
        cones = [f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{float(missing.loc[s].ra_deg):.8f},"
                 f"{float(missing.loc[s].dec_deg):.8f},0.000834))" for s in group]
        data = db.get("SIMBAD position", SIMBAD,
                      "SELECT oid,main_id,ra,dec,otype FROM basic WHERE " + " OR ".join(cones))
        if data is None:
            continue
        rows = data.to_dict("records")
        for star in group:
            t = missing.loc[star]
            hits = [r for r in rows if possep(num(r["ra"]), num(r["dec"]), num(t.ra_deg), num(t.dec_deg)) <= 3]
            if len(hits) == 1:
                best[star] = hits[0]
    return best


def pm_data(row):
    cov, pm = {}, {}
    for epoch in ["hip", "hg", "gaia"]:
        err = np.array([row[f"pmra_{epoch}_error"], row[f"pmdec_{epoch}_error"]])
        r = row[f"pmra_pmdec_{epoch}"]
        if not np.all(err > 0) or not abs(r) < 1:
            raise ValueError("invalid HGCA covariance")
        cov[epoch] = np.outer(err, err) * np.array([[1., r], [r, 1.]])
        pm[epoch] = np.array([row[f"pmra_{epoch}"], row[f"pmdec_{epoch}"]])
    delta = np.r_[pm["hip"] - pm["hg"], pm["gaia"] - pm["hg"]]
    covariance = np.block([[cov["hip"] + cov["hg"], cov["hg"]],
                           [cov["hg"], cov["gaia"] + cov["hg"]]])
    return delta, np.linalg.cholesky(covariance)


def mass(period, k, ecc, primary, sini=1.):
    if not (period > 0 and k > 0 and 0 <= ecc < 1 and primary > 0 and 0 < sini <= 1):
        return np.nan
    function = k**3 * (period * 86400) * (1 - ecc**2)**1.5 / (TWOPI * G.value * M_sun.value)
    def equation(m2):
        return (m2 * sini)**3 / (primary + m2)**2 - function
    upper = max(1., 3 * function / sini**3 + 3 * primary)
    return brentq(equation, 0., upper, xtol=1e-14) * (M_sun.value / M_jup.value)


def orbit(t, p, ecc, omega, phase, tref):
    mean = (phase + TWOPI * (t - tref) / p + np.pi) % TWOPI - np.pi
    anomaly = mean + 0.85 * ecc * np.sign(np.sin(mean))
    for _ in range(10):
        step = (anomaly - ecc * np.sin(anomaly) - mean) / (1 - ecc * np.cos(anomaly))
        anomaly -= step
        if np.max(np.abs(step)) < 1.e-10:
            break
    x, y = np.cos(anomaly) - ecc, np.sqrt(1 - ecc**2) * np.sin(anomaly)
    u, v = x * np.cos(omega) - y * np.sin(omega), x * np.sin(omega) + y * np.cos(omega)
    rv = u / (1 - ecc * np.cos(anomaly)) + ecc * np.cos(omega)
    return u, v, rv


def joint_fit(row, sig, epochs, primary, mass_err):
    """Fit RV and four correlated proper-motion differences, one dark companion."""
    from search import nightly
    data = nightly(epochs)
    times = data.bjd.to_numpy(dtype=float)
    y = data.rv.to_numpy(dtype=float)
    err = data.err.to_numpy(dtype=float)
    group = pd.factorize(data.era)[0]
    p0, k0, e0 = map(float, (sig.period_d, sig.k_ms, sig.ecc))
    if np.ptp(times) < 1.25 * p0 or p0 < 4 * 365.25 or k0 <= 0 or e0 >= .9:
        return {"joint_status":"insufficient long-orbit coverage"}
    tref = np.mean(times)
    omega = num(sig.omega_rad)
    phase = (TWOPI * (tref - num(sig.tperi_bjd)) / p0) % TWOPI
    centers = np.array([[row["epoch_ra_" + which], row["epoch_dec_" + which]]
                        for which in ("hip", "gaia")])
    # Fit uniform-window positions and slopes at each coordinate's central epoch.
    dt = np.linspace(-.5, .5, 25) * np.array([3.25, 34 / 12])[:, None, None]
    years = centers[:, :, None] + dt
    obs = (years - 2000.) * 365.25 + 2451545.
    delta, chol = pm_data(row)
    plx, plxerr = float(row["parallax_gaia"]), float(row["parallax_gaia_error"])
    if not (plx > 10 * plxerr > 0):
        return {"joint_status":"parallax too uncertain"}
    masks = [group == g for g in np.unique(group)]
    _, _, template = orbit(times, p0, e0, omega, phase, tref)
    residual = y - k0 * template
    sigmas = err.copy()
    for m in masks:
        mad = 1.4826 * np.median(abs(residual[m] - np.median(residual[m])))
        extra = max(0., mad**2 - np.median(err[m])**2)
        sigmas[m] = np.sqrt(err[m]**2 + extra)
    def errors(q):
        logp, logk, ecc, arg, m0, inc, node, parallax = q
        period, k = np.exp(logp), np.exp(logk)
        _, _, model = orbit(times, period, ecc, arg, m0, tref)
        residual = y - k * model
        for m in masks:
            residual[m] -= np.average(residual[m], weights=sigmas[m]**-2)
        u, v, _ = orbit(obs, period, ecc, arg, m0, tref)
        a = k * period * 86400 * np.sqrt(1 - ecc**2) / (TWOPI * np.sin(inc) * au.value) * parallax
        east = np.cos(node) * u - np.sin(node) * np.cos(inc) * v
        north = np.sin(node) * u + np.cos(node) * np.cos(inc) * v
        positions = a * np.stack([east[:, 0, :], north[:, 1, :]], axis=1)
        means = positions.mean(axis=-1)
        motions = (positions * dt).sum(axis=-1) / (dt**2).sum(axis=-1)
        hg = (means[1] - means[0]) / (centers[1] - centers[0])
        whitened = solve_triangular(chol, (motions - hg).ravel() - delta, lower=True)
        return np.r_[residual / sigmas, whitened, (parallax - plx) / plxerr]
    lower = [np.log(.7 * p0), np.log(.2 * k0), 0., -np.inf, -np.inf, .002, -np.inf, .01*plx]
    upper = [np.log(1.3 * p0), np.log(5 * k0), .95, np.inf, np.inf, np.pi-.002, np.inf, 2*plx]
    solutions = []
    for inc in [.45, 1.1, 2., 2.7]:
        for node in [0., np.pi]:
            start = [np.log(p0), np.log(k0), np.clip(e0, .001, .89), omega, phase, inc, node, plx]
            fit = least_squares(errors, start, bounds=(lower, upper), x_scale="jac", max_nfev=400)
            if fit.success:
                solutions.append(fit)
    if not solutions:
        return {"joint_status":"optimizer failed"}
    fit = min(solutions, key=lambda f:f.cost)
    q = fit.x
    m2 = mass(np.exp(q[0]), np.exp(q[1]), q[2], primary, np.sin(q[5]))
    chi_sq = np.dot(fit.fun[-5:-1], fit.fun[-5:-1])
    result = dict(joint_status="conditional mass only", joint_trial_mass_mjup=m2,
                  joint_inc_deg=np.degrees(q[5]), joint_period_d=np.exp(q[0]),
                  joint_pma_chi2=chi_sq, joint_rv_chi2=float(fit.fun[:-5] @ fit.fun[:-5]))
    if chi_sq > 13.3 or np.any(fit.active_mask) or sum(fit.fun**2) > 2 * max(1,len(fit.fun)-len(q)-len(masks)):
        return dict(result, joint_status="poor fit or boundary")
    if not np.isfinite(mass_err) or mass_err <= 0 or mass_err >= primary:
        return dict(result, joint_status="host-mass uncertainty unavailable")
    # Profile in physical mass. At edge-on orientation the first derivative of
    # mass with respect to inclination vanishes even for unconstrained inclinations.
    lower_profile = [lower[j] for j in (0,2,3,4,5,6,7)] + [max(.01,primary-5*mass_err)]
    upper_profile = [upper[j] for j in (0,2,3,4,5,6,7)] + [primary+5*mass_err]
    seed = np.r_[q[[0,2,3,4,5,6,7]],primary]
    masses = np.unique(np.r_[np.geomspace(max(.00001, m2/4),max(160.,4*m2),55), m2])
    costs = []
    for trial in masses:
        solar = trial * M_jup.value/M_sun.value
        def profile_residual(par):
            logp, ecc, arg, m0, inc, node, parallax, host = par
            function = (solar*np.sin(inc))**3/(host+solar)**2
            k = (TWOPI*G.value*M_sun.value*function/(86400*np.exp(logp)))**(1/3)/np.sqrt(1-ecc**2)
            return np.r_[errors([logp,np.log(k),ecc,arg,m0,inc,node,parallax]),(host-primary)/mass_err]
        sini = (np.sin(q[5])*m2/trial *
                ((primary+solar)/(primary + m2*M_jup.value/M_sun.value))**(2/3))
        angle = np.arcsin(np.clip(sini,.002,.999999))
        trial_costs = []
        for inc in (angle, np.pi-angle):
            start = seed.copy()
            start[4] = inc
            prof = least_squares(profile_residual,start,bounds=(lower_profile,upper_profile),
                                 x_scale="jac",max_nfev=250)
            if prof.success:
                trial_costs.append(2*prof.cost)
        if not trial_costs:
            return dict(result, joint_status="mass profile optimization failed")
        costs.append(min(trial_costs))
    inside = np.where(np.array(costs)<=min(2*fit.cost,min(costs))+3.84)[0]
    if len(inside)==0 or inside[0]==0 or inside[-1]==len(masses)-1:
        return dict(result, joint_status="unbounded mass profile")
    # Include one exterior grid step on each side (conservative rounding).
    lo, hi = masses[inside[0]-1], masses[inside[-1]+1]
    result.update(profile_mass_lo_mjup=lo,profile_mass_hi_mjup=hi)
    if (lo > m2 or hi < m2 or hi/lo>3 or hi>=80 or
            min(costs) < 2*fit.cost-2):
        return dict(result, joint_status="multimodal, broad or stellar-mass profile")
    return dict(result, joint_status="approximate dynamical mass", dynamical_mass_mjup=m2)



class HostIndex:
    """Name, Gaia source ID and conservative position index for companion catalogs."""
    def __init__(self, table, identifiers):
        from scipy.spatial import cKDTree
        self.names = defaultdict(set)
        self.rows = [] if table is None else table.to_dict("records")
        xyz, self.sky_rows = [], []
        self.tree = None
        for j, row in enumerate(self.rows):
            for field in identifiers:
                if row.get(field):
                    self.names[norm(row[field])].add(j)
            ra, dec = num(row.get("ra")), num(row.get("dec"))
            if np.isfinite(ra + dec):
                xyz.append(self.xyz(ra, dec))
                self.sky_rows.append(j)
        if xyz:
            self.tree = cKDTree(xyz)

    @staticmethod
    def xyz(ra, dec):
        a, d = np.radians([ra,dec])
        return np.array([np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)])

    def match(self, aliases, ra, dec):
        selected = set()
        for a in aliases:
            selected.update(self.names[norm(a)])
        if self.tree is not None and np.isfinite(ra + dec):
            # Nearby components are conservative possible conflicts; never use proximity alone for mass.
            hit = self.tree.query_ball_point(self.xyz(ra, dec), np.radians(4/3600))
            selected.update(self.sky_rows[i] for i in hit)
        return [self.rows[j] for j in sorted(selected)]


def match_period(period, reference, baseline):
    if not (period > 0 and reference > 0 and baseline > 0):
        return ""
    f = 1 / reference
    # Direct matches take priority. Limit loose one-Rayleigh windows at long periods.
    for kind, test in [("period", f), ("harmonic", 2*f), ("harmonic", f/2),
                       ("daily alias", abs(f-1)), ("daily alias", f+1),
                       ("sidereal-day alias", abs(f-1.002738)), ("sidereal-day alias", f+1.002738),
                       ("annual alias", abs(f-1/365.25)), ("annual alias", f+1/365.25)]:
        if test > 0 and abs(1/period-test) <= max(.012*test, min(1/baseline,.05*test)):
            return kind
    return ""


def save(out, name, rows, fields):
    frame = pd.DataFrame(rows)
    extras = [c for c in frame if c not in fields]
    frame.reindex(columns=fields+extras).to_csv(out / name, index=False)


def run(work, out):
    if os.environ.get("GITHUB_EVENT_NAME") != "workflow_dispatch":
        raise SystemExit("Archive vetting requires a workflow_dispatch job")
    db = Sources(work, out)
    targets = pd.read_csv(out/"targets.csv", dtype={"star": str}).set_index("star")
    stars = pd.read_csv(out/"stars.csv", dtype={"star": str}).set_index("star")
    signals = pd.read_csv(out/"signals.csv", dtype={"star": str})
    trend_sig = (pd.to_numeric(stars.slope_ms_yr,errors="coerce").abs() /
                 pd.to_numeric(stars.slope_err_ms_yr,errors="coerce").where(lambda x:x>0))
    trend_stars = set(stars.index[trend_sig >= 3])
    sample = set(signals.star) | trend_stars
    selected = targets.loc[targets.index.isin(sample)]
    identity = identify(db, selected)
    oids = {r["oid"] for r in identity.values()}
    ids, id_checked = batch_catalog(db, "SIMBAD aliases", SIMBAD,
        "SELECT oidref,id FROM ident WHERE oidref IN ({ids})", oids, "oidref")
    types, type_checked = batch_catalog(db, "SIMBAD types", SIMBAD,
        "SELECT oidref,otype FROM otypes WHERE oidref IN ({ids})", oids, "oidref")

    nasa = db.get("NASA pscomppars", NASA, "SELECT hostname,hd_name,hip_name,gaia_dr3_id,pl_name,"
                  "pl_orbper,st_rotp,st_mass,st_masserr1,st_masserr2,sy_snum,ra,dec FROM pscomppars", "csv")
    all_nasa = db.get("NASA all orbital solutions", NASA,
                      "SELECT DISTINCT hostname,pl_name,pl_orbper FROM ps WHERE pl_orbper IS NOT NULL", "csv")
    eu = db.get("exoplanet.eu", EU, "SELECT star_name,target_name,period,publication_status,ra,dec "
                "FROM exoplanet.epn_core")
    spores = db.get("SPORES-HWO II Table 6", VIZIER, 'SELECT "HD","Per","Class","BibCode" '
                    'FROM "J/AJ/170/343/table6"')
    trend_catalog = db.get("SPORES-HWO II Table 7", VIZIER,
                           'SELECT "HD" FROM "J/AJ/170/343/table7"')
    reported_trends = {"hd25673":"previously reported by Open-Space"}
    if trend_catalog is not None and "hd" in trend_catalog:
        for hd in trend_catalog.hd:
            reported_trends[norm("HD "+re.sub(r"^HD\s*","",str(hd),flags=re.I))] = "SPORES-HWO II Table 7"
    if spores is not None and "hd" in spores:
        spores["hd"] = ["HD " + re.sub(r"^HD\s*", "", str(v).strip(), flags=re.I) for v in spores.hd]
    references = defaultdict(list)
    def add_prior(star, period, label):
        references[norm(star)].append((period,label))
    for name, periods in {"HD 32564":[2.56,11.06,23.61,49.2], "HD 71334":[89.],
                          "MCC 150":[622.], "HD 207129":[1398.], "HD 125881":[41.]}.items():
        for p in periods:
            add_prior(name, p, "prior public HARPS search (rpinho/planet-hunt)")
    if spores is not None:
        for r in spores.to_dict("records"):
            add_prior(r["hd"], num(r["per"]), f"SPORES-HWO II {r['class']} {r['bibcode']}")
    primary_catalog = HostIndex(nasa, ["hostname","hd_name","hip_name","gaia_dr3_id"])
    eu_catalog = HostIndex(eu, ["star_name"])
    solutions = HostIndex(all_nasa, ["hostname"])
    aliases = {}
    gaia_ids = {}
    hip_ids = {}
    for star in sorted(sample):
        r = identity.get(star, {})
        choices = set(names(star))
        if r:
            choices.add(r["main_id"])
            choices.update(x["id"] for x in ids[r["oid"]])
        aliases[star] = choices
        gids = {m[1] for name in choices if (m := re.fullmatch(r"Gaia DR3\s+(\d+)", name, re.I))}
        hips = {m[1] for name in choices if (m := re.fullmatch(r"HIP\s*(\d+)", name, re.I))}
        if len(gids) == 1:
            gaia_ids[star] = gids.pop()
        if len(hips) == 1:
            hip_ids[star] = hips.pop()
    gids = set(gaia_ids.values())
    gaia, gaia_checked = batch_catalog(db, "Gaia binarity", GAIA,
        "SELECT source_id,non_single_star,ruwe FROM gaiadr3.gaia_source WHERE source_id IN ({ids})", gids,"source_id")
    rotation, _ = batch_catalog(db, "Gaia rotation", GAIA,
        "SELECT source_id,best_rotation_period FROM gaiadr3.vari_rotation_modulation WHERE source_id IN ({ids})",
        gids,"source_id")
    masses, _ = batch_catalog(db, "Gaia FLAME", GAIA,
        "SELECT source_id,mass_flame,mass_flame_lower,mass_flame_upper,flags_flame "
        "FROM gaiadr3.astrophysical_parameters WHERE source_id IN ({ids})", gids,"source_id")
    # Keep all associated archival bibliography matches for the stated title search.
    words = ("planet","brown dwarf","companio","radial veloc","spectroscopic binar",
             "astrometric orbit","rotational modulat","stellar rotation","doppler","harps")
    predicate = " OR ".join("LOWER(r.title) LIKE '%"+w+"%'" for w in words)
    papers, paper_checked = batch_catalog(db, "SIMBAD bibliography", SIMBAD,
        "SELECT h.oidref,r.bibcode,r.title,r.doi FROM has_ref AS h JOIN ref AS r ON r.oidbib=h.oidbibref "
        "WHERE h.oidref IN ({ids}) AND (" + predicate + ")", oids, "oidref", size=32)
    hgca = db.get("HGCA EDR3 (Brandt 2021)", HGCA, fmt="fits")
    by_gaia, by_hip = {}, {}
    if hgca is not None:
        by_gaia = {str(k): i for i,k in enumerate(hgca["gaia_source_id"])}
        by_hip = {str(k): i for i,k in enumerate(hgca["hip_id"])}
    acc_rows, lit_rows, vet = [], [], []
    epochs = None
    by_star = {k:g for k,g in signals.groupby("star", sort=False)}
    for star in sorted(sample):
        r = identity.get(star,{})
        oid = r.get("oid", "")
        choices = aliases[star]
        keyset = {norm(a) for a in choices}
        gid, hip = gaia_ids.get(star,""), hip_ids.get(star,"")
        if star in selected.index:
            t = selected.loc[star]
            ra, dec = num(t.ra_deg), num(t.dec_deg)
        else:
            ra, dec = np.nan, np.nan
        ns = primary_catalog.match(choices, ra, dec)
        eus = eu_catalog.match(choices, ra, dec)
        hostnames = choices | {x["hostname"] for x in ns}
        sols = solutions.match(hostnames, np.nan, np.nan)
        matched = [x for x in ns if any(norm(x.get(col, "")) in keyset
                   for col in ["hostname","hd_name","hip_name","gaia_dr3_id"])]
        m, me, msrc = np.nan, np.nan, ""
        rotations = [(num(item["st_rotp"]), "NASA") for item in matched if num(item["st_rotp"])>0]
        for item in rotation[gid]:
            if num(item["best_rotation_period"])>0:
                rotations.append((num(item["best_rotation_period"]),"Gaia DR3"))
        for item in matched:
            if num(item["st_mass"]) > 0:
                m, me, msrc = num(item["st_mass"]), np.nanmax([abs(num(item["st_masserr1"])),
                             abs(num(item["st_masserr2"]))]), "NASA pscomppars"
                break
        if not msrc and len(masses[gid]) == 1:
            item = masses[gid][0]
            if item["flags_flame"].strip() == "00" and num(item["mass_flame"]) > 0:
                m = num(item["mass_flame"])
                me = max(m-num(item["mass_flame_lower"]), num(item["mass_flame_upper"])-m)
                msrc = "Gaia DR3 FLAME (00)"
        binary_types = {x["otype"] for x in types[oid]} & {
            "**","**?","SB*","SB?","EB*","EB?","El*","Al*","bL*","WU*","Sy*","CV*"}
        binary = {("SIMBAD "+x) for x in binary_types}
        if any(num(x["sy_snum"])>1 for x in ns):
            binary.add("NASA multiple")
        if any(num(x["non_single_star"])>0 for x in gaia[gid]):
            binary.add("Gaia NSS")
        ruwe = max((num(x["ruwe"]) for x in gaia[gid]), default=np.nan)
        if ruwe > 1.4:
            binary.add("elevated Gaia RUWE")
        binary_searched = oid in type_checked and gid in gaia_checked and len(gaia[gid])==1 and bool(oid)
        bstatus = " | ".join(sorted(binary)) if binary else ("no catalog flag" if binary_searched else "incomplete")
        ids_found = {by_gaia[gid]} if gid in by_gaia else set()
        if hip in by_hip:
            ids_found.add(by_hip[hip])
        # Reject conflicting HIP or Gaia IDs, including HIP match with disagreeing DR3 ID.
        h = hgca[next(iter(ids_found))] if len(ids_found)==1 else None
        if h is not None and ((gid and str(h["gaia_source_id"])!=gid) or (hip and str(h["hip_id"])!=hip)):
            h = None
        pma_prob = np.nan
        slope = num(stars.loc[star,"slope_ms_yr"])
        se = num(stars.loc[star,"slope_err_ms_yr"])
        a = dict(star=star,slope_ms_yr=slope,slope_err_ms_yr=se,gaia_dr3_id=gid,hip=hip,
                 astrometry_check="no exact match" if hgca is not None else "catalog unavailable",
                 trend_prior=" | ".join(sorted({reported_trends[n] for n in keyset if n in reported_trends})))
        if h is not None:
            try:
                delta, chol = pm_data(h)
                z = solve_triangular(chol, delta, lower=True)
                pma_prob = chi2.sf(z @ z,4)
                a.update(astrometry_check="matched", pma_chi2_4dof=z@z,pma_p=pma_prob,
                         hip_minus_hg_ra_masyr=delta[0],hip_minus_hg_dec_masyr=delta[1],
                         gaia_minus_hg_ra_masyr=delta[2],gaia_minus_hg_dec_masyr=delta[3])
            except (ValueError, np.linalg.LinAlgError):
                a["astrometry_check"]="unusable covariance"
        if star in trend_stars:
            acc_rows.append(a)
        for paper in papers[oid]:
            lit_rows.append(dict(star=star,**{k:paper[k] for k in ["bibcode","title","doi"]}))
        if star not in by_star:
            continue
        g = by_star[star]
        for sig in g.itertuples(index=False):
            p = num(sig.period_d)
            baseline = num(stars.loc[star,"baseline_d"])
            prior, alias_matches, rot_hits = set(), set(), set()
            reported = [(num(x["pl_orbper"]), x["pl_name"]) for x in ns+sols]
            reported += [(num(x["period"]), "exoplanet.eu "+x["target_name"]) for x in eus]
            reported += [v for n in keyset for v in references[n]]
            for published, label in reported:
                kind = match_period(p, published, baseline)
                if kind=="period":
                    prior.add(f"{label} ({published:g} d)")
                elif kind:
                    alias_matches.add(f"{label} ({kind} of {published:g} d)")
            for prot, label in rotations:
                if any(abs(p/(ratio*prot)-1)<.12 for ratio in [.333333,.5,1,2]):
                    rot_hits.add(f"{label} {prot:g} d")
            rotation_check = " | ".join(sorted(rot_hits)) if rot_hits else (
                "clear of catalog periods" if rotations else "unavailable")
            minimum = mass(p,num(sig.k_ms),num(sig.ecc),m)
            reasons = []
            if not r or oid not in id_checked:
                reasons.append("unresolved identity or aliases")
            if nasa is None or all_nasa is None or eu is None or spores is None:
                reasons.append("prior catalog incomplete")
            if alias_matches:
                reasons.append("possible published harmonic or alias")
            if ns and not matched:
                reasons.append("nearby catalogued host needs identification")
            if any(label for pp,label in reported if not np.isfinite(pp)):
                reasons.append("known companion without comparable period")
            if bstatus!="no catalog flag":
                reasons.append("multiplicity unresolved")
            if not rotations:
                reasons.append("external rotation period unavailable")
            if oid not in paper_checked:
                reasons.append("bibliography unavailable")
            if not np.isfinite(minimum):
                reasons.append("host mass unavailable")
            if minimum>=80:
                reasons.append("minimum mass in stellar range")
            if rot_hits:
                reasons.append("rotation or harmonic overlap")
            if "hd25673" in keyset:
                reasons.append("published drift and existing repo study")
            if sig.rv_status != "pass":
                reasons.append("RV "+str(sig.rv_status))
            status = ("prior signal" if prior else "rejected" if
                      sig.rv_status=="rejected" or rot_hits or minimum>=80 else
                      "inconclusive" if reasons else "requires literature review")
            outrow = dict(star=star,signal_index=sig.signal_index,period_d=p,k_ms=sig.k_ms,
                          ecc=sig.ecc,fap=sig.fap,rv_status=sig.rv_status,identity=r.get("main_id",""),
                          aliases=" | ".join(sorted(choices)),gaia_dr3_id=gid,external_status=status,
                          reasons=" | ".join(reasons),rotation_check=rotation_check,
                          catalog_rotation_periods=" | ".join(f"{v:g} {s}" for v,s in sorted(set(rotations))),
                          binary_check=bstatus,ruwe=ruwe,host_mass_msun=m,host_mass_err_msun=me,
                          host_mass_source=msrc,minimum_mass_mjup=minimum,dynamical_mass_mjup=np.nan,
                          prior_signal=" | ".join(sorted(prior)),alias_matches=" | ".join(sorted(alias_matches)),
                          nasa_check=("source unavailable" if nasa is None or all_nasa is None else
                                      "host associated" if ns or sols else "no identified host"),
                          literature_check=("indexed titles retrieved" if oid in paper_checked else "unavailable"),
                          bibliography_matches=len(papers[oid]),joint_status="no qualifying single long orbit",
                          pma_p=pma_prob)
            others = [pp for pp,_ in reported if np.isfinite(pp) and match_period(p, pp,baseline)!="period"]
            other_rv = g[(g.signal_index!=sig.signal_index) & (g.fap<=.01)]
            if (h is not None and pma_prob<.001 and np.isfinite(minimum) and minimum<80 and
                    p>=4*365.25 and baseline>=1.25*p and num(sig.fap)<=.01 and num(sig.delta_bic)>=10 and
                    sig.rv_status!="rejected" and not others and not (binary-{"elevated Gaia RUWE"}) and binary_searched and other_rv.empty):
                if epochs is None:
                    epochs = pd.read_csv(work/"epochs.csv.gz",dtype={"star":str}).groupby("star",sort=False)
                try:
                    outrow.update(joint_fit(h,sig,epochs.get_group(star),m,me))
                except (ValueError,np.linalg.LinAlgError,KeyError,FloatingPointError) as exc:
                    outrow["joint_status"]="unusable fit: "+str(exc)[:180]
            vet.append(outrow)
    save(out, "vetting.csv", vet, ["star","signal_index","period_d","k_ms","ecc","fap","rv_status",
         "external_status","minimum_mass_mjup","dynamical_mass_mjup","reasons","identity",
         "rotation_check","binary_check","nasa_check","literature_check","joint_status"])
    save(out,"accelerations.csv",acc_rows,["star","slope_ms_yr","slope_err_ms_yr","astrometry_check",
         "trend_prior","pma_p","pma_chi2_4dof"])
    save(out,"literature.csv.gz",lit_rows,["star","bibcode","title","doi"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    run(args.work, args.out)

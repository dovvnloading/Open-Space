"""Check the shape predicted by Gaia at independent LAMOST exposure times."""
import argparse
import ast
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_MJD = 57388.5  # J2016.0, JD 2457389.0
AU_KM = 149597870.7
FLOOR = 1.0  # km/s, common nightly term


def write_csv(path, rows):
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def unpack(row, count, rng):
    kind = row["nss_solution_type"]
    if kind in ("SB1", "SB1C"):
        cols = ["period", "center_of_mass_velocity", "semi_amplitude_primary"]
        if kind == "SB1":
            cols += ["eccentricity", "arg_periastron"]
        cols += ["t_periastron"]
    elif kind == "AstroSpectroSB1" or kind.startswith("Orbital"):
        cols = ["ra", "dec", "parallax", "pmra", "pmdec"] + [c + "_thiele_innes" for c in "abfg"]
        if kind == "AstroSpectroSB1":
            cols += ["c_thiele_innes", "h_thiele_innes", "center_of_mass_velocity"]
        if kind in ("Orbital", "AstroSpectroSB1"):
            cols += ["eccentricity", "period"]
        else:
            cols += ["period", "eccentricity"]
        cols += ["t_periastron"]
    else:
        raise ValueError(kind)
    vals = np.array([float(row[k]) for k in cols])
    errors = np.array([float(row[k + "_error"]) for k in cols])
    n = len(cols)
    corr = np.eye(n)
    corr[np.tril_indices(n, -1)] = ast.literal_eval(row["corr_vec"])
    corr += np.tril(corr, -1).T
    cov = corr * np.outer(errors, errors)
    d = rng.multivariate_normal(vals, cov, count, check_valid="raise")
    return {k: np.r_[v, d[:, j]] for j, (k, v) in enumerate(zip(cols, vals))}


def predictions(row, times, rng):
    d = unpack(row, 4096, rng)
    p, tp = d["period"], d["t_periastron"]
    e = d.get("eccentricity", np.zeros_like(p))
    kind = row["nss_solution_type"]
    if kind in ("SB1", "SB1C"):
        k = d["semi_amplitude_primary"]
        omega = np.radians(d.get("arg_periastron", np.zeros_like(k)))
    elif kind == "AstroSpectroSB1":
        c, h = d["c_thiele_innes"], d["h_thiele_innes"]
        omega = np.arctan2(c, h)
        with np.errstate(invalid="ignore"):
            k = 2*np.pi*AU_KM*np.hypot(c,h)/(86400*p*np.sqrt(1-e**2))
    else:
        a, b, f, g = [d[i + "_thiele_innes"] for i in "abfg"]
        omega = (np.arctan2(b-f, a+g) - np.arctan2(b+f, a-g))/2
        # Projected size = alpha sin(i) = sqrt(u*v).
        u, v = np.hypot(a+g, b-f), np.hypot(a-g, b+f)
        with np.errstate(invalid="ignore"):
            k = 2*np.pi*AU_KM*np.sqrt(u*v)/(86400*p*d["parallax"]*np.sqrt(1-e**2))
    ok = (p > 0) & (e >= 0) & (e < 1) & (k > 0)
    if not ok[0]:
        raise ValueError("Invalid nominal orbit")
    p, e, tp, k, omega = [x[ok, None] for x in (p, e, tp, k, omega)]
    m = (2*np.pi*(times-REF_MJD-tp)/p + np.pi) % (2*np.pi) - np.pi
    # Bisection stays on the physical solution near e=1.
    lo, hi = m-e, m+e
    for _ in range(38):
        E = (lo+hi)/2
        above = E-e*np.sin(E) > m
        hi, lo = np.where(above,E,hi), np.where(above,lo,E)
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return k*(np.cos(nu+omega)+e*np.cos(omega))


def verified_rows():
    exposures = defaultdict(dict)
    matched = set()
    total = 0
    for receipt in json.loads((HERE/"raw/receipts.json").read_text()):
        data = gzip.decompress((HERE/"raw"/receipt["file"]).read_bytes())
        if hashlib.sha256(data).hexdigest() != receipt["sha256_csv"]:
            raise ValueError("LAMOST input hash mismatch")
        for r in csv.DictReader(data.decode("utf-8-sig").splitlines()):
            sid = r["gaia_source_id"]
            matched.add(sid)
            total += 1
            if r["coadd"]=="0" and r["lmjm"] and r["fibertype"] in ("Obj","F-std"):
                exposures[(sid, r["obsid"], r["lmjm"])][r["band"]] = r
    clean = defaultdict(list)
    for (sid, obs, minute), arms in sorted(exposures.items()):
        if set(arms) != {"B","R"}:
            continue
        b,r = arms["B"], arms["R"]
        if any(z["offsets"]!="0" or z["fibermask"]!="0" or float(z["snr"])<10
               or z["bad_b"]!="0" or z["bad_r"]!="0" for z in (b,r)):
            continue
        if any(z["rv_br_flag"]!="0" for z in (b,r)):
            continue
        y, sig = float(b["rv_br1"]), float(b["rv_br1_err"])
        if not (0 < sig < 100 and abs(y) < 450) or b['rv_br1'] != r['rv_br1']:
            continue
        clean[sid].append(dict(source_id=sid,obsid=obs,lmjm=minute,
                              mjd_utc_start=int(minute)/1440-1/3,obsdate=b['obsdate'],
                              rv_kms=y,err_kms=sig,snr_b=b['snr'],snr_r=r['snr']))
    return clean, dict(returned_rows=total, raw_source_matches=len(matched),
                       single_exposures=len(exposures), both_arm_clean_exposures=sum(map(len,clean.values())),
                       clean_sources=len(clean))


def nightly(rows, y, sigma, model):
    dates = sorted({r['obsdate'] for r in rows})
    obs, variance, models = [], [], []
    for date in dates:
        mask = [r['obsdate']==date for r in rows]
        w = sigma[mask]**-2
        obs.append(np.average(y[mask],weights=w))
        variance.append(1/w.sum()+FLOOR**2)
        models.append(np.average(model[:,mask],axis=1,weights=w))
    return dates,np.array(obs),np.array(variance),np.array(models).T


def main(nss):
    mass = list(csv.DictReader((ROOT/'gaia/mass_screen.csv').open()))
    mass = {(r['source_id'],r['nss_solution_type']):float(r['m2_screen_msun']) for r in mass}
    clean, counts = verified_rows()
    write_csv(HERE/'exposures.csv',[r for rows in clean.values() for r in rows])
    multi = {sid for sid,rows in clean.items() if len({r['obsdate'] for r in rows})>=2}
    counts['clean_multinight_sources'] = len(multi)
    counts['gaia_draws_requested'] = 4096
    totals, residuals, skipped = [], [], []
    with gzip.open(nss,'rt') as f:
        for row in csv.DictReader(line for line in f if not line.startswith('#')):
            sid, kind = row['source_id'],row['nss_solution_type']
            if sid not in multi or (sid,kind) not in mass:
                continue
            rows = clean[sid]
            times, y, sigma = [np.array([r[k] for r in rows]) for k in ('mjd_utc_start','rv_kms','err_kms')]
            try:
                model = predictions(row,times,np.random.default_rng(int(sid)))
            except (ValueError, np.linalg.LinAlgError) as exc:
                skipped.append([sid,kind,str(exc)])
                continue
            dates, v, var, m = nightly(rows,y,sigma,model)
            w = 1/var
            r = v-m
            r -= np.average(r,axis=1,weights=w)[:,None]
            x2 = np.sum(r*r*w,axis=1)
            if kind.startswith('Orbital'):
                r2 = v+m
                r2 -= np.average(r2,axis=1,weights=w)[:,None]
                x2_second = np.sum(r2*r2*w,axis=1)
                if x2_second[0]<x2[0]:
                    r[0]=r2[0]
                    m[0]*=-1
                x2=np.minimum(x2,x2_second)
            null = np.sum((v-np.average(v,weights=w))**2*w)
            dof = len(v)-1
            q = chi2.sf(x2,dof)
            totals.append(dict(source_id=sid,solution=kind,m2_screen=round(mass[sid,kind],4),n_exp=len(rows),
                n_nights=len(v),baseline_d=round(np.ptp(times),4),period_d=round(float(row['period']),6),
                measured_span_kms=round(np.ptp(v),3),nominal_model_span_kms=round(np.ptp(m[0]),3),
                chi2_constant=round(null,4),chi2_nominal=round(x2[0],4),dof=dof,
                accepted_gaia_draws=len(q)-1,q95_compatibility=round(float(np.quantile(q[1:],.95)),8),
                fraction_draws_p_gt_0p01=round(float((q[1:]>.01).mean()),6)))
            for d,obs,variance,pr,dr in zip(dates,v,var,m[0],r[0]):
                residuals.append(dict(source_id=sid,solution=kind,obsdate=d,rv_kms=round(obs,5),
                    err_kms=round(np.sqrt(variance),5),nominal_shifted_rv_kms=round(obs-dr,5),resid_kms=round(dr,5)))
    counts.update(evaluated_solutions=len(totals),skipped=skipped)
    write_csv(HERE/'shape_checks.csv',sorted(totals,key=lambda x:-x['m2_screen']))
    write_csv(HERE/'nightly_residuals.csv',residuals)
    (HERE/'depth.json').write_text(json.dumps(counts,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('nss',type=Path)
    main(p.parse_args().nss)

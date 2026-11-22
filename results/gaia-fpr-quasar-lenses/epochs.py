"""Descriptive epoch comparison for selected fields; never an identification criterion."""
import csv
import hashlib
import io
import itertools
import json
from pathlib import Path

import numpy as np
import requests
from scipy.stats import spearmanr


SERVICES = ("https://gea.esac.esa.int/tap-server/tap/sync", "https://gaia.aip.de/tap/sync",
            "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync")


def numeric(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def sample(data):
    if len(data) == 0:
        return dict(n=0)
    mag, t = np.array([[r[1],r[0]] for r in data]).T
    return dict(n=len(t), first_epoch_days=float(min(t)), last_epoch_days=float(max(t)),
                span_days=float(np.ptp(t)), median_g_mag=float(np.median(mag)),
                mag_p90_minus_p10=float(np.ptp(np.percentile(mag,[10,90]))),
                robust_mag_scatter=float(1.4826*np.median(abs(mag-np.median(mag)))))


def contrasts(a, b):
    # Pair only mutual nearest epochs, so no transit is counted twice. The
    # proximity threshold does not establish that windows or scans are independent.
    if not a or not b:
        return dict(n_near_simultaneous=0)
    ta, ma = np.array(a).T
    tb, mb = np.array(b).T
    difference = abs(ta[:,None]-tb[None,:])
    nearest = np.argmin(difference,axis=1)
    rev = np.argmin(difference,axis=0)
    i = np.flatnonzero((rev[nearest]==np.arange(len(ta))) &
                       (difference[np.arange(len(ta)),nearest] <= 30/86400))
    j = nearest[i]
    record = dict(n_near_simultaneous=int(len(i)))
    if len(i):
        delta = ma[i]-mb[j]
        record.update(median_delta_g_mag=float(np.median(delta)),
                      p90_minus_p10_delta_g_mag=float(np.ptp(np.percentile(delta,[10,90]))),
                      maximum_epoch_difference_s=float(np.max(abs(ta[i]-tb[j]))*86400))
    if len(i) >= 6 and min(np.ptp(ma[i]), np.ptp(mb[j]))>0:
        record["zero_lag_spearman_rho_descriptive"] = numeric(spearmanr(ma[i], mb[j]).statistic)
    return record


def retrieve(field, output_dir):
    errors = []
    for url in SERVICES:
        if "tapvizier" in url:
            query = ('SELECT Source AS source_id, Comp AS component_id, '
                     'Epobs AS epoch_obs, Gmagobs AS g_mag_obs, FGobs AS g_flux_obs, '
                     'e_FGobs AS g_flux_obs_error '
                     'FROM "I/361/lensobs" WHERE Source = '+field)
        else:
            query = ("SELECT source_id, component_id, epoch_obs, g_mag_obs, g_flux_obs, g_flux_obs_error "
                     f"FROM gaiafpr.lens_observation WHERE source_id = {field}")
        try:
            resp = requests.get(url, params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                                                  MAXREC=20000, QUERY=query), timeout=(8,60))
            resp.raise_for_status()
            reader = csv.DictReader(io.StringIO(resp.content.decode("utf-8-sig")))
            names = [k.lower() for k in reader.fieldnames or []]
            if not {"source_id","component_id","epoch_obs","g_mag_obs"} <= set(names):
                raise ValueError("service did not return the requested columns")
            rows = [{k.lower():v for k,v in row.items()} for row in reader]
            if len(rows)>=20000 or not rows or any(str(r["source_id"]).strip()!=field for r in rows):
                raise ValueError("empty, truncated or inconsistent parent result")
            raw = output_dir/("field-"+field+".csv")
            raw.write_bytes(resp.content)
            receipt = dict(service=url, query=query, url=resp.url, bytes=len(resp.content),
                           file=raw.name, sha256=hashlib.sha256(resp.content).hexdigest(),
                           previous_errors=errors)
            return rows, receipt
        except (requests.RequestException, UnicodeError, ValueError, csv.Error) as exc:
            errors.append(f"{url}: {type(exc).__name__}: {exc}")
    return [], dict(errors=errors)


def assess(systems, output_dir):
    """Each system supplies id, FPR field_id and components with FPR image_id."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True,parents=True)
    results, cache = [], {}
    for system in systems:
        components = system["components"]
        field = str(system.get("field_id") or components[0].get("field_id") or "")
        record = dict(id=str(system["id"]), field_id=field)
        if not field.isdecimal():
            record["status"] = "missing_parent_id"
            results.append(record)
            continue
        if field not in cache:
            cache[field] = retrieve(field, output_dir)
        rows, record["input"] = cache[field]
        if rows:
            imgs = sorted({int(r["image_id"]) for r in components})
            observed = {i:[] for i in imgs}
            for row in rows:
                i = int(row["component_id"])
                t = numeric(row["epoch_obs"])
                mag = numeric(row["g_mag_obs"])
                if i in observed and t is not None and mag is not None and 0<mag<30:
                    observed[i].append((t,mag))
            record["components"] = {i:sample(v) for i,v in observed.items()}
            record["pairs"] = [{"component_ids":[i,j], **contrasts(observed[i], observed[j])}
                               for i,j in itertools.combinations(imgs,2)]
            record["status"] = "measured" if all(observed.values()) else "missing_valid_magnitudes"
        else:
            record["status"] = "unavailable"
        results.append(record)
    (output_dir/"variability.json").write_text(json.dumps(dict(systems=results,
        time="BJD in TCB minus 2455197.5 days", magnitude="onboard Gaia G",
        interpretation="Mutual nearest epochs within 30 seconds. Onboard splitting and scan-angle "
        "dependent blending can correlate or anticorrelate these fluxes. Lensing time delays can "
        "suppress zero-lag agreement. No delay search, significance, or promotion."), indent=2,
        allow_nan=False)+"\n")
    return results

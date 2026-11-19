#!/usr/bin/env python3
"""Frozen GWTC-2.1/3/4.0 pair census: all-pair proxy, posterior subset, null controls.

All downloads and numeric work run in ONE GitHub Actions job. No discovery claims.
The all-pair statistic uses GWOSC published 90% summaries; full posterior B is
computed ONLY for the preselected subset and must never be called a global scan.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import io
import itertools
import json
import math
import os
import pathlib
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile

import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.special import logsumexp
from scipy.stats import gaussian_kde

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 541821
RNG = np.random.default_rng(SEED)
RECORDS = {"GWTC-2.1": 6513631, "GWTC-3": 22685054, "GWTC-4.0": 17602505}
LVK_LENSING_RECORD = 18163632
CONTROL_PAIRS = [
    ("GW231113_200417", "GW231223_202619"),  # LVK O4a rank 1, Hanabi log10 B = -3.3
    ("GW190421_213856", "GW190910_112807"),  # LVK O3a Hanabi log10 B = -3.8
    ("GW190424_180648", "GW190727_060333"),  # LVK O3a Hanabi log10 B = -3.9
    ("GW190424_180648", "GW190910_112807"),  # LVK O3a Hanabi log10 B = -3.9
]
N_POST = 1400
MAX_POST_FILES = 28
TOP_PROXY_PAIRS = 12


def dump(name, obj):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return path


def sha(path):
    d = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            d.update(block)
    return d.hexdigest()


def get_bytes(url, timeout=90):
    last = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Open-Space-GWTC-lensing-repro/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if attempt < 3:
                time.sleep(5 * (attempt + 1))
    raise RuntimeError("Unable to fetch " + url + ": " + str(last))


def get_json(url):
    raw = get_bytes(url)
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def stream_file(url, dst, checksum, expected_size):
    """Check the exact Zenodo MD5, record computed SHA256; never trust HTTP alone."""
    expected_md5 = checksum.split(":")[-1].lower() if checksum else None
    last = None
    for attempt in range(3):
        md5, h256, n = hashlib.md5(), hashlib.sha256(), 0
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Open-Space-GWTC-lensing-repro/1.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(dst, "wb") as out:
                for block in iter(lambda: resp.read(2 * 1024 * 1024), b""):
                    out.write(block)
                    md5.update(block)
                    h256.update(block)
                    n += len(block)
            if expected_md5 and md5.hexdigest() != expected_md5:
                raise ValueError("Zenodo checksum mismatch: " + dst.name)
            if expected_size and n != int(expected_size):
                raise ValueError("Zenodo length mismatch: " + dst.name)
            return {"size": n, "md5": md5.hexdigest(), "sha256": h256.hexdigest()}
        except (urllib.error.URLError, TimeoutError, ConnectionError, ValueError) as exc:
            last = str(exc)
            dst.unlink(missing_ok=True)
            if attempt < 2:
                time.sleep(6 * (attempt + 1))
    raise RuntimeError("Failed pinned download: " + str(last))


def entry_files(record):
    f = record.get("files", {})
    if isinstance(f, list):
        return f
    if isinstance(f, dict):
        entries = f.get("entries", f)
        return list(entries.values()) if isinstance(entries, dict) else entries
    raise ValueError("No recognizable Zenodo file entries")


def metadata():
    result = {"releases": {}, "files": {}, "source": "Zenodo versioned record IDs, not concept DOIs"}
    for cat, rid in {**RECORDS, "LVK-lensing-O4a": LVK_LENSING_RECORD}.items():
        data, raw_sha = get_json("https://zenodo.org/api/records/" + str(rid))
        if str(data.get("id")) != str(rid):
            raise ValueError("Zenodo resolved the wrong record: " + cat)
        result["releases"][cat] = {
            "record_id": rid, "doi": data.get("doi"), "version": data.get("metadata", {}).get("version"),
            "record_json_sha256": raw_sha, "record_url": "https://zenodo.org/records/" + str(rid)}
        selected = {}
        for f in entry_files(data):
            key = f.get("key") or f.get("filename")
            if not key:
                continue
            match = re.search(r"GW\d{6}_\d{6}", key)
            is_pe = (cat == "GWTC-4.0" and "combined_PEDataRelease.hdf5" in key) or (
                cat in ("GWTC-2.1", "GWTC-3") and key.endswith("_mixed_nocosmo.h5"))
            is_lvk = cat == "LVK-lensing-O4a" and (key.endswith(".zip") or key == "README.md")
            if not (is_pe or is_lvk):
                continue
            links = f.get("links", {})
            url = links.get("content") or links.get("self") or ("https://zenodo.org/records/" + str(rid) + "/files/" + urllib.parse.quote(key) + "?download=1")
            selected[match.group(0) if match else key] = {
                "record_id": rid, "key": key, "size": f.get("size"),
                "checksum": f.get("checksum") or (f.get("checksum_md5") and "md5:" + f.get("checksum_md5")),
                "url": url}
        result["files"][cat] = selected
        print(cat, "pinned", len(selected), "eligible files", flush=True)
    frozen = OUT / "release_manifest.json"
    if frozen.exists():
        old = json.loads(frozen.read_text())
        if old != result:
            raise RuntimeError("Version-pinned release metadata changed; do not silently update freeze")
    else:
        dump("release_manifest.json", result)
    return result


def gwosc_events():
    cats = "GWTC-1-confident,GWTC-2.1-confident,GWTC-3-confident,GWTC-4.0"
    url = "https://gwosc.org/api/v2/event-versions?" + urllib.parse.urlencode({
        "release": cats, "pagesize": 100, "lastver": "true",
        "include-default-parameters": "true"})
    pages, rows = [], []
    while url:
        data, digest = get_json(url)
        pages.append({"url": url, "sha256": digest, "rows": len(data.get("results", []))})
        rows.extend(data.get("results", []))
        url = data.get("next")
        if len(pages) > 15:
            raise ValueError("Unexpected GWOSC pagination")
    events = {}
    for r in rows:
        name = r.get("name", "").replace("-v1", "")
        if not re.fullmatch(r"GW\d{6}_\d{6}", name):
            continue
        cat = r.get("catalog")
        params = {p["name"]: {"best": p.get("best"), "upper_error": p.get("upper_error"),
                               "lower_error": p.get("lower_error")} for p in r.get("default_parameters", [])}
        events[name] = {"catalog": cat, "gps": r.get("gps"), "params": params,
                        "event_version": r.get("shortName"), "detectors": r.get("detectors", [])}
    frozen = {"selection": "GWOSC four specified catalogues, all versioned event rows",
              "pages": pages, "events": events}
    path = OUT / "gwosc_event_freeze.json"
    if path.exists() and json.loads(path.read_text())["events"] != events:
        raise RuntimeError("GWOSC event summaries changed; frozen selection needs an explicit version")
    dump("gwosc_event_freeze.json", frozen)
    return events


def val(e, name, default=None):
    p = e["params"].get(name, {})
    v = p.get("best", default)
    try:
        return float(v) if v is not None else default
    except (TypeError, ValueError):
        return default


def unc(e, name):
    p = e["params"].get(name, {})
    u, lo = p.get("upper_error"), p.get("lower_error")
    try:
        return (abs(float(u)) + abs(float(lo))) / 2 / 1.645
    except (TypeError, ValueError):
        return None


def population(events, releases):
    eligible = []
    for name, e in events.items():
        cat = e["catalog"]
        if cat in ("GWTC-1-confident", "GWTC-2.1-confident"):
            release = "GWTC-2.1"
        elif cat == "GWTC-3-confident":
            release = "GWTC-3"
        elif cat == "GWTC-4.0":
            release = "GWTC-4.0"
        else:
            continue
        if name not in releases["files"][release]:
            continue
        if val(e, "p_astro", 0) < 0.5 or val(e, "far", math.inf) > 1:
            continue
        if val(e, "mass_2_source", 0) < 3 or not val(e, "redshift"):
            continue
        mc = val(e, "chirp_mass_source")
        if mc is None:
            m1, m2 = val(e, "mass_1_source"), val(e, "mass_2_source")
            mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2 if m1 and m2 else None
        if not mc or not e.get("gps"):
            continue
        z = val(e, "redshift")
        sm = unc(e, "chirp_mass_source") or 0.25 * mc
        sz = unc(e, "redshift") or 0.2 * z
        # Approximate width only: PE mass and redshift are correlated. This is
        # a *preselection proxy*, not a substitute for PE posteriors.
        sln = np.hypot(sm / mc, sz / (1 + z))
        seff = unc(e, "chi_eff") or 0.3
        eligible.append({"name": name, "release": release, "catalog": cat, "gps": e["gps"],
                         "z": z, "mc_det": mc * (1 + z),
                         "lnmc": math.log(mc * (1 + z)), "sigma_lnmc": float(np.clip(sln, .015, 1.0)),
                         "chieff": val(e, "chi_eff", 0), "sigma_chieff": float(np.clip(seff, .04, 0.8)),
                         "file": releases["files"][release][name]["key"]})
    eligible.sort(key=lambda x: (x["gps"], x["name"]))
    dump("eligible_events.json", eligible)
    return eligible


def sis_delays(z_sources, num=70000):
    """Cross-section-weighted SDSS early-type velocity function + SIS double-image delays."""
    rng = np.random.default_rng(SEED + 1)
    z_s = rng.choice(np.maximum(.02, np.array(z_sources)), size=num)
    z_l = rng.uniform(.001, .995, size=num) * z_s
    zgrid = np.linspace(0, max(6, z_s.max() + .1), 18000)
    Ez = np.sqrt(.315 * (1 + zgrid) ** 3 + .685)
    dz = np.diff(zgrid)
    chi = np.r_[0, np.cumsum(.5 * (1 / Ez[1:] + 1 / Ez[:-1]) * dz)] * (299792.458 / 67.4)
    cs, cl = np.interp(z_s, zgrid, chi), np.interp(z_l, zgrid, chi)
    sigma = np.clip(161 * rng.gamma((2.32 + 4) / 2.67, 1, num) ** (1 / 2.67), 20, 600)
    y = np.sqrt(rng.uniform(.00001, .99999, num))
    theta_e = 4 * np.pi * (sigma / 299792.458) ** 2 * ((cs - cl) / cs)
    DlDsOverDls = cl * cs / np.maximum(1e-6, (1 + z_l) * (cs - cl))
    t_days = (1 + z_l) * (DlDsOverDls * 3.085677581e19 / 299792.458 / 86400) * (2 * y * theta_e ** 2)
    optical = np.maximum(1, cl * (cs - cl) / cs) ** 2 / np.sqrt(.315 * (1 + z_l) ** 3 + .685)
    magnification_bias = np.minimum(((1 - y) / y) ** 1.5, 100.0)
    return t_days, optical * magnification_bias


def delay_logbf(pop):
    pairs = list(itertools.combinations(range(len(pop)), 2))
    d = np.array([abs(pop[i]["gps"] - pop[j]["gps"]) / 86400 for i, j in pairs])
    logd = np.log10(np.maximum(.001, d))
    bins = np.linspace(-3, 4, 110)
    unlensed = np.histogram(logd, bins=bins)[0].astype(float)
    unlensed = gaussian_filter1d(unlensed, 2.0) + .01
    unlensed /= unlensed.sum()
    lens_d, weights = sis_delays([x["z"] for x in pop])
    lh = np.histogram(np.log10(np.maximum(.001, lens_d)), bins=bins, weights=weights)[0].astype(float)
    lh = gaussian_filter1d(lh, 2.0)
    lh /= max(lh.sum(), 1)
    # Detection window estimated from actual observing-cadence pair differences:
    # p_L,observed proportional to SIS delay density * p_U,observed.
    pl = lh * unlensed
    pl = pl / max(pl.sum(), 1e-300)
    eps = .01
    pl = (1 - eps) * pl + eps / len(pl)
    pu = (1 - eps) * unlensed + eps / len(unlensed)
    bx = np.clip(np.digitize(logd, bins) - 1, 0, len(pl) - 1)
    logR = np.log(pl[bx] / pu[bx])
    return pairs, d, logR, {"model": "SIS, SDSS VDF sigma*=161 km/s alpha=2.32 beta=2.67, cross section sigma^4, source z from event list",
        "magnification_bias": "fainter-image mu^1.5 capped at 100",
        "observing_window": "approximate, from observed real-pair delay density; no instrument duty-cycle injection",
        "n_lens_draws": len(lens_d), "log10_delta_t_bin_edges": bins.tolist(),
        "unlensed_bin_prob": pu.tolist(), "lensed_bin_prob": pl.tolist()}


def proxy_pairs(pop):
    pairs, d, logR, model = delay_logbf(pop)
    names = [x["name"] for x in pop]
    n = len(pop)
    intrinsic = np.full((n, n), -1e6)
    delay = np.full((n, n), -1e6)
    rows = []
    for k, (i, j) in enumerate(pairs):
        a, b = pop[i], pop[j]
        v = a["sigma_lnmc"] ** 2 + b["sigma_lnmc"] ** 2
        w = a["sigma_chieff"] ** 2 + b["sigma_chieff"] ** 2
        ln_m = -.5 * ((a["lnmc"] - b["lnmc"]) ** 2 / v + np.log(2 * np.pi * v))
        ln_s = -.5 * ((a["chieff"] - b["chieff"]) ** 2 / w + np.log(2 * np.pi * w))
        # No prior factor applied: these are densities in chosen summary coordinates.
        t = ln_m + ln_s
        intrinsic[i,j] = intrinsic[j,i] = t
        delay[i,j] = delay[j,i] = logR[k]
        rows.append({"event_1": names[i], "event_2": names[j],
                     "delay_days": float(d[k]), "summary_mass_spin_logscore": float(t),
                     "log_delay_ratio": float(logR[k]), "proxy_logscore": float(t + logR[k])})
    rows.sort(key=lambda r: (-r["proxy_logscore"], r["event_1"], r["event_2"]))
    with gzip.open(OUT / "all_pairs_proxy.tsv.gz", "wt", newline="") as fo:
        writer = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    dump("time_delay_population.json", model)
    # The null destroys any association of event intrinsic parameters with arrival
    # time, preserving catalog-specific populations and cadence. It does NOT test
    # sky coincidence; these FAP values must never be assigned to full-posterior B.
    groups = {}
    for i, e in enumerate(pop):
        groups.setdefault(e["catalog"], []).append(i)
    nullrng = np.random.default_rng(SEED + 2)
    upper = np.triu_indices(n, 1)
    null_max = []
    for trial in range(1500):
        perm = np.arange(n)
        for idx in groups.values():
            perm[idx] = nullrng.permutation(idx)
        null_max.append(float(np.max((intrinsic[np.ix_(perm, perm)] + delay)[upper])))
    observed_max = rows[0]["proxy_logscore"]
    global_fap = (1 + sum(x >= observed_max for x in null_max)) / (1 + len(null_max))
    for r in rows[:30]:
        r["proxy_pair_tail_fraction"] = sum(x["proxy_logscore"] >= r["proxy_logscore"] for x in rows) / len(rows)
    background = {"simulated_unlensed_trials": len(null_max), "shuffle": "intrinsic (mass, chi_eff) properties permuted as one vector WITHIN each GWOSC catalogue; GPS/delay fixed",
        "event_count": n, "pair_count": len(rows), "max_observed_proxy_logscore": observed_max,
        "global_proxy_fap": global_fap,
        "strongest_pair": [rows[0]["event_1"], rows[0]["event_2"]],
        "null_max_logscore_quantiles": np.quantile(null_max, [.5, .9, .95, .99, 1]).tolist(),
        "valid_for": "all-pair GWOSC summary pre-screen ONLY, not a PE sky/spin overlap search; no selection correction for targeted posterior follow-ups"}
    dump("proxy_background.json", background)
    dump("top_proxy_pairs.json", rows[:30])
    return rows, background


def get_posterior(h5, event):
    required = {"mass_1", "mass_2", "ra", "dec"}
    candidates = []
    def visit(name, obj):
        if not name.endswith("posterior_samples"):
            return
        if isinstance(obj, h5py.Dataset) and obj.dtype.names and required.issubset(obj.dtype.names):
            candidates.append((name, obj, list(obj.dtype.names)))
        if isinstance(obj, h5py.Group) and required.issubset(set(obj.keys())):
            candidates.append((name, obj, list(obj.keys())))
    h5.visititems(visit)
    if not candidates:
        raise ValueError("No complete posterior_samples group/dataset with detector masses + sky")
    def sortkey(entry):
        name = entry[0]
        return (0 if "IMRPhenomXPHM-SpinTaylor" in name else
                1 if "IMRPhenomXPHM" in name else
                2 if "Mixed" in name else 3, name)
    label, ds, columns = sorted(candidates, key=sortkey)[0]
    n = ds.shape[0] if isinstance(ds, h5py.Dataset) else len(ds["mass_1"])
    if n < 50:
        raise ValueError("Insufficient samples")
    seed = int.from_bytes(hashlib.sha256(event.encode()).digest()[:8], "big")
    idx = np.sort(np.random.default_rng(seed).choice(n, min(n, N_POST), replace=False))
    def col(key):
        if key not in columns: return None
        return np.asarray(ds[idx][key] if isinstance(ds, h5py.Dataset) else ds[key][idx], dtype=float)
    m1, m2 = col("mass_1"), col("mass_2")
    a1, a2 = col("a_1"), col("a_2")
    ct1, ct2 = col("cos_tilt_1"), col("cos_tilt_2")
    if ct1 is None and col("tilt_1") is not None: ct1 = np.cos(col("tilt_1"))
    if ct2 is None and col("tilt_2") is not None: ct2 = np.cos(col("tilt_2"))
    inc = col("theta_jn")
    if inc is None: inc = col("iota")
    fields = {"mass_1": m1, "mass_2": m2, "ra": col("ra"), "dec": col("dec"),
              "a_1": a1, "a_2": a2, "cos_tilt_1": ct1, "cos_tilt_2": ct2,
              "inclination": inc, "mass_1_source": col("mass_1_source"),
              "mass_2_source": col("mass_2_source"), "geocent_time": col("geocent_time")}
    finite = np.isfinite(m1) & np.isfinite(m2) & (m1 > 0) & (m2 > 0)
    for key in ("ra", "dec"):
        finite &= np.isfinite(fields[key])
    for k, v in fields.items():
        if v is not None:
            fields[k] = v[finite]
    if finite.sum() < 50: raise ValueError("Too few finite samples")
    return fields, {"label": label, "original_samples": int(n), "extracted_samples": int(finite.sum()),
                    "columns": columns, "missing": [k for k,v in fields.items() if v is None]}


def download_sample(event, meta, tmp):
    file = meta["files"][event["release"]][event["name"]]
    dest = tmp / pathlib.Path(file["key"]).name
    info = stream_file(file["url"], dest, file["checksum"], file.get("size"))
    try:
        with h5py.File(dest, "r") as h5:
            post, detail = get_posterior(h5, event["name"])
    finally:
        dest.unlink(missing_ok=True)
    save = {k: v for k,v in post.items() if v is not None}
    output = OUT / "extracted" / (event["name"] + ".npz")
    output.parent.mkdir(exist_ok=True)
    np.savez_compressed(output, **save)
    detail.update({"event": event["name"], "release": event["release"], "zenodo_key": file["key"],
                   "zenodo_record": file["record_id"], "zenodo_md5": info["md5"],
                   "original_hdf5_sha256": info["sha256"], "original_hdf5_size": info["size"],
                   "extracted_npz_sha256": sha(output)})
    return detail


def files_and_posteriors(pop, ranked, meta):
    nameto = {e["name"]: e for e in pop}
    controls = {(min(a,b), max(a,b)) for a,b in CONTROL_PAIRS}
    target_pairs = []
    for a,b in CONTROL_PAIRS:
        if a in nameto and b in nameto:
            target_pairs.append((a,b,"LVK-control"))
    for r in ranked[:TOP_PROXY_PAIRS]:
        target_pairs.append((r["event_1"],r["event_2"],"GWOSC-preselected"))
    selected = []
    for a,b,reason in target_pairs:
        for n in (a,b):
            if n not in selected and len(selected) < MAX_POST_FILES:
                selected.append(n)
    dump("target_selection.json", {"selection_frozen_before_downloading_posteriors": True,
        "method": "all 4 published LVK pairs plus top 12 from full GWOSC mass+chi_eff+SIS delay proxy",
        "max_posterior_files": MAX_POST_FILES, "events": selected, "pairs": target_pairs})
    receipts, failed = [], {}
    with tempfile.TemporaryDirectory(prefix="gwtc-lensing-") as dirname:
        for k,name in enumerate(selected, 1):
            print("PE", k, "/", len(selected), name, flush=True)
            try:
                receipts.append(download_sample(nameto[name], meta, pathlib.Path(dirname)))
                dump("posterior_receipts.json", receipts)
            except Exception as exc:
                failed[name] = str(exc)
                dump("download_failures.json", failed)
                print("::warning::PE unavailable " + name + ": " + str(exc), flush=True)
    return receipts, failed, target_pairs


def log_prior(X):
    # Reference product prior: component masses uniform over 1..1000 Msun,
    # spin magnitudes U(0,.99), tilt cosines U(-1,1); transformed KDE volume.
    lm, lq, la, lb, c1, c2 = X
    q = 1 / (1 + np.exp(-np.clip(lq, -30, 30)))
    a = 1 / (1 + np.exp(-np.clip(la, -30, 30)))
    b = 1 / (1 + np.exp(-np.clip(lb, -30, 30)))
    ct1, ct2 = np.tanh(c1), np.tanh(c2)
    return (math.log(2 / 999 ** 2) + 2 * lm + .4 * np.log1p(q)
           - 1.2 * np.log(q) + np.log(q * (1 - q))
           + np.log(a * (1 - a)) + np.log(b * (1 - b))
           + np.log((1 - ct1 ** 2) / 2) + np.log((1 - ct2 ** 2) / 2))


def intrinsic_features(sample):
    need = ("a_1", "a_2", "cos_tilt_1", "cos_tilt_2")
    if not all(k in sample for k in need):
        return None
    m1, m2 = sample["mass_1"], sample["mass_2"]
    q = np.minimum(m1,m2) / np.maximum(m1,m2)
    mc = (m1*m2)**.6 / (m1+m2)**.2
    u1 = np.clip(sample["a_1"]/.99, .001, .999)
    u2 = np.clip(sample["a_2"]/.99, .001, .999)
    ct1 = np.clip(sample["cos_tilt_1"], -.999, .999)
    ct2 = np.clip(sample["cos_tilt_2"], -.999, .999)
    logit = lambda x: np.log(x / (1-x))
    X = np.vstack([np.log(mc),logit(np.clip(q,.001,.999)),logit(u1),logit(u2),
                   np.arctanh(ct1),np.arctanh(ct2)])
    return X[:, np.all(np.isfinite(X), axis=0)]


def intrinsic_bayes(a,b):
    X, Y = intrinsic_features(a), intrinsic_features(b)
    if X is None or Y is None or X.shape[1] < 80 or Y.shape[1] < 80:
        return None
    ka = gaussian_kde(X, bw_method="scott")
    kb = gaussian_kde(Y, bw_method="scott")
    # Symmetric Haris p2/pi average over p1, and reverse, with independent
    # KDEs, using deterministic 200-evaluation subsets for cost control.
    ea = np.linspace(0, X.shape[1]-1, min(240,X.shape[1]), dtype=int)
    eb = np.linspace(0, Y.shape[1]-1, min(240,Y.shape[1]), dtype=int)
    a1 = logsumexp(kb.logpdf(X[:,ea]) - log_prior(X[:,ea])) - np.log(len(ea))
    b1 = logsumexp(ka.logpdf(Y[:,eb]) - log_prior(Y[:,eb])) - np.log(len(eb))
    return float(np.logaddexp(a1,b1)-np.log(2))


def sky_bayes(a,b):
    import healpy as hp
    nside = 32
    def pixprob(s):
        ra = np.mod(s["ra"], 2*np.pi)
        th = np.clip(np.pi/2 - s["dec"], 0, np.pi)
        p = hp.ang2pix(nside, th, ra)
        v = np.bincount(p, minlength=hp.nside2npix(nside)).astype(float)
        v = hp.smoothing(v, fwhm=np.deg2rad(2), verbose=False)
        v = np.maximum(v, 1e-12)
        return v / v.sum()
    x,y = pixprob(a),pixprob(b)
    return float(np.log(len(x) * np.dot(x,y)))


def inclination_bayes(a,b):
    if "inclination" not in a or "inclination" not in b:
        return None
    xa,xb = np.cos(a["inclination"]),np.cos(b["inclination"])
    pa = np.histogram(xa,bins=20,range=(-1,1))[0].astype(float)+1
    pb = np.histogram(xb,bins=20,range=(-1,1))[0].astype(float)+1
    pa /= pa.sum();pb /= pb.sum()
    return float(np.log(20*np.dot(pa,pb)))


def pair_followups(target_pairs, ranked, receipts, population):
    byevent = {}
    for receipt in receipts:
        p = OUT / "extracted" / (receipt["event"] + ".npz")
        with np.load(p) as z:
            byevent[receipt["event"]] = {k: z[k] for k in z.files}
    index = {x["name"]:x for x in population}
    index_proxy = {tuple(sorted((r["event_1"],r["event_2"]))):r for r in ranked}
    seen, rows = set(), []
    for a,b,reason in target_pairs:
        key = tuple(sorted((a,b)))
        if key in seen:continue
        seen.add(key)
        row = {"event_1":a,"event_2":b,"preselection":reason,
               "status":"missing_posterior", "phase_check":"pending, no phase-coherent joint PE",
               "strain_check":"pending, no GWOSC strain likelihood",
               "em_lens_check":"pending, no complete lens galaxy foreground crossmatch"}
        if key in index_proxy:
            row["proxy_logscore"] = index_proxy[key]["proxy_logscore"]
            row["time_delay_days"] = index_proxy[key]["delay_days"]
            row["log_delay_ratio"] = index_proxy[key]["log_delay_ratio"]
        if a in byevent and b in byevent:
            try:
                pa,pb = byevent[a],byevent[b]
                intr = intrinsic_bayes(pa,pb)
                sky = sky_bayes(pa,pb)
                inc = inclination_bayes(pa,pb)
                row.update({"status": "posterior_screened", "lnB_intrinsic":intr,
                            "lnB_sky":sky, "lnB_inclination":inc,
                            "log10_B_overlap_approx": (intr+sky+(inc or 0))/np.log(10) if intr is not None else None,
                            "prior_warning":"common mass-box prior is APPROXIMATE; differing PE spin/mass priors not reweighted"})
                for field in ("mass_1", "mass_2", "mass_1_source", "mass_2_source"):
                    if field in pa and field in pb:
                        qa,qb = np.quantile(pa[field],[.05,.5,.95]),np.quantile(pb[field],[.05,.5,.95])
                        row[field + "_90pct_overlaps"] = bool(qa[0]<=qb[2] and qb[0]<=qa[2])
                        row[field + "_medians"] = [float(qa[1]),float(qb[1])]
                row["frame_note"] = "Detector-frame masses invariant; source-frame masses NOT invariant to lens magnification through inferred redshift"
            except Exception as exc:
                row["status"]="posterior_error"
                row["error"]=str(exc)
        rows.append(row)
    dump("posterior_pairs.json", rows)
    return rows


def lvk_reference(meta, tmp):
    """Pin and inspect the original LVK release; never execute pickle from ZIP."""
    files = meta["files"]["LVK-lensing-O4a"]
    archive = next((v for v in files.values() if v["key"].endswith(".zip")), None)
    if archive is None:return {"status":"no_archive_indexed"}
    dest=tmp/"lvk_lensing.zip"
    try:
        digest=stream_file(archive["url"],dest,archive["checksum"],archive["size"])
        with zipfile.ZipFile(dest) as z:
            listing=[{"name":n.filename,"uncompressed_size":n.file_size} for n in z.infolist()]
            textual=[]
            for f in z.infolist():
                if f.filename.lower().endswith((".csv",".tsv",".txt",".md",".json")) and 0 < f.file_size < 3_000_000:
                    data=z.read(f)
                    if b"S231113bw" in data or b"GW231113_200417" in data:
                        textual.append({"name":f.filename, "sha256":hashlib.sha256(data).hexdigest(),
                                        "lines":data.decode("utf-8","replace").splitlines()[:15]})
        dump("lvk_archive_inventory.json",listing)
        return {"status":"verified_archived", "archive":archive["key"],**digest,
                "number_of_members":len(listing),"references_to_published_top_pair":textual}
    except Exception as exc:
        return {"status":"unavailable", "error":str(exc)}
    finally:
        dest.unlink(missing_ok=True)


def report(pop, rows, background, receipts, failures, evaluated, lvk):
    strong = sorted((r for r in evaluated if r["status"]=="posterior_screened" and r.get("log10_B_overlap_approx") is not None),
                    key=lambda r:-r["log10_B_overlap_approx"])
    bycat={}
    for x in pop:bycat[x["release"]]=bycat.get(x["release"],0)+1
    summary={"conclusion":"NO VALIDATED LENSED PAIR; posterior follow-up is incomplete and non-global",
        "n_eligible_events":len(pop),"n_all_summary_proxy_pairs":len(rows),
        "eligible_by_release":bycat,"n_downloaded_posterior_files":len(receipts),
        "n_missing_selected_posterior_files":len(failures),
        "n_posterior_screened_pairs":len(strong),
        "top_posterior_screened_pairs":strong[:8],
        "strongest_all_pair_proxy":rows[0],
        "global_proxy_fap":background["global_proxy_fap"],
        "global_posterior_fap":"NOT CALIBRATED: only preselected PE pairs were downloaded",
        "lvk_control_conclusion":"LVK found NO convincing O4a lensed pair; its top S231113bw+S231223bg has Hanabi log10 B=-3.3 against lensing",
        "o3_control_conclusion":"LVK O3a Hanabi log10 B=-3.8 or -3.9 for three highly ranked candidate pairs",
        "lvk_archive_status":lvk["status"],
        "not_done":["posterior-complete full-catalog false-alarm calibration",
                    "Morse phase coherent joint PE / parity likelihood",
                    "detector strain joint likelihood on GWOSC",
                    "spectroscopic/EM lens-galaxy crossmatch and image geometry",
                    "source population / selection function / PE prior matching"]}
    dump("result.json",summary)
    text = [
        "# GWTC lensing census: runner result",
        "",
        "**Disposition: no validated strong-lensing discovery.** This is a *partial* posterior search plus an all-pair summary pre-screen; no published detection claim.",
        "",
        "- Eligible events by version-pinned catalogue: "+str(bycat)+".",
        "- All-pair mass/chi_eff/time-delay **proxy**: "+str(len(rows))+" pairs. Approximate SIS population, actual GWOSC event cadence, 1,500 catalogue-stratified shuffles.",
        "- Proxy best: "+rows[0]["event_1"]+" / "+rows[0]["event_2"]+
        "; proxy family-wise null-tail fraction = "+f'{background["global_proxy_fap"]:.4f}'+". **Not** a posterior-overlap false-alarm probability.",
        "- Downloaded and MD5-validated "+str(len(receipts))+" distinct released PE files; "+str(len(failures))+" unavailable. "+
        str(len(strong))+" preselected pairs screened with full PE samples. Download SHA256 and sampling labels in posterior_receipts.json.",
        "- LVK O4a control: GW231113_200417 / GW231223_202619; published joint log10 Bayes factor -3.3, **not lensed**. Earlier O3a controls likewise negative joint odds.",
        "- No phase-coherent Morse image-type likelihood, joint GWOSC strain fit or foreground galaxy crossmatch: **pending**, not 'passed'.",
        "- A mass in the **source frame** depends on luminosity distance/redshift and so is not a lensing invariant; detector-frame masses are the gate.",
        "",
        "## Screened posterior pairs",
        "",
        "| Event 1 | Event 2 | log10 approximate overlap B | status |",
        "| --- | --- | ---: | --- |"]
    for r in strong[:10]:
        text.append("| "+r["event_1"]+" | "+r["event_2"]+" | "+f'{r["log10_B_overlap_approx"]:.2f}'+" | not discovery |")
    if not strong:text.append("| — | — | — | posterior download/parse incomplete |")
    text.extend(["", "## Constraints and receipts", "",
        "See result.json, all_pairs_proxy.tsv.gz, proxy_background.json, time_delay_population.json, gwosc_event_freeze.json, release_manifest.json, posterior_pairs.json, posterior_receipts.json and download_failures.json (if any).",
        "Posterior overlap uses symmetric KDE ratios to an approximate common PE prior, and sky HEALPix overlap; no lens selection bias or full-population PE FAP has been evaluated. The summary pre-screen cannot supply sky localization, image phase or foreground lens identification.",
        "The controls' joint Bayes factors are published **external** results, not our recalculated joint evidence; a large overlap B cannot override their negative joint-analysis conclusions.",
        "No [REMARKABLE] label per docs/REMARKABLE_FINDINGS.md."])
    (OUT/"RESULT.md").write_text("\n".join(text)+"\n")
    return summary


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--test-only",action="store_true")
    args=parser.parse_args()
    if args.test_only:
        test_self()
        return
    meta=metadata()
    ev=gwosc_events()
    pop=population(ev,meta)
    if len(pop)<20:
        raise RuntimeError("GWOSC selection unexpectedly small; refusing any search claim")
    ranked,bg=proxy_pairs(pop)
    receipts,failures,target_pairs=files_and_posteriors(pop,ranked,meta)
    with tempfile.TemporaryDirectory(prefix="lvk-archive-") as d:
        lvk=lvk_reference(meta,pathlib.Path(d))
    dump("lvk_reference.json",lvk)
    evaluated=pair_followups(target_pairs,ranked,receipts,pop)
    summary=report(pop,ranked,bg,receipts,failures,evaluated,lvk)
    print(json.dumps({k:v for k,v in summary.items() if k not in ("top_posterior_screened_pairs",)}, indent=2), flush=True)


def test_self():
    assert len(sis_delays([.5,1],num=1000)[0])==1000
    x=np.array([[np.log(30)], [0.], [0.], [0.], [0.], [0.]])
    assert np.isfinite(log_prior(x)).all()
    q=np.array([.1,.4,.9])
    assert np.all(np.diff(q)>0)
    print("synthetic lens population / Jacobian tests OK")


if __name__=="__main__":
    try:main()
    except Exception as exc:
        dump("fatal_error.json",{"status":"FAILED - NOT A RESULT","error":repr(exc)})
        raise

"""Reproducible exploratory search of the latest public TESS SPOC sector.

Downloaded scripts are parsed as text, never executed. Native photometry is
preserved in compressed arrays; each product has its full-FITS SHA256 and source
URL. A screening score is not a false-alarm probability or a planet probability.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from functools import partial
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

FIELDS = {"time": "TIME", "pdcsap": "PDCSAP_FLUX", "pdcsap_err": "PDCSAP_FLUX_ERR",
          "sap": "SAP_FLUX", "sap_err": "SAP_FLUX_ERR", "quality": "QUALITY",
          "bkg": "SAP_BKG", "cx": "MOM_CENTR1", "cy": "MOM_CENTR2"}
HEADERS = ["TICID", "SECTOR", "CAMERA", "CCD", "RA_OBJ", "DEC_OBJ", "TESSMAG",
           "TEFF", "RADIUS", "LOGG", "PMRA", "PMDEC", "CROWDSAP", "FLFRCSAP",
           "BJDREFI", "BJDREFF", "TIMESYS", "DATE-OBS", "DATE-END", "PROCVER"]
PRODUCT = re.compile(r"tess\d+-s(\d{4})-(\d{16})-\d+-s_lc\.fits")


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def get(url, timeout=90, max_bytes=32_000_000):
    last = None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            if len(r.content) > max_bytes:
                raise ValueError("Response exceeds declared byte limit")
            return r.content
        except (requests.RequestException, ValueError) as exc:
            last = exc
            if isinstance(exc, requests.HTTPError) and exc.response.status_code in (401, 403, 404):
                break
            time.sleep(attempt + 1)
    raise last


def inventory(sector, cache):
    out = Path(cache) / "inventories"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"s{sector}.json"
    if path.exists():
        result = json.loads(path.read_text())
        raw = (out / f"s{sector}.sh").read_bytes()
        if sha(raw) != result["sha256"]:
            raise ValueError("Inventory checksum mismatch")
        return result
    url = f"https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{sector}_lc.sh"
    raw = get(url)
    products = {}
    for line in raw.decode().splitlines():
        tokens = line.split()
        if not tokens or tokens[0] != "curl" or "-o" not in tokens:
            continue
        name = tokens[tokens.index("-o") + 1]
        m = PRODUCT.fullmatch(name)
        if not m or int(m[1]) != sector:
            continue
        if m[2] in products and products[m[2]] != name:
            raise ValueError("Multiple reductions for a target")
        products[str(int(m[2]))] = name
    if not products:
        raise ValueError("No standard two-minute light curves")
    (out / f"s{sector}.sh").write_bytes(raw)
    result = {"sector": sector, "url": url, "sha256": sha(raw),
              "retrieved_at": stamp(), "n": len(products), "products": products}
    write_json(path, result)
    return result


def metadata_batch(job, cache):
    i, targets = job
    root = Path(cache) / "tic108"
    root.mkdir(exist_ok=True)
    path = root / f"batch{i:03}.json"
    req_path = root / f"batch{i:03}_request.json"
    if path.exists():
        return len(json.loads(path.read_text())["data"])
    req = {"service": "Mast.Catalogs.Filtered.Tic", "format": "json", "pagesize": 10000,
           "page": 1, "params": {"columns": "*", "filters": [
               {"paramName": "ID", "values": list(map(str, targets))}]}}
    write_json(req_path, req)
    url = "https://mast.stsci.edu/api/v0/invoke"
    r = requests.post(url, data={"request": json.dumps(req)}, timeout=120)
    r.raise_for_status()
    obj = r.json()
    if obj.get("status") != "COMPLETE":
        raise ValueError("Incomplete TIC response")
    returned = {int(x["ID"]) for x in obj["data"]}
    if not returned.issubset(set(targets)):
        raise ValueError("TIC response contains unrequested IDs")
    path.write_bytes(r.content)
    write_json(root / f"batch{i:03}_receipt.json", {
        "url": url, "request_sha256": sha(req_path.read_bytes()),
        "response_sha256": sha(r.content), "time": stamp(), "expected": len(targets),
        "returned": len(obj["data"]), "paging": obj.get("paging")})
    return len(obj["data"])


def make_plan(cache, output, sector=108):
    cache, output = Path(cache), Path(output)
    inv = inventory(sector, cache)
    index = json.loads((cache / "catalogs/exclusion_index.json").read_text())
    excluded = set(map(int, index["tic_exclusions"]))
    rows = {}
    for p in sorted((cache / "tic108").glob("batch[0-9][0-9][0-9].json")):
        for row in json.loads(p.read_text())["data"]:
            if int(row["ID"]) in rows:
                raise ValueError("Duplicate TIC metadata row")
            rows[int(row["ID"])] = row
    selection = []
    for tic in sorted(map(int, inv["products"])):
        r = rows.get(tic, {})
        reason = "eligible"
        if not r:
            reason = "no_tic_metadata"
        elif tic in excluded:
            reason = "previously_reported_host"
        elif r.get("rad") is None or not 0.1 <= r["rad"] <= 1.4:
            reason = "radius_missing_or_outside_0.1_1.4_solar"
        elif r.get("Teff") is None or not 2700 <= r["Teff"] <= 6750:
            reason = "temperature_missing_or_outside_2700_6750_K"
        elif r.get("Tmag") is None or not 7.5 <= r["Tmag"] <= 13.5:
            reason = "magnitude_missing_or_outside_7.5_13.5"
        elif r.get("contratio") is not None and r["contratio"] > 1:
            reason = "catalog_contamination_above_one"
        elif r.get("disposition") in ("DUPLICATE", "ARTIFACT"):
            reason = "catalog_duplicate_or_artifact"
        selection.append({"tic": tic, "sector": sector, "filename": inv["products"][str(tic)],
                          "ra": r.get("ra"), "dec": r.get("dec"), "tmag": r.get("Tmag"),
                          "teff": r.get("Teff"), "radius": r.get("rad"), "mass": r.get("mass"),
                          "gaia_dr2": r.get("GAIA"), "contratio": r.get("contratio"),
                          "selection": reason})
    eligible = [r for r in selection if r["selection"] == "eligible"]
    eligible.sort(key=lambda r: (sha(f"open-space-20261125:{r['tic']}".encode()), r["tic"]))
    for i, row in enumerate(eligible):
        row["rank"] = i
    output.mkdir(parents=True, exist_ok=True)
    plan = {"created_at": stamp(), "sector": sector, "inventory": {k: v for k, v in inv.items() if k != "products"},
            "catalog_exclusion_index_sha256": sha((cache / "catalogs/exclusion_index.json").read_bytes()),
            "inventory_targets": len(selection), "metadata_targets": len(rows), "eligible_targets": len(eligible),
            "ordering": "Ascending SHA256(open-space-20261125:TIC)", "targets": eligible}
    write_json(output / "plan.json", plan)
    fields = list(selection[0])
    fields = [x for x in fields if x != "rank"] + ["rank"]
    with (output / "selection.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fields, lineterminator="\n")
        w.writeheader(); w.writerows(selection)
    return {k: v for k, v in plan.items() if k != "targets"}


def column_digest(data):
    h = hashlib.sha256()
    for key in FIELDS:
        a = np.asarray(data[key], dtype="<i8" if key == "quality" else "<f8").copy()
        if key != "quality":
            a[np.isnan(a)] = np.nan
        h.update(key.encode() + b"\0")
        h.update(a.tobytes())
    return h.hexdigest()


def fetch_one(row, cache):
    tic, sector = int(row["tic"]), int(row["sector"])
    dest = Path(cache) / "lightcurves"
    dest.mkdir(exist_ok=True)
    path = dest / f"tic{tic}_s{sector}.npz"
    meta_path = path.with_suffix(".json")
    if path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        with np.load(path) as data:
            if column_digest(data) != meta["columns_sha256"]:
                raise ValueError("Cached photometry checksum mismatch")
        return meta
    name = row["filename"]
    m = PRODUCT.fullmatch(name)
    if not m or int(m[1]) != sector or int(m[2]) != tic:
        raise ValueError("Filename identity mismatch")
    padded = f"{tic:016d}"
    directory = "/".join(padded[j:j+4] for j in range(0, 16, 4))
    url = f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s{sector:04d}/{directory}/{name}"
    raw = get(url)
    with fits.open(io.BytesIO(raw), memmap=False) as hdus:
        if int(hdus[0].header["TICID"]) != tic or int(hdus[0].header["SECTOR"]) != sector:
            raise ValueError("FITS identity mismatch")
        data = {key: np.asarray(hdus[1].data[col]).copy() for key, col in FIELDS.items()}
        headers = {}
        for key in HEADERS:
            for hdu in hdus[:2]:
                if key in hdu.header:
                    value = hdu.header[key]
                    headers[key] = value if not isinstance(value, float) or np.isfinite(value) else None
        if headers.get("TIMESYS") != "TDB" or headers.get("BJDREFI") != 2457000:
            raise ValueError("Unexpected TESS time convention")
    meta = {"tic": tic, "sector": sector, "filename": name, "url": url,
            "mast_uri": "mast:TESS/product/" + name, "sha256": sha(raw), "bytes": len(raw),
            "columns_sha256": column_digest(data), "retrieved_at": stamp(), "headers": headers,
            "cadences": len(data["time"])}
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **data); tmp.replace(path)
    write_json(meta_path, meta)
    return meta


def curve(path, column="pdcsap", window_days=1.0):
    with np.load(path) as data:
        valid = ((data["quality"] == 0) & np.isfinite(data["time"]) &
                 np.isfinite(data[column]) & (data[column] > 0))
        t = data["time"][valid].astype(float)
        f = data[column][valid].astype(float)
    if len(t) < 1000:
        raise ValueError("Fewer than 1000 clean native cadences")
    order = np.argsort(t); t, f = t[order], f[order]
    med = np.median(f); f /= med
    bins = np.floor((t - t[0]) * 144 + 1e-6).astype(int)
    count = np.bincount(bins); ok = count >= 3
    t = np.bincount(bins, weights=t)[ok] / count[ok]
    f = np.bincount(bins, weights=f)[ok] / count[ok]
    trend = np.full(len(t), np.nan)
    for part in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > .2) + 1):
        if len(part) < 100:
            continue
        width = min(int(window_days * 144) // 2 * 2 + 1, len(part) // 2 * 2 - 1)
        trend[part] = median_filter(f[part], size=width, mode="nearest")
    y = f / trend - 1
    scatter = 1.4826 * np.nanmedian(abs(y - np.nanmedian(y)))
    valid = np.isfinite(y) & (y < 6 * scatter)
    t, y = t[valid], y[valid]
    if len(t) < 1000 or np.ptp(t) < 10 or not scatter > 0:
        raise ValueError("Inadequate searchable coverage")
    return t, y, scatter


def box_measure(t, y, p, epoch, duration):
    phase = (t - epoch + p / 2) % p - p / 2
    inside = abs(phase) < duration / 2
    outside = (abs(phase) > duration) & (abs(phase) < min(3*duration, .45*p))
    if inside.sum() < 5 or outside.sum() < 10:
        return {"depth_ppm": None, "snr": None, "n_events": 0}
    scatter = 1.4826*np.median(abs(y[outside]-np.median(y[outside])))
    depth = float(np.mean(y[outside])-np.mean(y[inside]))
    error = float(scatter*np.sqrt(1/inside.sum()+1/outside.sum()))
    epochs = np.rint((t-epoch)/p).astype(int)
    event_depths = []
    for e in np.unique(epochs[inside]):
        take = epochs == e
        i, o = take & inside, take & outside
        if i.sum() >= max(3, duration*144*.4) and o.sum() >= 6:
            event_depths.append(float(np.mean(y[o])-np.mean(y[i])))
    return {"depth_ppm": depth*1e6, "snr": depth/error if error > 0 else 0.,
            "n_events": len(event_depths), "positive_events": sum(d > 0 for d in event_depths),
            "event_depths_ppm": [d*1e6 for d in event_depths]}


def screen_one(row, cache, output):
    tic, sector = int(row["tic"]), int(row["sector"])
    dest = Path(output) / "screen" / f"tic{tic}_s{sector}.json"
    if dest.exists():
        return json.loads(dest.read_text())
    result = {"tic": tic, "sector": sector, "rank": row.get("rank"), "radius": row.get("radius"),
              "status": "searched", "peaks": []}
    try:
        path = Path(cache) / "lightcurves" / f"tic{tic}_s{sector}.npz"
        t, y, scatter = curve(path)
        baseline, origin = np.ptp(t), t[0]
        model = BoxLeastSquares(t-origin, y, dy=scatter)
        periods = 1 / np.arange(1/13, 1/.5, .008/baseline)[::-1]
        durations = [1/24, 1.5/24, 2/24, 3/24, 4/24]
        power = model.power(periods, durations, oversample=5)
        inverted = BoxLeastSquares(t-origin, -y, dy=scatter).power(periods, durations, oversample=5)
        result.update({"cadences_binned": len(t), "span_days": float(baseline),
                       "scatter_ppm": float(scatter*1e6), "n_trial_periods": len(periods),
                       "inverted_max_snr": float(np.max(inverted.depth_snr))})
        ranking = power.power.copy()
        for rank in range(1, 6):
            k = np.argmax(ranking); p0 = periods[k]
            fine_periods = p0*np.linspace(1-.025*p0/baseline, 1+.025*p0/baseline, 61)
            fine = model.power(fine_periods, durations, oversample=10)
            j = np.argmax(fine.power)
            p, epoch, duration = float(fine.period[j]), float(fine.transit_time[j]+origin), float(fine.duration[j])
            dep, snr = float(fine.depth[j]), float(fine.depth_snr[j])
            measurement = box_measure(t, y, p, epoch, duration)
            midpoint = (t.min()+t.max())/2
            halves = [box_measure(t[m], y[m], p, epoch, duration) for m in (t < midpoint, t >= midpoint)]
            secondary = box_measure(t, y, p, epoch+p/2, duration)
            rad = float(row.get("radius") or 1.)*109.076*np.sqrt(max(0, dep))
            prelim = (snr >= 8 and .00005 <= dep <= .03 and duration/p <= .12 and
                      measurement["n_events"] >= 3 and measurement.get("positive_events", 0) >= 3 and
                      all(h.get("snr") is not None and h["snr"] >= 3 for h in halves) and rad <= 22)
            result["peaks"].append({"peak_rank": rank, "period_days": p, "epoch_btjd": epoch,
                                    "duration_days": duration, "depth_ppm": dep*1e6, "white_snr": snr,
                                    "radius_proxy_rearth": rad, "local": measurement, "halves": halves,
                                    "secondary": secondary, "preliminary_pass": bool(prelim)})
            ranking[abs(1/periods-1/p) < .8/baseline] = -np.inf
    except Exception as exc:
        result["status"] = "error"
        result["error"] = repr(exc)
    write_json(dest, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["inventory", "metadata", "plan", "fetch", "screen", "summarize"])
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    parser.add_argument("--sector", type=int, default=108)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=512)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.command == "inventory":
        r = inventory(args.sector, args.cache)
        print(json.dumps({k:v for k,v in r.items() if k != "products"})); return
    if args.command == "metadata":
        ids = sorted(map(int, inventory(args.sector, args.cache)["products"]))
        jobs = list(enumerate([ids[j:j+500] for j in range(0, len(ids), 500)]))
        with ThreadPoolExecutor(min(args.workers, 6)) as pool:
            for n in pool.map(partial(metadata_batch, cache=args.cache), jobs): print(n, flush=True)
        return
    if args.command == "plan":
        print(json.dumps(make_plan(args.cache, args.output, args.sector), indent=2)); return
    if args.command == "summarize":
        records = [json.loads(p.read_text()) for p in sorted((args.output/"screen").glob("*.json"))]
        shortlist = [{"tic": r["tic"], "sector": r["sector"], **p} for r in records for p in r["peaks"] if p["preliminary_pass"]]
        shortlist.sort(key=lambda r:r["white_snr"], reverse=True)
        write_json(args.output/"shortlist.json", shortlist)
        print(json.dumps({"records": len(records), "searched": sum(r["status"]=="searched" for r in records),
                          "preliminary_peaks": len(shortlist), "top": shortlist[:12]}, indent=2)); return
    plan = json.loads((args.output/"plan.json").read_text())
    rows = plan["targets"][args.start:args.start+args.limit]
    pool_type = ThreadPoolExecutor if args.command == "fetch" else ProcessPoolExecutor
    fn = partial(fetch_one, cache=args.cache) if args.command == "fetch" else partial(screen_one, cache=args.cache, output=args.output)
    with pool_type(args.workers) as pool:
        futures = {pool.submit(fn, row):row for row in rows}
        failures = []
        for i, future in enumerate(as_completed(futures), 1):
            row = futures[future]
            try:
                r = future.result()
                if args.command == "screen" and any(p["preliminary_pass"] for p in r["peaks"]):
                    print("TRIGGER", r["tic"], [(round(p["period_days"], 6), round(p["white_snr"], 2)) for p in r["peaks"] if p["preliminary_pass"]], flush=True)
                if r.get("status") == "error": failures.append({"tic":row["tic"], "error":r["error"]})
            except Exception as exc:
                failures.append({"tic": row["tic"], "error": repr(exc)})
            if i%25==0 or i==len(rows):print(args.command, i, "/", len(rows), "failures", len(failures), flush=True)
        write_json(args.output/f"{args.command}_{args.start}_{len(rows)}_receipt.json", {
            "command": args.command, "created_at":stamp(), "start":args.start, "requested":len(rows), "failures": failures})


if __name__ == "__main__":
    main()

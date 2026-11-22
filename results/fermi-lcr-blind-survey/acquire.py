"""Read the public Fermi-LAT Light Curve Repository, inside GitHub Actions."""
import argparse
import concurrent.futures
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import gzip
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import threading
import time
from urllib.parse import quote, urlencode

import numpy as np
import requests

BASE = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/"
CADENCES = {"3day": "daily", "weekly": "weekly", "monthly": "monthly"}
THREAD = threading.local()


def slug(name):
    return name.replace(" ", "_")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def session():
    if not hasattr(THREAD, "session"):
        THREAD.session = requests.Session()
        THREAD.session.headers["User-Agent"] = "Open-Space-LCR-research/1.0 (github.com/dovvnloading/Open-Space)"
    return THREAD.session


class ArchiveDeferred(Exception):
    """The published Retry-After exceeds this job's request budget."""


def retry_after(response):
    value = response.headers.get("Retry-After", "").strip()
    if not value:
        return None
    try:
        return max(0., float(value))
    except ValueError:
        try:
            when = parsedate_to_datetime(value)
            return max(0., (when - datetime.now(timezone.utc)).total_seconds())
        except (TypeError, ValueError, OverflowError):
            return None


def get(url, **kwargs):
    for attempt in range(3):
        # Keep one session per worker; don't rotate identities or substitute endpoints on denial.
        wait = getattr(THREAD, "next_request", 0.) - time.monotonic()
        if wait > 120:
            raise ArchiveDeferred(f"Server requested another {wait:.0f}s before any request")
        if wait > 0:
            time.sleep(wait)
        try:
            r = session().get(url, timeout=(12, 65), **kwargs)
            THREAD.next_request = time.monotonic() + 2.
            retry = retry_after(r)
            # Retry an explicit temporary refusal, including a rate-limited 403.
            # Plain Forbidden is recorded and respected.
            transient = r.status_code in (408, 425, 429, 500, 502, 503, 504)
            if r.status_code == 403 and retry is not None:
                transient = True
            if transient:
                delay = max(retry or 0., min(90, 10 * 3**attempt))
                THREAD.next_request = time.monotonic() + delay
                if attempt < 2 and delay <= 120:
                    continue
            if r.status_code == 403:
                THREAD.next_request = max(THREAD.next_request, time.monotonic() + 5.)
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout):
            THREAD.next_request = time.monotonic() + min(90, 10 * 3**attempt)
            if attempt == 2:
                raise


def norm(k):
    return re.sub("[^a-z0-9]", "", str(k).lower())


def field(row, *keys):
    bykey = {norm(k): v for k, v in row.items()}
    return next((bykey[norm(k)] for k in keys if norm(k) in bykey), None)


def number(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def json_doc(content):
    return json.loads(content.decode("utf-8-sig").strip())


class CadenceButtons(HTMLParser):
    def __init__(self):
        super().__init__()
        self.current = None
        self.options = {}

    def handle_starttag(self, tag, attrs):
        if tag == "button":
            self.current = dict(attrs).get("id", "")

    def handle_data(self, text):
        if self.current:
            self.options[self.current] = self.options.get(self.current, "") + text

    def handle_endtag(self, tag):
        if tag == "button":
            self.current = None


def index(out):
    """Select the published LCR population from the 4FGL-DR2 table used by NASA's own UI."""
    out.mkdir(parents=True, exist_ok=True)
    rawdir = out / "bootstrap"
    rawdir.mkdir(exist_ok=True)
    catalog_url = BASE + "queryDB.php?typeOfRequest=SourceList&catalog=4FGL"
    source_page = BASE + "source.html?source_name=4FGL%20J2253.9%2B1609"
    docs = {}
    contents = {}
    for name, url in (("landing", BASE), ("source_ui", source_page), ("dr2", catalog_url)):
        r = get(url)
        path = rawdir / (name + ".gz")
        save_response(path, r.content)
        docs[name] = {"url": r.url, "bytes": len(r.content), "sha256": sha(r.content), "file": str(path.relative_to(out))}
        contents[name] = r.content
    ui = contents["landing"].decode(errors="replace")
    excerpts = [line.strip()[:1600] for line in ui.splitlines() if re.search(r"21\.6|21\.7|1,?525", line, re.I)]
    save_json(out / "selection-evidence.json", {"responses": docs, "selection_excerpts": excerpts})
    # DR2 has ten annual intervals. 21.67 is the published 99% threshold on nine dof.
    threshold, inclusive = 21.67, False
    for line in excerpts:
        comparison = re.search(r"Variability_Index.{0,15}?(>=?)\s*(21\.6[0-9]+)", line, re.I)
        if comparison:
            threshold = float(comparison.group(2))
            inclusive = comparison.group(1) == ">="
            break
    items = json_doc(contents["dr2"])
    roster = []
    for row in items:
        val = number(row["Variability_Index"])
        if val is None or val < threshold or (val == threshold and not inclusive):
            continue
        n = row["Source_Name"].strip()
        if not n.startswith("4FGL J"):
            raise ValueError(f"Unexpected catalog identifier {n!r}")
        roster.append({"name": n, "assoc": (row.get("ASSOC1") or "").strip(),
                       "ra": float(row["RAJ2000"]), "dec": float(row["DEJ2000"]), "meta": row})
    roster.sort(key=lambda rec: rec["name"])
    if not 1200 <= len(roster) <= 2200 or len(set(r["name"] for r in roster)) != len(roster):
        raise ValueError(f"Unexpected DR2 variability-selected membership {len(roster)}")
    save_json(out / "catalog.json", roster)
    origin = {**docs["dr2"], "n_dr2": len(items), "n_selected": len(roster),
              "selection": f"4FGL-DR2 Variability_Index {'>=' if inclusive else '>'} {threshold}"}
    save_json(out / "catalog_origin.json", origin)
    parser = CadenceButtons()
    parser.feed(contents["source_ui"].decode(errors="replace"))
    cadences = dict(CADENCES)
    for key, label in parser.options.items():
        stripped = re.sub(r"[^a-z0-9]", "", label.lower())
        if stripped in ("3day", "3days"):
            cadences["3day"] = key
        elif stripped in ("weekly", "7days", "7day"):
            cadences["weekly"] = key
        elif stripped in ("monthly", "30days", "30day"):
            cadences["monthly"] = key
    capture = {"endpoint": BASE + "queryDB.php", "cadences": cadences,
               "source_page_cadence_controls": {k: v.strip() for k, v in parser.options.items()
                                               if k in cadences.values()}}
    save_json(out / "endpoints.json", capture)
    sample = next(x for x in roster if x["name"] == "4FGL J2253.9+1609")
    for cadence in cadences:
        alternatives = {"3day": ("daily", "3days", "3day"), "monthly": ("monthly",),
                        "weekly": ("weekly",)}
        for token in dict.fromkeys((cadences[cadence], *alternatives[cadence])):
            capture["cadences"][cadence] = token
            example = download(sample, cadence, out / "smoke", capture)
            if (example["status"] == "ok" and example["n_flux"] >= 10 and example["n_err"] >= 10
                    and example["last_mjd"] - example["first_mjd"] > 3652):
                break
            print("Archive schema diagnostic", cadence, example, flush=True)
            raw = out / "smoke" / example.get("raw", "unavailable")
            if raw.is_file():
                try:
                    obj = json_doc(gzip.decompress(raw.read_bytes()))
                    print({k: v[:1] for k, v in obj.items() if isinstance(v, list)}, flush=True)
                except (AttributeError, TypeError, ValueError):
                    pass
        else:
            print(f"LCR {cadence}: control unavailable; retain per-source acquisition coverage.", flush=True)
        save_json(out / ("example-" + cadence + ".json"), example)
    save_json(out / "endpoints.json", capture)
    print(f"Indexed {len(roster)} sources selected by {origin['selection']}", flush=True)
    return roster


def to_mjd(value):
    # Known numeric axes: MJD, JD, Fermi Mission Elapsed Time (seconds since UTC 2001-01-01),
    # Unix seconds, and JS milliseconds. MET is smoothly increasing across UTC leap seconds.
    t = number(value)
    if t is None:
        return None
    if 51000 < t < 90000:
        return t
    if 2450000 < t < 2500000:
        return t - 2400000.5
    if t > 1e11:
        return 40587.0 + t / 86400000
    if t > 1.1e9:
        return 40587.0 + t / 86400
    if 2e8 < t < 1.1e9:
        approx = 51910.0 + t / 86400
        # UTC leap seconds introduced since the MET epoch (last IERS leap 2016-12-31).
        return approx - sum(approx >= d for d in (53736, 54832, 56109, 57204, 57754)) / 86400
    return None


def series(data, *keys):
    seq = field(data, *keys)
    values = {}
    if isinstance(seq, list):
        for row in seq:
            if isinstance(row, (list, tuple)) and len(row) >= 2:
                t = to_mjd(row[0])
                vals = tuple(number(v) for v in row[1:])
                if t is not None and vals[0] is not None:
                    values[t] = vals
    return values


def parse_curve(data, cadence):
    if not isinstance(data, dict):
        raise ValueError("response is not a light curve object")
    f = series(data, "flux", "photon_flux", "flux_values")
    e = series(data, "flux_error", "flux_errors", "fluxerror", "photon_flux_error")
    u = series(data, "flux_upper_limits", "flux_upper_limit", "upper_limits", "upper_limit", "flux_ul95", "flux_ul")
    s = series(data, "ts", "test_statistic")
    c = series(data, "fit_convergence", "fit_status", "convergence", "return_code")
    times = sorted(set(f) | set(u) | set(s) | set(c))
    if not times or not (f or u):
        raise ValueError("unrecognized or empty flux/upper-limit schema; keys " + ",".join(data.keys())[:500])
    # Detections and upper limits occupy separate LCR series at the selected TS threshold.
    vals = {k: np.full(len(times), np.nan, dtype=float) for k in ("flux", "err", "ul", "ts", "fit")}
    for i, t in enumerate(times):
        ts = s.get(t, (np.nan,))[0]
        fit = c.get(t, (np.nan,))[0]
        vals["ts"][i], vals["fit"][i] = ts, fit
        if np.isfinite(fit) and fit != 0:
            continue
        val = f.get(t, (np.nan,))[0]
        upper = u.get(t, (np.nan,))[0]
        if val > 0 and (not np.isfinite(ts) or ts >= 4):
            vals["flux"][i] = val
            er = e.get(t)
            if er is not None:
                if len(er) >= 2 and er[1] is not None:
                    # The FSSC errorbar series uses [time, low, high].
                    err = (er[1] - er[0]) / 2 if er[0] <= val <= er[1] else np.nan
                else:
                    err = abs(er[0])
                if err > 0:
                    vals["err"][i] = err
        if upper > 0:
            vals["ul"][i] = upper
    arr = np.array(times, dtype=float)
    width = float(np.median(np.diff(arr))) if len(arr) > 1 else None
    if len(arr) > 10 and not ({"3day": (2.5, 4), "weekly": (6, 8), "monthly": (26, 33)}[cadence][0] <= width
                              <= {"3day": (2.5, 4), "weekly": (6, 8), "monthly": (26, 33)}[cadence][1]):
        # Very sparse arrays can jump over empty bins. Verify that those gaps lie on the requested grid.
        allowed = {"3day": (3,), "weekly": (7,), "monthly": (30, 30.4375, 28)}[cadence]
        deltas = np.diff(arr)
        if not any(np.mean(abs(deltas / k - np.rint(deltas / k)) < .035) > .85 and width >= k for k in allowed):
            raise ValueError(f"observed time spacing {width:g} d disagrees with cadence {cadence}")
    return {"t": arr, **vals}


def query_url(capture, name, cadence):
    return capture["endpoint"] + "?" + urlencode({"typeOfRequest": "lightCurveData", "source_name": name,
                "cadence": capture["cadences"][cadence], "flux_type": "photon", "index_type": "fixed",
                "ts_min": "4"}, quote_via=quote)


def save_response(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Fixed gzip timestamp keeps compressed snapshots reproducible.
    path.write_bytes(gzip.compress(data, compresslevel=5, mtime=0))


def download(source, cadence, out, capture):
    name = source["name"]
    url = query_url(capture, name, cadence)
    rec = {"name": name, "cadence": cadence, "url": url, "status": "error"}
    try:
        r = get(url)
        raw = Path("raw") / slug(name) / (cadence + ".json.gz")
        save_response(out / raw, r.content)
        rec.update(raw=str(raw), bytes=len(r.content), sha256=sha(r.content))
        data = json_doc(r.content)
        vals = parse_curve(data, cadence)
        path = Path("curves") / slug(name) / (cadence + ".npz")
        (out / path).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / path, **vals)
        rec.update(status="ok", curve=str(path), n=len(vals["t"]),
                   n_flux=int(np.isfinite(vals["flux"]).sum()), n_ul=int(np.isfinite(vals["ul"]).sum()),
                   n_err=int(np.isfinite(vals["err"]).sum()), n_ts=int(np.isfinite(vals["ts"]).sum()),
                   n_fit=int(np.isfinite(vals["fit"]).sum()),
                   n_failed_fit=int(np.sum(np.isfinite(vals["fit"]) & (vals["fit"] != 0))),
                   first_mjd=vals["t"][0], last_mjd=vals["t"][-1])
    except requests.HTTPError as exc:
        response = exc.response
        rec.update(status="http_error", http_status=response.status_code, error=str(exc)[:500],
                   response_sha256=sha(response.content), response_bytes=len(response.content),
                   retry_after=response.headers.get("Retry-After"))
    except ArchiveDeferred as exc:
        rec.update(status="deferred", error=str(exc))
    except (requests.RequestException, ValueError, TypeError, KeyError, UnicodeError) as exc:
        rec["error"] = str(exc)[:500]
    return rec


def control(index_path, out):
    """Check the monthly product without making an accessible-source decision for the whole roster."""
    out.mkdir(parents=True, exist_ok=True)
    roster = json.loads(index_path.read_text())
    sample = next(x for x in roster if x["name"] == "4FGL J2253.9+1609")
    capture = json.loads((index_path.parent / "endpoints.json").read_text())
    rec = download(sample, "monthly", out, capture)
    save_json(out / "monthly-control.json", rec)
    print(f"Monthly control for 3C 454.3: {rec['status']}, valid flux bins {rec.get('n_flux', 0)}", flush=True)



def fetch(index_path, out, shard, shards, workers):
    if not 0 <= shard < shards:
        raise ValueError("invalid shard")
    roster = json.loads(index_path.read_text())[shard::shards]
    capture = json.loads((index_path.parent / "endpoints.json").read_text())
    (out / "curves").mkdir(parents=True, exist_ok=True)
    (out / "raw").mkdir(parents=True, exist_ok=True)
    jobs = [(src, cadence) for src in roster for cadence in ("monthly", "weekly", "3day")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool, (out / "manifest.jsonl").open("w") as manifest:
        futures = [pool.submit(download, src, cadence, out, capture) for src, cadence in jobs]
        done = ok = 0
        for f in concurrent.futures.as_completed(futures):
            rec = f.result()
            manifest.write(json.dumps(rec, allow_nan=False) + "\n")
            manifest.flush()
            done += 1
            ok += rec["status"] == "ok"
            if done % 50 == 0 or done == len(futures):
                print(f"shard {shard}: {done}/{len(futures)} responses, {ok} validated", flush=True)
    # Return all acquisition failures as measured coverage information. Caller checks completeness.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    ix = sub.add_parser("index")
    ix.add_argument("--out", required=True, type=Path)
    diag = sub.add_parser("control")
    diag.add_argument("--index", required=True, type=Path)
    diag.add_argument("--out", required=True, type=Path)
    ft = sub.add_parser("fetch")
    ft.add_argument("--index", required=True, type=Path)
    ft.add_argument("--out", required=True, type=Path)
    ft.add_argument("--shard", required=True, type=int)
    ft.add_argument("--shards", required=True, type=int)
    ft.add_argument("--workers", default=2, type=int)
    args = parser.parse_args()
    if args.cmd == "index":
        index(args.out)
    elif args.cmd == "control":
        control(args.index, args.out)
    else:
        fetch(args.index, args.out, args.shard, args.shards, args.workers)


if __name__ == "__main__":
    main()

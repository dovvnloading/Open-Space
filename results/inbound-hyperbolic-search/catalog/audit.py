"""Fetch or replay the complete full-precision JPL SBDB e >= 1 selection."""
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASE = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
FIELDS = (
    "spkid,pdes,full_name,kind,class,orbit_id,epoch,e,q,i,om,w,tp,sigma_e,sigma_q,"
    "source,soln_date,data_arc,first_obs,last_obs,n_obs_used,two_body,A1,A2,A3,DT,S0"
)
URL = BASE + "?" + urlencode({
    "fields": FIELDS, "full-prec": "true",
    "sb-cdata": json.dumps({"AND": ["e|GE|1"]}, separators=(",", ":")), "sort": "spkid"
})
GM_SUN = 132712440041.93938  # km^3 / s^2, JPL DE440
AU_KM = 149597870.7
ISOS = {"1I", "2I", "3I"}


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def numeric(val):
    return None if val in ("", None) else float(val)


def process(row):
    e = float(row["e"])
    q = float(row["q"])
    if e < 1 or q <= 0:
        raise ValueError(f"invalid selection/elements for {row['spkid']}")
    speed = math.sqrt(GM_SUN / AU_KM * (e - 1) / q)
    sigma = numeric(row["sigma_e"])
    z = (e - 1) / sigma if sigma and sigma > 0 else None
    title = row["full_name"].strip()
    known = next((code for code in sorted(ISOS)
                  if title.startswith(code + "/") or row["pdes"] == code), "")
    ng = ",".join(field for field in ("A1", "A2", "A3", "DT", "S0")
                  if numeric(row[field]) not in (0, None))
    return {**row, "full_name": title, "nominal_osculating_v_km_s": speed,
            "formal_heliocentric_z_e": z, "nonzero_model_fields": ng, "known_iso": known}


def summarize(out, raw):
    snapshot = json.loads(raw)
    if snapshot["signature"]["version"] != "1.0":
        raise ValueError("unexpected JPL query API version")
    if snapshot["fields"] != FIELDS.split(","):
        raise ValueError("field names/order changed")
    rows = [process(dict(zip(snapshot["fields"], r, strict=True))) for r in snapshot["data"]]
    if len(rows) != int(snapshot["count"]) or len({r["spkid"] for r in rows}) != len(rows):
        raise ValueError("incomplete or duplicated query results")
    rows.sort(key=lambda r: (-r["nominal_osculating_v_km_s"], int(r["spkid"])))
    with (out / "catalogue.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    for r in rows[:30]:
        print("RANK", r["spkid"], r["full_name"], r["nominal_osculating_v_km_s"],
              r["formal_heliocentric_z_e"], r["data_arc"], r["n_obs_used"], r["nonzero_model_fields"])
    counts = collections.Counter
    status = {
        "total": len(rows),
        "parabolic_exact": sum(float(r["e"]) == 1 for r in rows),
        "hyperbolic_nominal": sum(float(r["e"]) > 1 for r in rows),
        "formal_sigma_e_positive": sum(r["formal_heliocentric_z_e"] is not None for r in rows),
        "formal_z_e_ge_5": sum((r["formal_heliocentric_z_e"] or 0) >= 5 for r in rows),
        "has_nonzero_model_field": sum(bool(r["nonzero_model_fields"]) for r in rows),
        "by_class": dict(counts(r["class"] or "(missing)" for r in rows)),
        "by_source": dict(counts(r["source"] or "(missing)" for r in rows)),
        "by_kind": dict(counts(r["kind"] or "(missing)" for r in rows)),
        "solution_dates_missing": sum(not r["soln_date"] for r in rows),
        "latest_last_obs": max(r["last_obs"] or "" for r in rows),
        "latest_solution_Pacific": max(r["soln_date"] or "" for r in rows),
        "arcs_days_le_30": sum(r["data_arc"] is not None and float(r["data_arc"]) <= 30 for r in rows),
        "arc_missing": sum(r["data_arc"] is None for r in rows),
        "n_obs_missing": sum(r["n_obs_used"] is None for r in rows),
        "known_iso_controls": [r["full_name"] for r in rows if r["known_iso"]],
        "nominal_osculating_speed_ge_km_s": {
            str(v): sum(r["nominal_osculating_v_km_s"] >= v for r in rows) for v in (1, 3, 5, 10, 20)
        },
        "noncontrol_speed_ge_5_and_formal_z_ge_5": [
            r["full_name"] for r in rows if not r["known_iso"]
            and r["nominal_osculating_v_km_s"] >= 5 and (r["formal_heliocentric_z_e"] or 0) >= 5
        ],
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "output_csv_sha256": hashlib.sha256((out / "catalogue.csv").read_bytes()).hexdigest()
    }
    write_json(out / "summary.json", status)
    return status


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--fetch", action="store_true", help="fetch live snapshot, otherwise verify and replay saved bytes")
    args = p.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    raw_path = out / "sbdb.json.gz"
    receipt_path = out / "receipt.json"
    if args.fetch:
        request = Request(URL, headers={"User-Agent": "Open-Space-inbound-audit/1.0 (+https://github.com/dovvnloading/Open-Space)"})
        with urlopen(request, timeout=120) as resp:
            raw = resp.read()
            receipt = {"url": resp.geturl(), "http_status": resp.status,
                       "http_date": resp.headers.get("Date"), "accessed_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
        receipt.update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        raw_path.write_bytes(gzip.compress(raw, mtime=0))
        write_json(receipt_path, receipt)
    else:
        raw = gzip.decompress(raw_path.read_bytes())
        receipt = json.loads(receipt_path.read_text())
        if receipt["sha256"] != hashlib.sha256(raw).hexdigest() or receipt["bytes"] != len(raw):
            raise ValueError("snapshot checksum or byte size mismatch")
    print(json.dumps(summarize(out, raw), indent=2))


if __name__ == "__main__":
    main()

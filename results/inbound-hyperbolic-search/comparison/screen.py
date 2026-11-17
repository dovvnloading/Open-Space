"""Freeze the rolling confirmation-page catalogue and high-e bound-orbit comparison."""
import base64
import csv
import gzip
import hashlib
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "evidence")
OUT.mkdir(exist_ok=True, parents=True)
RECEIPTS = {}


def get(name, url):
    request = Request(url, headers={"User-Agent": "Open-Space-public-astrometry-research/1"})
    with urlopen(request, timeout=120) as response:
        raw = response.read()
        RECEIPTS[name] = dict(request=url, url=response.url, status=response.status,
                              http_date=response.headers.get("Date"),
                              last_modified=response.headers.get("Last-Modified"),
                              retrieved_utc=datetime.now(timezone.utc).isoformat(),
                              sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw))
    (OUT / (name + ".gz")).write_bytes(gzip.compress(raw, mtime=0))
    return raw


def write_rows(name, header, rows):
    with (OUT/name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rolling():
    url = "https://get-cp-orbits.jostjahn.de/CUSTOMOLD.ORB"
    lines = get("neocp_rolling.orb", url).decode("latin1").splitlines()
    rows = []
    for line in lines:
        if not line.strip() or line.lstrip().startswith("Object "):
            continue
        rows.append({"trksub": line[:7].strip(), "epoch": line[20:25].strip(),
                     "e": float(line[70:79]), "line": line.rstrip()})
    if not rows:
        raise RuntimeError("empty rolling catalogue")
    write_rows("neocp_screen.csv", ["trksub","epoch","e","line"], rows)
    return {"rows": len(rows), "unique_ids": len({r["trksub"] for r in rows}),
            "e_ge_1": [r["trksub"] for r in rows if r["e"] >= 1],
            "largest_e": sorted(rows, key=lambda r: -r["e"])[:3],
            "source": "Get-CP nominal solutions for current and recently removed NEOCP/PCCP objects; rolling archive, finite time window"}


def short_arc():
    fields = "spkid,pdes,full_name,e,sigma_e,q,data_arc,first_obs,last_obs,n_obs_used,source,orbit_id"
    where = {"AND": ["e|GE|0.95","e|LT|1"]}
    url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urlencode(
        {"sb-cdata": json.dumps(where, separators=(",",":")), "fields": fields, "full-prec":"true"})
    data = json.loads(get("sbdb_near_bound.json", url))
    assert fields.split(",") == data["fields"] and int(data["count"]) == len(data["data"])
    all_rows = [dict(zip(data["fields"], row)) for row in data["data"]]
    def selected(r):
        return r["source"] in {"MPC:mpo","MPC:mp1"} or (r["data_arc"] is not None and float(r["data_arc"]) <= 30)
    rows = [r for r in all_rows if selected(r)]
    z_rows = []
    for r in rows:
        r["z_to_e1"] = None
        if r["sigma_e"] is not None and float(r["sigma_e"]) > 0:
            r["z_to_e1"] = (1 - float(r["e"])) / float(r["sigma_e"])
            if r["z_to_e1"] <= 3:
                z_rows.append(r)
    write_rows("short_arc_screen.csv", data["fields"]+["z_to_e1"], rows)
    return {"all_near_bound": len(all_rows), "single_opposition_source_or_arc_le_30_days": len(rows),
            "with_e_uncertainty": sum(r["z_to_e1"] is not None for r in rows),
            "within_3_formal_sigma_of_parabola": len(z_rows),
            "recent_within_3_sigma_last_obs_ge_2025": [r["pdes"] for r in z_rows if (r["last_obs"] or "") >= "2025-01-01"],
            "recent_selected_last_obs_ge_2025": sum((r["last_obs"] or "") >= "2025-01-01" for r in rows)}


results = {}
for name, run in [("rolling", rolling), ("short_arc", short_arc)]:
    try:
        results[name] = run()
    except Exception as exc:
        results[name] = {"error": str(exc)}
for name, value in [("screen_summary.json", results), ("receipts.json", RECEIPTS)]:
    (OUT/name).write_text(json.dumps(value, indent=2, sort_keys=True)+"\n")
print(json.dumps(results))
for path in sorted(OUT.iterdir()):
    if path.suffix == ".gz" and path.stat().st_size > 200000:
        continue
    encoded = base64.b64encode(path.read_bytes()).decode()
    print("RESULT", path.name, encoded)
if any("error" in r for r in results.values()):
    raise SystemExit(1)

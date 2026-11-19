#!/usr/bin/env python3
"""Reproducible CDS COSMOS2025 LePHARE high-z *triage*; NOT a joint posterior.

Downloads and SHA256-hashes the entire CDS-fixed-width table and its ReadMe,
parses field offsets from the official byte-by-byte layout, verifies 784016
physical rows, and reports every z>6 source with a potentially large mass.
No per-source SED, photo-z PDF, cosmic-baryon exceedance or LCDM inference
is asserted by this catalog-only script.

Usage: python results/cosmos2025-baryon-bound/screen.py --out results/cosmos2025-baryon-bound/run
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import pathlib
import re
import traceback
import urllib.request

BASES = (
    "https://cdsarc.cds.unistra.fr/ftp/J/A+A/704/A339/",
    "https://cdsarc.cds.unistra.fr/ftp/J/A%2BA/704/A339/",
    "https://cdsarc.u-strasbg.fr/ftp/J/A+A/704/A339/",
)
RELEASE = "CDS J/A+A/704/A339, VizieR inserted 2026-01-06, updated 2026-03-02; not verified against IAP v1.1 (2026-05-04)"
N_EXPECTED = 784016
COLUMNS = ["zfinal", "type", "zpdf_med", "zpdf_l68", "zpdf_u68", "mass_med", "mass_l68", "mass_u68",
           "mass_minchi2", "chi2_best", "chi_star", "chi2_agn", "flag_chandra"]
BINS = [(6., 7.), (7., 8.), (8., 9.), (9., 10.), (10., 12.), (12., 15.)]

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def get(urls, dest):
    failures = []
    for url in urls:
        try:
            h = hashlib.sha256()
            size = 0
            req = urllib.request.Request(url, headers={"User-Agent": "Open-Space-COSMOS2025-reproducibility/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as f:
                while True:
                    block = r.read(1024 * 1024)
                    if not block:
                        break
                    f.write(block)
                    h.update(block)
                    size += len(block)
            if size < 100:
                raise RuntimeError("too-short response")
            return dict(url=url, downloaded_utc=utc(), size_bytes=size, sha256=h.hexdigest())
        except Exception as exc:
            failures.append(dict(url=url, error=repr(exc)))
            dest.unlink(missing_ok=True)
    raise RuntimeError("All CDS download endpoints failed: " + json.dumps(failures))

def fields_from_readme(readme):
    anchor = re.search(r"Byte-by-byte Description of file:\s*lephare\.dat", readme, re.I)
    if not anchor:
        raise ValueError("No lephare.dat layout in ReadMe (wrong version?)")
    rest = readme[anchor.end():]
    end = re.search(r"Byte-by-byte Description of file:", rest, re.I)
    if end:
        rest = rest[:end.start()]
    # Units and labels are single tokens. Parse byte ranges as 1-based inclusive.
    definition = re.compile(r"^\s*(\d+)(?:-\s*(\d+))?\s+([AIFED]\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)\s*(.*)$")
    found = []
    for line in rest.splitlines():
        m = definition.match(line)
        if m:
            start = int(m.group(1)); endbyte = int(m.group(2) or m.group(1))
            if endbyte > 388:
                continue
            found.append(dict(start=start - 1, end=endbyte, fmt=m.group(3), unit=m.group(4),
                              label=m.group(5), description=m.group(6)))
        elif found and line.strip():
            found[-1]["description"] += " " + line.strip()
    if len(found) < 38:
        raise ValueError("Unexpectedly few lephare fields: %d; first=%s" % (len(found), found[:4]))
    mapping = {}
    for col in COLUMNS:
        aliases = [col, col.replace("_", "-")]
        matches = [f for f in found if any(re.search(r"\(" + re.escape(a) + r"\)", f["description"], re.I)
                                         for a in aliases)]
        if not matches:
            matches = [f for f in found if f["label"].lower().replace("-", "_") == col]
        if col == "type" and not matches:
            matches = [f for f in found if f["label"].lower() == "type"]
        if len(matches) == 1:
            mapping[col] = matches[0]
    mandatory = ["zfinal", "type", "zpdf_u68", "mass_med", "mass_l68", "mass_u68"]
    missing = [x for x in mandatory if x not in mapping]
    if missing:
        raise ValueError("Could not uniquely map essential %s; 43-field layout: %s" %
                         (missing, [(x["label"], x["start"], x["description"][:120]) for x in found]))
    return found, mapping

def parse_field(line, field):
    raw = line[field["start"]:field["end"]].strip()
    if not raw or raw in ("--", "-", "null", "NaN"):
        return None
    try:
        num = float(raw.replace("D", "E"))
        return num if math.isfinite(num) and num not in (-99., -99.9, -999., -999.9, -9999.) else None
    except ValueError:
        return None

def triage(readme, source, out):
    full, cols = fields_from_readme(readme)
    (out / "columns.json").write_text(json.dumps({"all_fields": full, "mapped": cols}, indent=2) + "\n")
    n = 0
    z6 = 0
    ztail6 = 0
    z_gt15 = 0
    counts = {}
    rows = []
    bad = 0
    openfile = gzip.open if source.name.endswith(".gz") else open
    with openfile(source, "rt", encoding="ascii", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            n += 1
            vals = {k: parse_field(line, fld) for k, fld in cols.items()}
            z, zu = vals.get("zfinal"), vals.get("zpdf_u68")
            mass, mu = vals.get("mass_med"), vals.get("mass_u68")
            typ = vals.get("type")
            if z is not None and z > 6 and typ == 0:
                z6 += 1
            if zu is not None and zu > 6 and typ == 0:
                ztail6 += 1
            if z is not None and z > 15 and typ == 0:
                z_gt15 += 1
            if z is None or mass is None or typ is None:
                bad += 1
            if typ != 0 or (z is None or z <= 6) and (zu is None or zu <= 6):
                continue
            item = {"source_row_1based": n}
            item.update(vals)
            rows.append(item)
            for lo, hi in BINS:
                if z is not None and lo < z <= hi:
                    key = "%g<z<=%g" % (lo, hi)
                    binitem = counts.setdefault(key, {"n_zfinal": 0, "n_mass_l68_over_10": 0,
                                                      "max_mass_med": None, "max_mass_u68": None})
                    binitem["n_zfinal"] += 1
                    if vals.get("mass_l68") is not None and vals["mass_l68"] >= 10:
                        binitem["n_mass_l68_over_10"] += 1
                    for name, attr in (("max_mass_med", "mass_med"), ("max_mass_u68", "mass_u68")):
                        value = vals.get(attr)
                        if value is not None and (binitem[name] is None or value > binitem[name]):
                            binitem[name] = value
    if n != N_EXPECTED:
        raise RuntimeError("Row count %d != documented 784016; refusing partial scan" % n)
    # A redshift upper 68% interval is NOT P(z>6). Preserve all qualifying
    # objects for later full posterior matching, then cap compact repo output.
    rows.sort(key=lambda r: r.get("mass_u68") if r.get("mass_u68") is not None else -1e6, reverse=True)
    fields = ["source_row_1based"] + [c for c in COLUMNS if c in cols]
    # Keep ALL marginal high-z triage rows, not only the top-ranked shortlist.
    # gzip mtime=0 makes this derived evidence stable across reruns.
    import io
    with (out / "highz_triage_all.csv.gz").open("wb") as fileobj:
        with gzip.GzipFile(filename="",fileobj=fileobj,mode="wb",mtime=0) as compressed:
            with io.TextIOWrapper(compressed,encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=fields,extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
    digest=hashlib.sha256((out / "highz_triage_all.csv.gz").read_bytes()).hexdigest()
    with (out / "highz_triage_top250.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows[:250])
    return {"row_count": n, "galaxy_zfinal_gt6": z6,
            "galaxy_zpdf_u68_gt6": ztail6, "galaxy_zfinal_gt15": z_gt15, "galaxy_zfinal_or_u68_gt6": len(rows),
            "unusable_rows": bad, "bin_summary": counts,
            "triage_all": "highz_triage_all.csv.gz",
            "triage_all_sha256": digest,
            "triage_top250": "highz_triage_top250.csv",
            "selection_warning": "Univariate zPDF 68% upper bound only: NOT a full joint redshift-mass posterior, NOT a test of exceedance."}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    report = {"started_utc": utc(), "release": RELEASE, "expected_rows": N_EXPECTED, "status": "started",
              "selection": "Type=0 AND (zfinal>6 OR zPDF-u68>6), ranked by marginal mass-u68"}
    try:
        rpath = out / "ReadMe"
        cpath = out / "lephare.dat.gz"
        receipt1 = get([b + "ReadMe" for b in BASES], rpath)
        report["readme"] = receipt1
        try:
            receipt2 = get([b + "lephare.dat.gz" for b in BASES], cpath)
        except Exception:
            cpath = out / "lephare.dat"
            receipt2 = get([b + "lephare.dat" for b in BASES], cpath)
        report["lephare"] = receipt2
        report["catalogue"] = triage(rpath.read_text(encoding="latin-1"), cpath, out)
        report["status"] = "catalogue-triage-complete; astrophysical checks pending"
    except Exception:
        report["status"] = "FAILED (incomplete data; no scientific conclusion)"
        report["error"] = traceback.format_exc()[-3500:]
    report["ended_utc"] = utc()
    (out / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    # Raw catalogue and ReadMe are intentionally not retained in the git tree;
    # compact JSON provenance and small CSV are the public durable receipts.
    if report["status"].startswith("FAILED"):
        raise SystemExit(2)

if __name__ == "__main__":
    main()

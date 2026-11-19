#!/usr/bin/env python3
"""Reconstruct the fixed Auger 2024 pilot from the official Zenodo summary.zip.

Requires only the Python 3 standard library and working outbound HTTPS.
Does not use GitHub Actions. Stops on byte-level mirror mismatch by default.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

DOI = "10.5281/zenodo.10488964"
URL = "https://zenodo.org/records/10488964/files/summary.zip?download=1"
EXPECTED_ZIP_MD5 = "8f882fb04dcdc3d267b72ea62c78c49c"  # Zenodo metadata, verify on use
MIRROR_GIT_SHA1 = {
    "dataSummarySD1500.csv": "5b085d1dd826e4bb0bd0439b4ac93578b05f349a",
    "dataSummaryInclined.csv": "71fc0b9f6bc6dd3c93ddb4e6872b6bccffc58bd8",
    "dataSummarySD750.csv": "a45be721265ded02f95b4d222bfc339f5f1ea406",
}

def git_blob_oid(blob):
    return hashlib.sha1(("blob %s\0" % len(blob)).encode() + blob).hexdigest()

def number(value):
    x = float(value)
    import math
    if not math.isfinite(x):
        raise ValueError("Nonfinite")
    return x

def canonical(x):
    return format(x, ".15g")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, help="Use a predownloaded official ZIP instead of HTTPS")
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "inputs" / "auger_mirror_sd_ge8.csv")
    parser.add_argument("--allow-mirror-mismatch", action="store_true", help="Permit new official bytes after careful provenance review")
    args = parser.parse_args()
    data = args.zip.read_bytes() if args.zip else urllib.request.urlopen(URL, timeout=120).read()
    md5 = hashlib.md5(data).hexdigest()
    if md5 != EXPECTED_ZIP_MD5 and not args.allow_mirror_mismatch:
        raise ValueError(f"Official zip MD5 mismatch: {md5} (expected {EXPECTED_ZIP_MD5})")
    manifest = {"doi": DOI, "zip_md5": md5, "zip_sha256": hashlib.sha256(data).hexdigest(), "samples": {}}
    events, seen = [], set()
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        for sample in ("SD1500", "Inclined", "SD750"):
            name = "dataSummary" + sample + ".csv"
            candidates = [p for p in archive.namelist() if p.endswith("/" + name) or p == name]
            if len(candidates) != 1:
                raise ValueError(f"Expected exactly one member {name}: {candidates}")
            raw = archive.read(candidates[0])
            sha = git_blob_oid(raw)
            if sha != MIRROR_GIT_SHA1[name] and not args.allow_mirror_mismatch:
                raise ValueError(f"Original {name} mismatches Git mirror blob: {sha}")
            rowcount = selected = duplicates = 0
            for row in csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))):
                rowcount += 1
                try:
                    e = number(row["sd_energy"])
                    ra, dec, zen = [number(row[k]) for k in ("sd_ra", "sd_dec", "sd_theta")]
                    if not (e >= 8 and 0 <= ra < 360 and -90 <= dec <= 90 and 0 <= zen <= 80):
                        continue
                    selected += 1
                    if sample == "SD750":
                        continue
                    sid = row["sdid"]
                    if sid in seen:
                        duplicates += 1
                        continue
                    seen.add(sid)
                    events.append((sample, sid, row["gpstime"], canonical(e), canonical(ra), canonical(dec), canonical(zen)))
                except (ValueError, KeyError, TypeError):
                    continue
            manifest["samples"][sample] = {"member": candidates[0], "git_blob_sha1": sha, "sha256": hashlib.sha256(raw).hexdigest(), "rows": rowcount, "above8": selected, "duplicate_above8": duplicates}
    events.sort(key=lambda v: v[1])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("sample","sdid","gpstime","energy_eev","ra_deg","dec_deg","zenith_deg"))
        writer.writerows(events)
    manifest["selected_ledger"] = {"rows":len(events),"sha256":hashlib.sha256(args.out.read_bytes()).hexdigest(),"path":str(args.out)}
    print(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()

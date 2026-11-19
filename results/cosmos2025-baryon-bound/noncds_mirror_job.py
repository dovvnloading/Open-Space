#!/usr/bin/env python3
"""COSMOS2025 non-CDS official-source discovery, frozen FITS receipt and conservative triage.

This intentionally NEVER submits publisher CAPTCHA/forms, invents an IRSA/Zenodo
mirror, or treats fixed-z LePHARE mass intervals as a joint p(M*,z). All downloaded
data live outside the git tree. A full object-wise significance claim is not made.
"""
from __future__ import annotations
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import socket
import sys
import traceback
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urljoin, urlsplit
from urllib.request import Request, urlopen

MAX_HTML = 2_000_000
MAX_SOURCE_BYTES = 7_000_000_000
EXPECTED_ROWS = 784_016
UA = "Open-Space-COSMOS2025-nonCDS-release-audit/2.0 (public astronomy data)"
PAGES = [
    ("iap_catalogue", "https://cosmos2025.iap.fr/catalog.html"),
    ("iap_download_form", "https://cosmos2025.iap.fr/catalog_download.php"),
    ("caltech_dr1", "https://cosmos.astro.caltech.edu/page/cosmosweb-dr"),
    ("irsa_cosmos", "https://irsa.ipac.caltech.edu/data/COSMOS/overview.html"),
    ("irsa_tap_tables", "https://irsa.ipac.caltech.edu/TAP/tables"),
    ("zenodo_cosmos2025", "https://zenodo.org/api/records?q=COSMOS2025&size=25"),
    ("zenodo_mastercatalog", "https://zenodo.org/api/records?q=COSMOSWeb_mastercatalog&size=25"),
    ("publisher_paper", "https://www.aanda.org/articles/aa/full_html/2025/12/aa55799-25/aa55799-25.html"),
    ("github_control", "https://github.com/"),
]
ALLOWED_PUBLISHER_HOSTS = {"cosmos2025.iap.fr", "cosmos.astro.caltech.edu"}
EXPECTED_PHOT = "PHOTOMETRY"
EXPECTED_PHYS = "LEPHARE"


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def request(url, timeout=40):
    return urlopen(Request(url, headers={"User-Agent": UA, "Accept": "*/*"}), timeout=timeout)


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.has_form = False
        self.has_captcha = False
    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if tag == "form":
            self.has_form = True
        for field in ("href", "data-url"):
            if d.get(field):
                self.links.append(d[field])
        if "captcha" in " ".join(str(v) for v in d.values()).lower():
            self.has_captcha = True


def probe(entry):
    label, url = entry
    out = {"id": label, "url": url, "started_utc": now()}
    try:
        host = urlsplit(url).hostname
        out["host"] = host
        addr = socket.getaddrinfo(host, 443, type=socket.SOCK_STREAM)
        out["dns"] = {"status": "ok", "addresses": sorted({a[4][0] for a in addr})[:8]}
    except Exception as exc:
        out["dns"] = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc)}
        out["http"] = {"status": "not_attempted_no_dns"}
        return out, b""
    try:
        with request(url) as r:
            # Bounded metadata sample; no credentials, cookies or form submission.
            sample = r.read(MAX_HTML + 1)
            out["http"] = {"status": r.status, "final_url": r.url,
                           "content_type": r.headers.get("Content-Type"),
                           "content_length": r.headers.get("Content-Length"),
                           "sample_bytes": len(sample), "sample_sha256": hashlib.sha256(sample).hexdigest(),
                           "sample_truncated": len(sample) > MAX_HTML}
        return out, sample[:MAX_HTML]
    except Exception as exc:
        out["http"] = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc)}
        if isinstance(exc, HTTPError):
            out["http"]["code"] = exc.code
        return out, b""


def fits_link(url):
    return bool(re.search(r"\.(?:fits|fit|fits\.gz|fit\.gz)(?:\?.*)?$", url, re.I))


def publisher_links(page_records):
    links = []
    for rec, data in page_records:
        if not data or not rec["id"] in ("iap_catalogue", "caltech_dr1", "publisher_paper"):
            continue
        parser = LinkParser()
        parser.feed(data.decode("utf-8", "replace"))
        rec["page_metadata"] = {"form_present": parser.has_form,
                                "captcha_field_present": parser.has_captcha,
                                "links_seen": len(parser.links)}
        for href in parser.links:
            target = urljoin(rec["url"], href)
            if urlsplit(target).scheme != "https" or not fits_link(target):
                continue
            # A link directly on the official IAP/Caltech release is a provenance
            # pointer, including if its downloadable file is hosted elsewhere.
            if rec["id"] in ("iap_catalogue", "caltech_dr1"):
                links.append((target, "direct_href_from_" + rec["id"]))
    return links


def zenodo_links(page_records):
    links = []
    for rec, data in page_records:
        if not rec["id"].startswith("zenodo_") or not data:
            continue
        try:
            obj = json.loads(data)
            hits = obj.get("hits", {}).get("hits", [])
            rec["zenodo_records_examined"] = len(hits)
            rec["zenodo_rejected"] = []
            for record in hits:
                meta = record.get("metadata") or {}
                title = str(meta.get("title", ""))
                names = " ".join(str(x.get("name", "")) for x in meta.get("creators", []))
                rid = record.get("id")
                if "cosmos2025" not in title.lower().replace("-", "") and "cosmos-web" not in title.lower():
                    rec["zenodo_rejected"].append({"id": rid, "reason": "not_catalogue_title"})
                    continue
                if not any(x in names.lower() for x in ("shuntov", "akins", "casey", "kartaltepe", "paquereau")):
                    rec["zenodo_rejected"].append({"id": rid, "reason": "not_team_authored"})
                    continue
                for file in record.get("files", []):
                    key = file.get("key", "")
                    direct = (file.get("links") or {}).get("self") or ""
                    # No SMF summary files, thumbnail FITS, morphology-only files,
                    # or unrelated COSMOS generations accepted as the primary source.
                    if "mastercatalog" in key.lower() and fits_link(key) and direct.startswith("https://"):
                        links.append((direct, "team_Zenodo_record_" + str(rid) + "/" + key))
                    else:
                        rec["zenodo_rejected"].append({"id": rid, "file": key, "reason": "not_full_mastercatalog"})
        except Exception as exc:
            rec["zenodo_parse_failure"] = repr(exc)
    return links


def download_verified_candidate(url, source, tmpdir, log):
    entry = {"url": url, "provenance": source, "started_utc": now()}
    dest = tmpdir / "COSMOS2025_source.download"
    fitsfile = tmpdir / "COSMOS2025_source.fits"
    dest.unlink(missing_ok=True)
    fitsfile.unlink(missing_ok=True)
    try:
        count = 0
        sha = hashlib.sha256()
        with request(url, timeout=70) as response, dest.open("wb") as target:
            entry["http_status"] = response.status
            entry["final_url"] = response.geturl()
            entry["content_type"] = response.headers.get("Content-Type")
            length = response.headers.get("Content-Length")
            entry["content_length_header"] = length
            if length and int(length) > MAX_SOURCE_BYTES:
                raise ValueError("source exceeds 7GB safety cap")
            while True:
                b = response.read(1024 * 1024)
                if not b:
                    break
                count += len(b)
                if count > MAX_SOURCE_BYTES:
                    raise ValueError("source exceeded 7GB cap while downloading")
                target.write(b)
                sha.update(b)
        entry["raw_bytes"] = count
        entry["raw_sha256"] = sha.hexdigest()
        if count < 100_000:
            raise ValueError("too small to be full 784016-row FITS catalogue")
        if dest.open("rb").read(2) == b"\x1f\x8b":
            with gzip.open(dest, "rb") as src, fitsfile.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
        else:
            shutil.copyfile(dest, fitsfile)
        entry["fits_bytes"] = fitsfile.stat().st_size
        entry["fits_sha256"] = digest(fitsfile)
        from astropy.io import fits
        with fits.open(fitsfile, memmap=True, lazy_load_hdus=True) as hdus:
            h = [{"index": i, "name": ext.name,
                  "n_rows": int(len(ext.data)) if isinstance(ext, (fits.BinTableHDU, fits.TableHDU)) else None,
                  "n_cols": len(ext.columns) if isinstance(ext, (fits.BinTableHDU, fits.TableHDU)) else None}
                 for i, ext in enumerate(hdus)]
            entry["fits_extensions"] = h
            valid = (len(h) >= 3 and h[1]["n_rows"] == EXPECTED_ROWS and
                     h[2]["n_rows"] == EXPECTED_ROWS and
                     EXPECTED_PHOT in h[1]["name"].upper() and
                     EXPECTED_PHYS in h[2]["name"].upper() and
                     h[1]["n_cols"] >= 200 and h[2]["n_cols"] >= 35)
            if not valid:
                raise ValueError("FITS is not authoritative 784016-row PHOTOMETRY + LEPHARE multi-extension release")
        entry["status"] = "validated_full_cosmos2025_catalogue"
        entry["finished_utc"] = now()
        log.append(entry)
        return fitsfile, entry
    except Exception as exc:
        entry["status"] = "failed_validation_or_transfer"
        entry["error_type"] = type(exc).__name__
        entry["error"] = str(exc)
        entry["finished_utc"] = now()
        log.append(entry)
        dest.unlink(missing_ok=True)
        fitsfile.unlink(missing_ok=True)
        return None, None


def fieldmap(columns):
    names = {re.sub("[^a-z0-9]", "", n.lower()): n for n in columns}
    wanted = ["zfinal", "zpdfu68", "massmed", "massl68", "massu68", "type",
              "chi2agn", "chistar", "flagchandra"]
    return {key: names[key] for key in wanted if key in names}


def triage_fits(path, out):
    import numpy as np
    from astropy.io import fits
    with fits.open(path, memmap=True) as h:
        p, z = h[1], h[2]
        fields = fieldmap(z.columns.names)
        if not all(f in fields for f in ("zfinal", "zpdfu68", "massmed", "massu68", "type")):
            raise ValueError("Unknown LEPHARE FITS field schema: " + repr(z.columns.names))
        photo = fieldmap(p.columns.names)
        # The strict criterion is joint p(M*,z) from an *independent* SED run;
        # this is only deliberately broad screening to prepare source vetting.
        cols = [fields[k] for k in fields]
        summary = {"source_rows": len(z.data), "physics_fields": fields,
                   "photometry_columns": list(p.columns.names),
                   "mass_note": "fixed-redshift marginal mass quantiles, not a mass-redshift posterior",
                   "selection_note": "triage: Type in {0,2} AND (zfinal>6 OR zPDF_upper68>6); incomplete for secondary high-z modes",
                   "all_marginal_highz_rows": 0, "galaxy_point_z_gt6": 0,
                   "agn_triage": 0, "bin_counts_point_z": {}, "top250": []}
        best = []
        bins = [(6, 7), (7, 8), (8, 9), (9, 10), (10, 12), (12, 15)]
        for start in range(0, len(z.data), 20000):
            end = min(start + 20000, len(z.data))
            arrays = {}
            for key in fields:
                try:
                    arrays[key] = np.array(z.data[fields[key]][start:end], dtype=float)
                except (TypeError, ValueError):
                    arrays[key] = np.full(end - start, np.nan)
            zs, zu, mass, ty = [arrays[k] for k in ("zfinal", "zpdfu68", "massu68", "type")]
            valid = np.isin(ty, [0., 2.]) & ((zs > 6) | (zu > 6))
            ix = np.where(valid)[0]
            summary["all_marginal_highz_rows"] += len(ix)
            summary["galaxy_point_z_gt6"] += int(np.sum((ty == 0) & (zs > 6)))
            summary["agn_triage"] += int(np.sum(valid & (ty == 2)))
            for low, hi in bins:
                key = f"{low}<z<={hi}"
                summary["bin_counts_point_z"][key] = summary["bin_counts_point_z"].get(key, 0) + int(np.sum(valid & (zs > low) & (zs <= hi)))
            for i in ix:
                r = {"source_row_1based": int(start + i + 1)}
                for key in fields:
                    v = arrays[key][i]
                    r[key] = float(v) if np.isfinite(v) else None
                best.append(r)
            # Retain top 250 plus a little slack after each 20k chunk, by
            # M_u68; this cut is NOT valid for a final statistical inference.
            if len(best) > 1000:
                best.sort(key=lambda r: r.get("massu68") if r.get("massu68") is not None else -1000, reverse=True)
                best = best[:250]
        best.sort(key=lambda r: r.get("massu68") if r.get("massu68") is not None else -1000, reverse=True)
        best = best[:250]
        assoc = {re.sub("[^a-z0-9]", "", n.lower()): n for n in p.columns.names}
        for r in best:
            idx = r["source_row_1based"] - 1
            for key in ("id", "radeg", "dedeg", "tile", "warnflag", "flagchandra"):
                if key in assoc:
                    val = p.data[assoc[key]][idx]
                    r["phot_" + key] = val.decode(errors="replace") if isinstance(val, bytes) else val.item() if hasattr(val, "item") else str(val)
        summary["top250"] = best
    out.mkdir(parents=True, exist_ok=True)
    (out / "fits_triage.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    return {"all_marginal_highz_rows": summary["all_marginal_highz_rows"],
            "galaxy_point_z_gt6": summary["galaxy_point_z_gt6"],
            "top_count": len(best), "triage_path": "fits_triage.json",
            "status": "preliminary_only_missing_joint_posteriors_and_independent_sed_refits"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tmp", type=Path, required=True)
    ap.add_argument("--official-url", default="", help="optional explicit publisher-provided URL (must be HTTPS FITS)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.tmp.mkdir(parents=True, exist_ok=True)
    record = {"started_utc": now(), "runner": os.getenv("RUNNER_NAME", "local"),
              "github_run_id": os.getenv("GITHUB_RUN_ID"),
              "required_release": "COSMOS2025 official IAP v1.1 2026-05-04",
              "expected_rows": EXPECTED_ROWS,
              "refuse_cds": True, "publisher_form_submission_attempted": False,
              "status": "STARTED", "source_bytes_verified": 0,
              "source_sha256": None, "checks": [], "download_attempts": [],
              "scientific_disposition": "not_measured"}
    code_hash = digest(Path(__file__))
    record["script_sha256"] = code_hash
    try:
        with ThreadPoolExecutor(max_workers=9) as pool:
            checked = list(pool.map(probe, PAGES))
        record["checks"] = [p[0] for p in checked]
        candidate_links = publisher_links(checked) + zenodo_links(checked)
        # IRSA TAP table-list is a discovery probe, not proof of a full
        # 784016-row multi-extension COSMOS2025 catalogue.
        for r, b in checked:
            if r["id"] == "irsa_tap_tables" and b:
                r["mentions_cosmos2025"] = "cosmos2025" in b.lower().decode("utf-8", "replace")
            if r["id"] == "iap_download_form" and b:
                parsed = LinkParser()
                parsed.feed(b.decode("utf-8", "replace"))
                r["requires_form"] = parsed.has_form
                r["captcha_detected"] = parsed.has_captcha or b"minicaptcha" in b.lower()
        if args.official_url:
            host = urlsplit(args.official_url).hostname
            if urlsplit(args.official_url).scheme != "https" or not fits_link(args.official_url) or (
                host not in ALLOWED_PUBLISHER_HOSTS and host != "zenodo.org"
            ):
                raise ValueError("Explicit URL must be HTTPS FITS on preapproved official publisher or Zenodo host")
            candidate_links.insert(0, (args.official_url, "operator_provided_approved_host_not_independently_authenticated"))
        seen = set()
        candidates = [(u, p) for u, p in candidate_links if not (u in seen or seen.add(u))]
        record["candidate_urls"] = [{"url": u, "reason": p} for u, p in candidates]
        for url, provenance in candidates[:8]:
            source, entry = download_verified_candidate(url, provenance, args.tmp, record["download_attempts"])
            if source is None:
                continue
            record["source_url"] = entry["final_url"]
            record["source_sha256"] = entry["raw_sha256"]
            record["source_bytes_verified"] = entry["raw_bytes"]
            record["fits_sha256"] = entry["fits_sha256"]
            try:
                record["preliminary_triage"] = triage_fits(source, args.out)
                record["status"] = "CATALOGUE_VERIFIED_BUT_NO_JOINT_POSTERIOR_OR_REFITS"
            except Exception as exc:
                record["status"] = "CATALOGUE_VERIFIED_BUT_TRIAGE_FAILED"
                record["triage_error"] = repr(exc)
            break
        if not record["source_sha256"]:
            record["status"] = "BLOCKED_NO_VERIFIED_NONCDS_CATALOGUE"
            record["blocked_reason"] = (
                "No official full-COSMOS2025 FITS link found/reachable without publisher form, "
                "or every discovered candidate failed transfer/784016-row FITS validation; "
                "see per-host HTTP and exact transfer errors. Not a scientific null.")
    except Exception:
        record["status"] = "ERROR_NO_SCIENCE"
        record["error"] = traceback.format_exc()[-4000:]
    finally:
        record["finished_utc"] = now()
        # No source FITS and no gigantic published tables committed to git.
        for file in args.tmp.glob("COSMOS2025_source.*"):
            file.unlink(missing_ok=True)
        files = [p for p in args.out.iterdir() if p.is_file() and p.name != "mirror_receipt.json"]
        record["derived_sha256"] = {p.name: digest(p) for p in sorted(files)}
        (args.out / "mirror_receipt.json").write_text(json.dumps(record, indent=2, ensure_ascii=False, default=str) + "\n")
    print(json.dumps({"status": record["status"], "source": record["source_url"] if "source_url" in record else None,
                      "sha256": record["source_sha256"], "candidate_links": len(record.get("candidate_urls", [])),
                      "host_checks": len(record["checks"])}, indent=2), flush=True)
    # This job can never assert the requested discovery without independently
    # fitted joint posterior + full source/systematic vetting.
    return 2 if record["status"] != "CATALOGUE_VERIFIED_BUT_NO_JOINT_POSTERIOR_OR_REFITS" else 0


if __name__ == "__main__":
    sys.exit(main())

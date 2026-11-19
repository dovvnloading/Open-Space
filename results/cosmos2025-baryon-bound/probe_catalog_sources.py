#!/usr/bin/env python3
"""Bounded, auditable public COSMOS2025 route probe; never bypass IAP access form.

Run locally (not GitHub Actions):
  python probe_catalog_sources.py --out host_probe.json
HTTP HEAD only after DNS works. No datasets are downloaded or claimed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import socket
import ssl
import urllib.error
import urllib.parse
import urllib.request

# "mirror discovery" probes are not assertions that those sites host COSMOS2025.
ROUTES = [
    ("IAP official catalogue page", "https://cosmos2025.iap.fr/catalog.html", "publisher"),
    ("IAP official request form", "https://cosmos2025.iap.fr/catalog_download.php", "publisher gated form"),
    ("COSMOS-Web public release portal", "https://cosmos.astro.caltech.edu/page/cosmosweb-dr", "publisher pointing to IAP"),
    ("IRSA paper bibliography", "https://irsa.ipac.caltech.edu/bibdata/2025/S/2025A%26A...704A.339S.html", "possible data-pointer; not verified catalogue mirror"),
    ("IRSA TAP capability", "https://irsa.ipac.caltech.edu/TAP/capabilities", "candidate query service; COSMOS2025 table not verified"),
    ("Zenodo records API search", "https://zenodo.org/api/records?q=COSMOS2025&size=1", "candidate mirror discovery; not verified mirror"),
    ("CDS published catalogue", "https://cdsarc.cds.unistra.fr/ftp/J/A%2BA/704/A339/ReadMe", "published older catalogue"),
    ("MAST bibliography", "https://archive.stsci.edu/mastbibref.php?bibcode=2025A%26A...704A.339S", "possible data-pointer; not verified catalogue mirror"),
    ("GitHub connectivity control", "https://github.com/", "network control only"),
]

def probe(item):
    label, url, role = item
    host = urllib.parse.urlsplit(url).hostname
    record = {"label": label, "url": url, "role": role, "hostname": host}
    try:
        results = socket.getaddrinfo(host, 443, type=socket.SOCK_STREAM)
        ips = sorted({r[4][0] for r in results})
        record["dns"] = {"status": "resolved", "addresses": ips[:6]}
    except OSError as exc:
        record["dns"] = {"status": "failed", "exception": type(exc).__name__, "errno": getattr(exc, "errno", None), "message": str(exc)}
        record["https"] = {"status": "not_attempted_dns_failed"}
        return record
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Open-Space-COSMOS2025-source-route-probe/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=4, context=ssl.create_default_context()) as response:
            record["https"] = {"status": "response", "code": response.status,
                               "content_type": response.headers.get("Content-Type"),
                               "content_length": response.headers.get("Content-Length"),
                               "final_url": response.geturl()}
    except Exception as exc:
        detail = {"status": "failed", "exception": type(exc).__name__, "message": str(exc)}
        if isinstance(exc, urllib.error.HTTPError):
            detail["http_code"] = exc.code
        record["https"] = detail
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="host_probe.json")
    args = parser.parse_args()
    with ThreadPoolExecutor(max_workers=9) as pool:
        checks = list(pool.map(probe, ROUTES))
    accessible = [x for x in checks if x["https"].get("status") == "response" and 200 <= x["https"].get("code", 0) < 400]
    receipt = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "execution_environment": "local ChatGPT execution container (not GitHub Actions)",
        "request_type": "DNS A/AAAA resolution followed by HTTPS HEAD only on success; no dataset payload fetch",
        "catalogue_download_attempted": False,
        "catalogue_bytes_received": 0,
        "catalogue_sha256": None,
        "catalogue_rows_verified": 0,
        "access_result": "HOST_CHECKS_NO_HTTP_ROUTE" if not accessible else "SOME_HTTP_ENDPOINTS_ACCESSIBLE_NOT_CATALOGUE_PROOF",
        "important_caveat": "DNS failure of the sandbox does not establish publisher outage; IAP uses a request form and is not an open direct FITS endpoint; IRSA/Zenodo not confirmed mirrors.",
        "checks": checks,
        "script_sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(receipt, f, indent=2, ensure_ascii=False)
        f.write("\n")
    for x in checks:
        print(x["hostname"], x["dns"]["status"], x["https"]["status"])
    print("RESULT:", receipt["access_result"])

if __name__ == "__main__":
    main()

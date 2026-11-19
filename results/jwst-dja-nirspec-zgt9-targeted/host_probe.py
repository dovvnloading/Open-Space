#!/usr/bin/env python3
"""Bounded direct DNS / HTTPS HEAD probes. Does not download archive FITS."""
import concurrent.futures as cf
import datetime as dt
import json
import platform
import subprocess
from urllib.parse import urlsplit

URLS = {
 "zenodo-summary": "https://zenodo.org/records/15472354/files/dja_msaexp_emission_lines_v4.4.csv.gz?download=1",
 "dja-s3-summary-path": "https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.4.csv.gz",
 "dja-s3-virtual-host": "https://msaexp-nirspec.s3.amazonaws.com/extractions/dja_msaexp_emission_lines_v4.4.csv.gz",
 "dja-documentation": "https://dawn-cph.github.io/dja/general/api_summary/",
 "mast-api": "https://mast.stsci.edu/api/v0/invoke",
 "mast-archive": "https://archive.stsci.edu/",
}

def call(argv, max_seconds):
    try:
        p = subprocess.run(argv, text=True, capture_output=True, timeout=max_seconds)
        return {"returncode": p.returncode, "stdout": p.stdout.strip()[:300],
                "stderr": p.stderr.strip()[:300]}
    except subprocess.TimeoutExpired:
        return {"returncode": 124, "stdout": "",
                "stderr": f"timeout after {max_seconds}s"}

def probe(item):
    name, url = item
    host = urlsplit(url).hostname
    dns = call(["getent", "ahostsv4", host], 4)
    https = call(["curl", "--noproxy", "*", "--head", "--location",
                  "--silent", "--show-error", "--connect-timeout", "3",
                  "--max-time", "7", "--output", "/dev/null", "--write-out",
                  "http_code=%{http_code} remote_ip=%{remote_ip} err=%{errormsg}", url], 9)
    return {"name": name, "url": url, "hostname": host,
            "dns": dns, "https_head": https}

if __name__ == "__main__":
    with cf.ThreadPoolExecutor(max_workers=len(URLS)) as executor:
        results = list(executor.map(probe, URLS.items()))
    print(json.dumps({
        "checked_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "environment": "ChatGPT compute container (NOT GitHub Actions/web crawler)",
        "os": platform.system() + " " + platform.release(),
        "probe": "getent ahostsv4 + direct curl --noproxy '*' HTTPS HEAD",
        "results": results,
        "limitation": "No raw source catalog or spectra received; no astrophysical inference."
    }, indent=2))

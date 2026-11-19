#!/usr/bin/env python3
"""Bounded HEASARC XMM PPS path/product probe; not a sky screen."""
import datetime as dt
import hashlib
import html
import json
import re
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse, unquote
import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
HOST = "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/rev0/"
S = requests.Session()
S.headers["User-Agent"] = "Open-Space-DR14-QPE-HEASARC/1.0 (research receipts)"

def listing(obsid, retries=2):
    url = HOST + obsid + "/PPS/"
    trials = []
    for k in range(retries):
        t0 = time.monotonic()
        trial = {"attempt": k+1, "url": url}
        try:
            with S.get(url, timeout=(12, 22)) as response:
                trial.update(http_status=response.status_code, content_type=response.headers.get("Content-Type"),
                             bytes=len(response.content), elapsed_s=round(time.monotonic()-t0, 2))
                response.raise_for_status()
                body = response.text[:4000000]
                hrefs = re.findall(r'href\s*=\s*["\x27]([^"\x27]+)', body, flags=re.I)
                entries=[]
                for ref in hrefs:
                    full = urljoin(url, html.unescape(ref))
                    parsed = urlparse(full)
                    if parsed.netloc.lower() != "heasarc.gsfc.nasa.gov":
                        continue
                    if not parsed.path.startswith("/FTP/xmm/data/rev0/"+obsid+"/PPS/"):
                        continue
                    name = unquote(parsed.path.rsplit("/", 1)[-1])
                    if re.fullmatch(r"[A-Za-z0-9_.+-]{4,200}", name):
                        entries.append(name)
                trial["entries"] = len(set(entries))
                trial["pps_src_tsr"] = len([x for x in entries if "PN" in x.upper() and "SRCTSR" in x.upper()])
                trial["pps_fbk_tsr"] = len([x for x in entries if "PN" in x.upper() and "FBKTSR" in x.upper()])
                trial["pn_examples"] = [x for x in entries if x.upper().startswith("P"+obsid+"PN")][:25]
                trial["html_sha256"] = hashlib.sha256(response.content).hexdigest()
                if not entries:
                    trial["short_body"] = body[:300]
                trials.append(trial)
                return trials
        except Exception as e:
            trial["error"]=repr(e)[:240]
            trial["elapsed_s"]=round(time.monotonic()-t0,2)
            trials.append(trial)
            if k+1 < retries: time.sleep(1+k)
    return trials

def main():
    cfg = json.loads((BASE/"control_batch.json").read_text())
    assert cfg["stage"] == "heasarc_probe"
    obsids=["0123700101", "0872392901", "0005010301", "0851180401"]
    rec={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "host":HOST,"purpose":"US HEASARC mirror PPS directory and source/background filename probe",
         "obsids":{},"complete_survey":False,"pn_curves_screened":0}
    for obsid in obsids:
        rec["obsids"][obsid]=listing(obsid)
        (OUT/"heasarc_probe_01.json").write_text(json.dumps(rec,indent=2)+"\n")
        print("HEASARC_PROBE",obsid,json.dumps(rec["obsids"][obsid][-1])[:1000],flush=True)

if __name__=="__main__": main()

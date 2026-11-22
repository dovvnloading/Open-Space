"""Capture the official website's source-table and light-curve request schemas."""
import base64
import hashlib
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from selenium import webdriver
from bs4 import BeautifulSoup

root = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/"
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
log = []
def save():
    (out / "schema.json").write_text(json.dumps(log, indent=2))
def capture(url, payload):
    data = payload.encode() if isinstance(payload, str) else payload
    text = data.decode(errors="replace")
    rec = {"url": url, "nbytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    try:
        j = json.loads(text)
        def head(val):
            if isinstance(val,list):
                return [head(val[0])] if val else []
            if isinstance(val,dict):
                return {k:head(v) for k,v in list(val.items())[:90]}
            return str(val)[:220]
        rec["schema"] = head(j)
        rec["rows_or_keys"] = len(j)
    except (ValueError,TypeError):
        rec["excerpts"] = [l.strip()[:1300] for l in text.splitlines() if re.search(r'queryDB|ajax|\.json|source_name|\.php|\.csv|time_interval|cadence|dataSource|source\.html|<script.*src',l,re.I)][:120]
    log.append(rec)
    save()

pages = [root, urljoin(root,"source.html?source_name=4FGL%20J1555.7%2B1111")]
seen = set()
for u in pages.copy():
    try:
        r=requests.get(u,timeout=25)
        capture(r.url,r.content)
        soup=BeautifulSoup(r.content,"html.parser")
        scripts=[urljoin(u,s.get("src","")) for s in soup.select("script[src]")]
        for src in scripts:
            if not src.startswith(root) or src in seen or re.search(r"jquery|plotly|bootstrap|aladin|datatables|\.min\.js",src,re.I):continue
            seen.add(src)
            p=requests.get(src,timeout=15)
            capture(p.url,p.content)
    except Exception as e:
        log.append({"failure":str(e)[:350],"url":u});save()
try:
    opt=webdriver.ChromeOptions()
    for arg in ("--headless=new","--no-sandbox","--disable-dev-shm-usage"):opt.add_argument(arg)
    opt.set_capability("goog:loggingPrefs",{"performance":"ALL"})
    driver=webdriver.Chrome(options=opt)
    driver.set_page_load_timeout(35)
    driver.execute_cdp_cmd("Network.enable",{})
    for url in pages:
        try:driver.get(url)
        except Exception:pass
        time.sleep(5)
        for entry in driver.get_log("performance"):
            ev=json.loads(entry["message"])["message"]
            if ev["method"]!="Network.responseReceived":continue
            par=ev["params"]; u=par["response"]["url"]
            if not u.startswith(root) or par["type"] not in ("XHR","Fetch"):continue
            try:
                b=driver.execute_cdp_cmd("Network.getResponseBody",{"requestId":par["requestId"]})
                capture(u,base64.b64decode(b["body"]) if b["base64Encoded"] else b["body"])
            except Exception as e:log.append({"url":u,"failure":str(e)[:350]})
        save()
    driver.quit()
except Exception as e:log.append({"browser_error":str(e)[:450]});save()

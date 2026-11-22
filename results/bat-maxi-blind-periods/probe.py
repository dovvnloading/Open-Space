"""Record the public archive entry points and product examples on an Actions runner."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urljoin
import hashlib
import json
import re

import requests
from bs4 import BeautifulSoup

out = Path(__file__).resolve().parent / "archive-layout.json"
bat = "https://swift.gsfc.nasa.gov/results/transients/"
maxi = "https://maxi.riken.jp/"
urls = [bat, bat + "weak/", bat + "Weak.html", bat + "Crab/", maxi, maxi + "top/slist.html",
        maxi + "star_data/J0534+220/J0534+220.html", maxi + "pubdata/v7l.rkn/J0534+220/gsc_1day_all.dat",
        "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/", "https://maxi.riken.jp/pubdata/v7l.rkn/"]
def inspect(url):
    try:
        r = requests.get(url, timeout=(8, 22), headers={"User-Agent":"Open-Space astronomical period survey"})
        record = dict(url=r.url, status=r.status_code, size=len(r.content), sha256=hashlib.sha256(r.content).hexdigest())
        if "html" in r.headers.get("Content-Type","") or "<html" in r.text[:1500].lower():
            doc = BeautifulSoup(r.content,"html.parser")
            links = [(a.get_text(" ",strip=True),urljoin(r.url,a["href"])) for a in doc.select("a[href]")]
            record["links"] = (links[:12] + [v for v in links if re.search(r"weak|catalog|source list|all sources|monitor.*source|\.fits|1day|[.]lc\b|Crab|slist|1-day|daily|pubdata|\.cat",str(v),re.I)][:55])
            record["rows"] = [" ".join(tr.stripped_strings)[:400] for tr in doc.select("tr")[:8]]
        else:
            record["sample"]=r.text.splitlines()[:10]
        return record
    except Exception as exc:
        return dict(url=url,error=repr(exc))
def capture(urls):
    with ThreadPoolExecutor(max_workers=6) as pool:
        return list(pool.map(inspect,urls))
records=capture(urls)
more = []
for k in (0,4,5):
    for label,href in records[k].get("links",[]):
        if re.search(r"Crab|weak|slist|all sources|source list",label,re.I) and href not in urls+more:
            more.append(href)
records += capture(more[:10])
out.write_text(json.dumps(records,indent=2)+"\n")

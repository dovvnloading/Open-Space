"""Inventory the release documentation and table metadata inside GitHub Actions."""
import concurrent.futures as cf
import hashlib
import json
import os
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import Request, urlopen

OUT = Path("results/xmmsl2-erass1-rass/diagnostics")
OUT.mkdir(parents=True, exist_ok=True)
ERO = "https://erosita.mpe.mpg.de/dr1/"
HEA = "https://heasarc.gsfc.nasa.gov/"
VIZ = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TAP = HEA + "xamin/vo/tap/sync"


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.words = [], []
        self.hidden = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag in ("script", "style"):
            self.hidden += 1
        if tag in ("form", "a", "iframe"):
            self.links.append(attrs.get("href", attrs.get("action", attrs.get("src", ""))))

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.hidden = max(0, self.hidden - 1)

    def handle_data(self, data):
        if self.hidden == 0 and data.strip():
            self.words.append(data.strip())


def tap(base, query):
    return base + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))


def get(item):
    label, url = item
    rec = dict(label=label, url=url)
    try:
        with urlopen(Request(url, headers={"User-Agent": "Open-Space scientific release inventory"}), timeout=25) as response:
            data = response.read(1_500_001)
            rec.update(status=response.status, final=response.url, bytes_read=len(data), sha256=hashlib.sha256(data).hexdigest())
        body = data.decode("utf8", "replace")
        if "<html" in body.lower() or "<!doctype html" in body.lower():
            p = Page()
            p.feed(body)
            rec.update(links=list(dict.fromkeys(urljoin(url, l) for l in p.links if l))[:260], text=" ".join(p.words)[:18000])
        else:
            rec["text"] = body[:18000]
    except Exception as exc:
        rec["error"] = str(exc)
    return rec


queries = {
    "ero": ERO, "erodata": ERO+"AllSkySurveyData_dr1/", "erocatalogues": ERO+"AllSkySurveyData_dr1/Catalogues_dr1/",
    "erodat": ERO+"erodat/", "eroul": ERO+"erodat/upperlimit/", "erosky": ERO+"erodat/skyview/",
    "rass_archive": HEA+"FTP/rosat/data/pspc/processed_data/rass/", "rass_root": HEA+"FTP/rosat/data/pspc/processed_data/",
    "rass_maps": HEA+"FTP/rosat/data/pspc/rass/", "rass_images": HEA+"docs/rosat/survey/rass.html",
    "rass_info": "https://www.xray.mpe.mpg.de/rosat/survey/rass-3/",
    "skyview_cnt": "https://skyview.gsfc.nasa.gov/current/cgi/survey.pl?RASS-Cnt+Broad",
    "skyview_exp": "https://skyview.gsfc.nasa.gov/current/cgi/survey.pl?RASS-Exp+Broad",
    "skyview_surveys":"https://skyview.gsfc.nasa.gov/current/surveys.html",
    "hea_tables":tap(TAP,"SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%erass%' OR table_name LIKE '%2rxs%' OR table_name LIKE '%slew%'"),
    "hea_slew":tap(TAP,"SELECT column_name,unit,description FROM TAP_SCHEMA.columns WHERE table_name='xmmslewcln'"),
    "viz_tables":tap(VIZ,"SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%682/A34%' OR table_name LIKE '%588/A103%'"),
    "ned":tap("https://ned.ipac.caltech.edu/tap/sync","SELECT TOP 1 * FROM NEDTAP.objdir"),
    "tns":"https://www.wis-tns.org/search?format=csv&ra=15&decl=0&radius=30&coords_unit=arcsec",
    "slew_doc":HEA+"W3Browse/xmm-newton/xmmslewcln.html",
    "rass_doc":HEA+"W3Browse/rosat/rass2rxs.html",
}
with cf.ThreadPoolExecutor(max_workers=7) as pool:
    rows = list(pool.map(get, queries.items()))
    seen=set(queries.values())
    for depth in range(2):
        extra = []
        for row in rows:
            for link in row.get("links", []):
                if link in seen or not link.startswith(ERO) or urlsplit(link).fragment:
                    continue
                if re.search(r"upper|limit|sensitiv|skytile|exposure|[Mm]aps|catalog|download|[Aa]ccess|[Ff]lare", link) and not re.search(r"\.fits|\.gz|\.zip|\.pdf", link, re.I):
                    seen.add(link)
                    extra.append(("ero_link",link))
        if extra:
            rows.extend(pool.map(get, extra[:30]))
for i, row in enumerate(rows):
    (OUT/f"{i:02d}_{row['label']}.json").write_text(json.dumps(row,indent=2)+"\n")
(OUT/"index.json").write_text(json.dumps({"run":os.environ.get("GITHUB_RUN_ID"), "documents":[(r["label"],r["url"],r.get("status"), r.get("error")) for r in rows]},indent=2)+"\n")

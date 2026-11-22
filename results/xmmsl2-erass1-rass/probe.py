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
    return base + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable" if base==TAP else "csv", QUERY=query))


def get(item):
    label, url = item
    rec = dict(label=label, url=url)
    try:
        with urlopen(Request(url, headers={"User-Agent": "Open-Space scientific release inventory"}), timeout=25) as response:
            data = response.read(1_500_001)
            rec.update(status=response.status, final=response.url, bytes_read=len(data), sha256=hashlib.sha256(data).hexdigest())
        body = data.decode("utf8", "replace")
        if re.search(r"<(?:html|head|body)\b",body,re.I):
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
    "ero_api": "https://erosita.mpe.mpg.de/erodat/apis/", "ero_ul_doc": ERO+"AllSkySurveyData_dr1/UpperLimitServer_dr1/",
    "xmm_archive":HEA+"FTP/xmm/data/catalogues/","archive_xmm2":"https://www.cosmos.esa.int/web/xmm-newton/xsa",
    "rass_release":HEA+"FTP/rosat/data/pspc/processed_data/rass/release/","rass_pspc":HEA+"FTP/rosat/data/pspc/",
    "rass_mpe":"https://www.mpe.mpg.de/ROSAT/2RXS/", "rass_home":HEA+"FTP/rosat/",
    "rass_metadata":"https://skyview.gsfc.nasa.gov/current/cgi/survey.pl?survey=RASS-Cnt%20Broad",
    "rass_exposure":"https://skyview.gsfc.nasa.gov/current/cgi/survey.pl?survey=RASS-Exp%20Broad",
    "xmm2_cds":"https://cdsarc.cds.unistra.fr/ftp/IX/53/", "rass_cds":"https://cdsarc.cds.unistra.fr/ftp/J/A+A/588/A103/",
    "xmm2_user":"https://www.cosmos.esa.int/web/xmm-newton/xmmsl2-ug",
    "xmm2_slew":"https://www.cosmos.esa.int/web/xmm-newton/epic-slew-catalogue",
    "nxsa":"https://nxsa.esac.esa.int/catalogues/",
    "rass_caldb":HEA+"FTP/caldb/data/rosat/pspc/cpf/matrices/",
    "rass_indices":tap(TAP,"SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%rass%' OR table_name LIKE 'rosmaster'"),
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
                if re.search(r"upper|limit|sensitiv|skytile|exposure|maps|download|access", link, re.I) and not re.search(r"\.fits|\.gz|\.zip|\.pdf", link, re.I):
                    seen.add(link)
                    extra.append(("ero_link",link))
        if depth == 0:
            release = next((r for r in rows if r["label"]=="rass_release"),{})
            samples = [l for l in release.get("links",[]) if re.search(r"/rs\d+[^/]*/$", l)]
            if samples:
                extra.insert(0, ("rass_field", samples[len(samples)//2]))
        if extra:
            rows.extend(pool.map(get, extra[:30]))
for i, row in enumerate(rows):
    (OUT/f"{i:02d}_{row['label']}.json").write_text(json.dumps(row,indent=2)+"\n")
(OUT/"index.json").write_text(json.dumps({"run":os.environ.get("GITHUB_RUN_ID"), "documents":[(r["label"],r["url"],r.get("status"), r.get("error")) for r in rows]},indent=2)+"\n")

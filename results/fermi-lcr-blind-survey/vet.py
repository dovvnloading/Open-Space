"""Collect review evidence for surviving Fermi LCR leads, from Actions only.

Records the actual queries, status codes and payload hashes, distinguishes a
failed lookup from a successful empty result, and never assigns novelty from
an empty search. A person must inspect close periods, aliases and light curves.
"""
import hashlib
import io
import json
import math
import os
import threading
import re
import shutil
import subprocess
import tempfile
import unicodedata
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote, unquote, urljoin, urlparse

import requests
from astropy.io import fits
from bs4 import BeautifulSoup


ADS = "https://api.adsabs.harvard.edu/v1/search/query"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
ATEL = "https://www.astronomerstelegram.org/"
CROSSREF = "https://api.crossref.org/works"
ARXIV = "https://export.arxiv.org/api/query"
DR4 = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/"
HEADERS = {"User-Agent": "Open-Space-Fermi-LCR/1.0 (astronomical research; github.com/dovvnloading/Open-Space)"}
TOPIC = re.compile(r"period|oscillat|supermassive.{0,20}binar|quasi.?period", re.I)


def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).lower().replace("−", "-").replace("–", "-")
    return re.sub(r"[^a-z0-9+\-]", "", s)


def dedup(seq):
    return list(dict.fromkeys(str(s).strip() for s in seq if s and str(s).strip()))


def get(url, audit, params=None, headers=None, limit=24_000_000, timeout=18):
    """Only public response metadata are archived, never an authorization header."""
    record = {"url": url, "params": params or {}}
    audit.append(record)
    try:
        with requests.get(url, params=params, headers={**HEADERS, **(headers or {})},
                          stream=True, timeout=(7, timeout)) as res:
            record.update(status=res.status_code, final_url=res.url)
            res.raise_for_status()
            data = bytearray()
            for chunk in res.iter_content(65536):
                data.extend(chunk)
                if len(data) > limit:
                    raise ValueError("response exceeded size limit")
            record.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
            return bytes(data)
    except (requests.RequestException, ValueError) as err:
        record["error"] = str(err)[:300]
        return None


def get_json(url, audit, **kwargs):
    data = get(url, audit, **kwargs)
    if data is None:
        return None
    try:
        return json.loads(data)
    except (ValueError, UnicodeError):
        audit[-1]["error"] = "response is not JSON"
        return None


def text_value(x):
    return x.decode("ascii", "replace").strip() if isinstance(x, bytes) else str(x).strip()


def dr4_catalog(audit):
    # Discover the current version in the 14-year (DR4) release only.
    page = get(DR4, audit)
    urls = re.findall(r'''href=["']([^"']*gll_psc_v\d+\.fit[s]?)["']''',
                      (page or b"").decode(errors="replace"), flags=re.I)
    files = dedup(["gll_psc_v35.fit"] + [x.rsplit("/", 1)[-1] for x in urls])
    for name in files:
        raw = get(DR4 + name, audit, limit=80_000_000)
        if not raw:
            continue
        try:
            with fits.open(io.BytesIO(raw), memmap=False) as ff:
                tab = ff["LAT_Point_Source_Catalog"].data
                keys = ["Source_Name", "ASSOC1", "ASSOC2", "ASSOC_TEV", "ASSOC_FGL",
                        "CLASS1", "Flags", "RAJ2000", "DEJ2000",
                        "Conf_95_SemiMajor", "ASSOC_PROB_BAY", "ASSOC_PROB_LR"]
                rows = []
                for row in tab:
                    values = {k: text_value(row[k]) for k in keys if k in tab.names}
                    rows.append(values)
                return rows
        except (OSError, KeyError, ValueError) as err:
            audit[-1]["error"] = "unreadable DR4 catalog: " + str(err)[:150]
    return []


def angle(ra1, dec1, ra2, dec2):
    a, d, b, e = map(math.radians, [ra1, dec1, ra2, dec2])
    cs = math.sin(d)*math.sin(e) + math.cos(d)*math.cos(e)*math.cos(a-b)
    return math.degrees(math.acos(max(-1, min(1, cs))))


def identify(lead, table):
    for row in table:
        if norm(row["Source_Name"]) == norm(lead["name"]):
            return dict(row, match="exact 4FGL name")
    if lead.get("ra") is None or lead.get("dec") is None:
        return None
    try:
        ra, dec = float(lead["ra"]), float(lead["dec"])
    except (TypeError, ValueError):
        return None
    matches = []
    for row in table:
        sep = angle(ra, dec, float(row["RAJ2000"]), float(row["DEJ2000"]))
        if sep <= float(row.get("Conf_95_SemiMajor", "nan")):
            matches.append(dict(row, sep_deg=sep, match="inside DR4 95% semimajor radius"))
    return matches[0] if len(matches) == 1 else None


def tap(query, audit):
    obj = get_json(SIMBAD, audit, params={"request": "doQuery", "lang": "adql",
                                           "format": "json", "query": query})
    if not isinstance(obj, dict) or "data" not in obj:
        return None
    cols = [m["name"] for m in obj["metadata"]]
    return [dict(zip(cols, row)) for row in obj["data"]]


def simbad(aliases, audit):
    expr = ",".join("'" + s.replace("'", "''") + "'" for s in aliases)
    ids = tap(
        "SELECT DISTINCT TOP 20 b.oid AS oid,b.main_id,b.ra,b.dec,b.otype,ids.ids "
        "FROM ident AS i JOIN basic AS b ON b.oid=i.oidref "
        "JOIN ids ON ids.oidref=b.oid WHERE i.id IN (" + expr + ")", audit)
    if not ids:
        return {"matched": ids, "photometry": None, "references": None}
    keys = ",".join(str(int(z["oid"])) for z in ids)
    phot = tap(
        "SELECT TOP 150 oidref,filter,flux,bibcode FROM flux WHERE oidref IN (" + keys + ")", audit)
    refs = tap(
        "SELECT DISTINCT TOP 600 r.bibcode,r.title,r.journal,r.year "
        "FROM has_ref AS h JOIN ref AS r ON h.oidbib=r.oidbib WHERE h.oidref IN ("
        + keys + ") ORDER BY r.year DESC", audit)
    return {"matched": ids, "photometry": phot, "references": refs}


def atel_archive(audit):
    found = {}
    for params in ({"atels": "all"}, {"all": ""}, {"atels": ""}, {}):
        page = get(ATEL, audit, params=params)
        if not page:
            continue
        html = BeautifulSoup(page, "html.parser")
        for link in html.select('a[href*="read="]'):
            num = re.search(r"[?&]read=(\d+)", link["href"])
            if not num:
                continue
            title = link.get_text(" ", strip=True)
            if len(title) < 12:
                par = link.find_parent(["tr", "li", "p"])
                if par:
                    title = par.get_text(" ", strip=True)[:500]
            if len(title) >= 12:
                found[int(num[1])] = title
        if len(found) > 10000:
            break
    return found


def arxiv_ids(title, doi, audit):
    params = {"search_query": 'ti:"' + re.sub(r'["\n]', " ", title) + '"', "max_results": 5}
    data = get(ARXIV, audit, params=params, timeout=12)
    if not data:
        return []
    try:
        root = ET.fromstring(data)
    except ET.ParseError:
        return []
    tags = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}
    result = []
    for ent in root.findall("a:entry", tags):
        paper_title = ent.findtext("a:title", "", tags)
        paper_doi = ent.findtext("x:doi", "", tags)
        if norm(title) == norm(paper_title) or (doi and doi.lower() == paper_doi.lower()):
            result.append(ent.findtext("a:id", "", tags).split("/abs/")[-1])
    return [s for s in result if re.fullmatch(r"[\d.]+(?:v\d+)?", s)]


def pdf_text(raw):
    if raw is None or b"%PDF-" not in raw[:1024] or not shutil.which("pdftotext"):
        return ""
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "article.pdf"
        source.write_bytes(raw)
        run = subprocess.run(["pdftotext", "-layout", str(source), "-"], capture_output=True,
                             timeout=30, check=False)
    if run.returncode:
        return ""
    return run.stdout[:4_000_000].decode(errors="replace")


def topic_papers(audit):
    papers = {}
    for params in [
        {"query.author": "P D Peñil", "query.bibliographic": "gamma ray periodicity blazars"},
        {"query.bibliographic": "systematic survey periodicities gamma ray Fermi LAT blazars"},
    ]:
        obj = get_json(CROSSREF, audit, params={**params, "rows": 24})
        for w in (obj or {}).get("message", {}).get("items", []):
            title = BeautifulSoup(" ".join(w.get("title", [])), "html.parser").get_text(" ", strip=True)
            if not TOPIC.search(title) or not re.search(r"blazar|fermi|gamma|active galactic", title, re.I):
                continue
            doi = w.get("DOI", "").lower()
            if not doi or doi in papers:
                continue
            author = [a.get("family", "") for a in w.get("author", [])]
            penil = any(norm(a) == "penil" for a in author)
            if not penil and sum(not p["penil"] for p in papers.values()) >= 6:
                continue
            papers[doi] = {"title": title, "doi": doi, "authors": author[:8], "penil": penil,
                           "links": [x.get("URL", "") for x in w.get("link", [])]}
    cache = []
    for doi, p in papers.items():
        # OpenAlex ties repository copies to the DOI. Limit network destinations
        # to publication archives. Where no accessible manuscript exists, keep
        # the bibliographic record and the coverage gap.
        oa = get_json("https://api.openalex.org/works/https://doi.org/" + doi, audit, timeout=10)
        p["arxiv_ids"] = []
        pdf_urls = []
        for loc in (oa or {}).get("locations", []):
            for url in [loc.get("landing_page_url") or "", loc.get("pdf_url") or ""]:
                m = re.search(r"arxiv\.org/(?:abs|pdf)/([0-9]+\.[0-9]+)", url)
                if m:
                    p["arxiv_ids"].append(m[1])
                elif url == loc.get("pdf_url"):
                    pdf_urls.append(url)
        if not p["arxiv_ids"]:
            p["arxiv_ids"] = arxiv_ids(p["title"], doi, audit)
        urls = [f"https://arxiv.org/pdf/{s}" for s in dedup(p["arxiv_ids"])] + pdf_urls + p.pop("links")
        text = ""
        for url in dedup(urls):
            host = (urlparse(url).hostname or "").lower()
            if not any(host == h or host.endswith("." + h) for h in
                       ["arxiv.org", "iop.org", "aanda.org", "oup.com", "springer.com", "mdpi.com"]):
                continue
            raw = get(url, audit, timeout=12)
            try:
                text = pdf_text(raw)
            except (subprocess.TimeoutExpired, OSError):
                text = ""
            if len(text) > 500:
                p["full_text_url"] = url
                break
        p["text_characters"] = len(text)
        cache.append((p, text))
    return cache


_UI_LOCK = threading.Lock()
_UI_DRIVER = None


def public_ads(link):
    """Read the public search UI through its usual rendered results page."""
    global _UI_DRIVER
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
    except ImportError:
        return {"status": "public search browser unavailable", "search_url": link, "records": []}
    with _UI_LOCK:
        try:
            if _UI_DRIVER is None:
                opts = webdriver.ChromeOptions()
                for flag in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage"):
                    opts.add_argument(flag)
                _UI_DRIVER = webdriver.Chrome(options=opts)
                _UI_DRIVER.set_page_load_timeout(20)
            _UI_DRIVER.get(link)
            def read_screen(driver):
                body = driver.find_element(By.TAG_NAME, "body").text
                records = {}
                for a in driver.find_elements(By.CSS_SELECTOR, 'a[href*="/abs/"]'):
                    url = a.get_attribute("href") or ""
                    match = re.search(r"/abs/([^/#?]+)", unquote(url))
                    if match and a.text.strip():
                        bib = match[1]
                        data = records.setdefault(bib, {"bibcode": bib, "title": a.text.strip()})
                        if len(a.text.strip()) > len(data["title"]):
                            data["title"] = a.text.strip()
                if records:
                    return {"status": "public ADS results page", "records": list(records.values())}
                if re.search(r"no results (?:were )?found|your search returned (?:no|0) results", body, re.I):
                    return {"status": "public ADS search returned no results", "records": []}
                return False
            result = WebDriverWait(_UI_DRIVER, 10, poll_frequency=0.6).until(read_screen)
            return dict(result, search_url=link, coverage="first rendered results page")
        except Exception as err:
            # Selenium errors vary by browser version; leave the query as unverified.
            return {"status": "public ADS search unavailable", "detail": type(err).__name__,
                    "search_url": link, "records": []}


def ads(query, audit, token):
    link = "https://ui.adsabs.harvard.edu/search/q=" + quote(query, safe="") + "&sort=date%20desc"
    if token:
        obj = get_json(ADS, audit, params={"q": query, "fl": "bibcode,title,year", "rows": 100},
                       headers={"Authorization": "Bearer " + token})
        if obj and isinstance(obj.get("response"), dict):
            resp = obj["response"]
            return {"status": "queried", "found": resp.get("numFound"), "search_url": link,
                    "records": [{"bibcode": x.get("bibcode"), "title": (x.get("title") or [""])[0],
                                 "year": x.get("year")} for x in resp.get("docs", [])]}
    return public_ads(link)


def dr4_notes(audit):
    page = get(DR4, audit)
    if not page:
        return []
    html = BeautifulSoup(page, "html.parser")
    texts = [{"source": DR4, "text": html.get_text(" ", strip=True)}]
    links = [a for a in html.select("a[href]") if re.search(
        r"note|caveat|documentation|data release 4", a.get_text(" ", strip=True), re.I)]
    urls = ["https://arxiv.org/pdf/2307.12546"]
    for a in links:
        url = urljoin(DR4, a["href"])
        if url.startswith(DR4) and len(urls) < 5:
            urls.append(url)
    for url in dedup(urls):
        raw = get(url, audit)
        if not raw:
            continue
        if b"%PDF-" in raw[:1024]:
            try:
                content = pdf_text(raw)
            except (OSError, subprocess.TimeoutExpired):
                content = ""
        elif b"<html" in raw[:3000].lower():
            content = BeautifulSoup(raw, "html.parser").get_text(" ", strip=True)
        else:
            content = ""
        if content:
            texts.append({"source": url, "text": content})
    return texts


def counterpart_bands(counterpart, dr):
    result = []
    for row in counterpart["matched"] or []:
        ids = (row.get("ids") or "").split("|")
        radio = [s for s in ids if re.match(
            r"NVSS|VLASS|FIRST|PKS|CRATES|GB6|GB |PMN|AT20G|SUMSS|TXS|RFC|CLASS |S[45] ", s)]
        optical = [s for s in ids if re.match(r"Gaia|SDSS|QSO B|2QZ|6dF|LAMOST|USNO|APM |ZTF", s)]
        phot = [r for r in counterpart["photometry"] or []
                if str(r.get("oidref")) == str(row["oid"]) and str(r.get("filter")) in
                {"u", "g", "r", "i", "z", "G", "U", "B", "V", "R", "I", "BP", "RP"}]
        r = {"main_id": row["main_id"], "type": row["otype"], "radio_ids": radio,
             "optical_ids": optical, "optical_photometry": phot}
        if dr and row.get("ra") is not None and row.get("dec") is not None:
            r["dr4_separation_deg"] = angle(float(dr["RAJ2000"]), float(dr["DEJ2000"]),
                                             float(row["ra"]), float(row["dec"]))
        result.append(r)
    return result


def primary_aliases(lead, dr):
    aliases = lead.get("aliases") or []
    if isinstance(aliases, str):
        aliases = [aliases]
    seq = [lead["name"], lead.get("assoc"), *aliases]
    if dr:
        seq += [dr.get("ASSOC1"), dr.get("ASSOC2"), dr.get("ASSOC_FGL"), dr.get("Source_Name")]
    return [s for s in dedup(seq) if len(s) > 3 and s.lower() not in ("none", "unknown", "unassociated")]


def review_one(item, catalog, telegrams, papers, notes):
    audit = []
    dr = identify(item, catalog)
    aliases = primary_aliases(item, dr)
    counterpart = simbad(aliases, audit)
    for source in counterpart["matched"] or []:
        id_list = str(source.get("ids") or "").split("|")
        aliases += [source.get("main_id")]
        for prefix in ("PKS", "B2 ", "S5 ", "S4 ", "BL Lac", "1ES", "3C ", "PG ", "Mrk", "OJ ", "TXS"):
            aliases += [s for s in id_list if s.startswith(prefix)][:2]
    aliases = dedup(aliases)
    # Variants such as '4FGL J1555.7+1111' and 'J1555.7+1111' both occur.
    matchers = [(s, norm(s)) for s in aliases]
    matchers += [(s, norm(s.split(" ", 1)[1])) for s in aliases
                 if re.match(r"^[1-4]FGL J", s) and " " in s]
    def hits(txt):
        x = norm(txt)
        return dedup([s for s, n in matchers if len(n) >= 4 and n in x])
    q = "(" + " OR ".join('"' + s.replace('"', '').replace('\\', '') + '"' for s in aliases[:12]) + ")"
    tquery = q + ' AND (periodic* OR oscillat* OR binar* OR flare* OR outburst* OR "state change" OR "changing look")'
    token = os.environ.get("ADS_API_TOKEN") or os.environ.get("ADS_TOKEN")
    ads_result = ads(tquery, audit, token)
    ads_atel = ads(q + " AND bibstem:ATel", audit, token)
    allrefs = counterpart.get("references")
    cds_refs = [{"bibcode": z.get("bibcode"), "title": z.get("title"), "year": z.get("year"),
                 "ads_url": "https://ui.adsabs.harvard.edu/abs/" + str(z.get("bibcode")) + "/abstract"}
                for z in allrefs or [] if re.search(
                    r"period|oscillat|flar|outburst|transition|state|chang|binary|binaries",
                    z.get("title") or "", re.I) or "ATel" in str(z.get("bibcode"))]
    # Public publication lookup remains available if ADS credentials are absent.
    crossrefs = []
    for alias in aliases[:4]:
        obj = get_json(CROSSREF, audit, params={"query.bibliographic":
                       alias + " gamma ray periodicity flaring high low state transition", "rows": 12})
        for x in (obj or {}).get("message", {}).get("items", []):
            title = BeautifulSoup(" ".join(x.get("title") or []), "html.parser").get_text(" ", strip=True)
            if hits(title):
                crossrefs.append({"title": title, "doi": x.get("DOI"), "matched": hits(title)})
    paper_hits = []
    for p, txt in papers:
        found = hits(txt)
        if found:
            lines = [line.strip()[:220] for line in txt.splitlines() if hits(line)][:4]
            paper_hits.append({"doi": p["doi"], "title": p["title"], "aliases": found,
                               "table_or_text_lines": lines})
    atels = [{"number": k, "title": title, "url": ATEL + f"?read={k}"}
             for k, title in telegrams.items() if hits(title)]
    note_hits = [{"url": n["source"], "aliases": hits(n["text"])} for n in notes if hits(n["text"])]
    bands = counterpart_bands(counterpart, dr)
    gaps = []
    if not any(b["radio_ids"] or b["optical_ids"] or b["optical_photometry"] for b in bands):
        gaps.append("no identified optical or radio catalog cross-match recovered")
    if not any(n["source"] != DR4 for n in notes):
        gaps.append("4FGL-DR4 release manuscript/notes unavailable")
    if any("unavailable" in z["status"] for z in (ads_result, ads_atel)):
        gaps.append("ADS search unavailable; SIMBAD-linked ADS bibcodes and Crossref are partial")
    if any(z.get("coverage") for z in (ads_result, ads_atel)):
        gaps.append("public ADS query limited to first rendered search page; inspect pagination")
    if not dr:
        gaps.append("4FGL-DR4 identification unavailable or ambiguous")
    if not counterpart["matched"]:
        gaps.append("no exact SIMBAD association/identification recovered")
    if allrefs is None or len(allrefs) >= 600:
        gaps.append("CDS object bibliography unavailable or truncated")
    if not any(p["penil"] and txt for p, txt in papers):
        gaps.append("no searchable Peñil et al. manuscript retrieved")
    if any(not txt for p, txt in papers):
        gaps.append("some comparison manuscripts lacked searchable full text")
    if len(telegrams) < 10000:
        gaps.append("ATel title archive coverage limited; consult object bibliography and archive")
    return {"name": item["name"], "lead": item, "status": "incomplete" if gaps else "requires_review",
            "aliases": aliases, "dr4": dr, "counterpart": counterpart,
            "ads": ads_result, "ads_atel": ads_atel, "simbad_ads_evidence": cds_refs,
            "dr4_note_hits": note_hits, "optical_radio": bands,
            "crossref": crossrefs, "atel_title_hits": atels, "periodicity_table_hits": paper_hits,
            "coverage_gaps": gaps, "receipts": audit}


def collect(leads, output_dir):
    """Write exact per-lead bibliographic evidence; no call if there are no leads."""
    dest = Path(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    if not leads:
        (dest / "index.json").write_text("[]\n")
        return []
    audit = []
    catalog = dr4_catalog(audit)
    notes = dr4_notes(audit)
    telegrams = atel_archive(audit)
    papers = topic_papers(audit)
    summary = {"dr4_rows": len(catalog), "atel_title_count": len(telegrams),
               "atel_number_range": [min(telegrams), max(telegrams)] if telegrams else None,
               "papers": [p for p, _ in papers], "dr4_notes": [n["source"] for n in notes],
               "receipts": audit}
    (dest / "shared-evidence.json").write_text(json.dumps(summary, indent=2) + "\n")
    def one(item):
        try:
            return review_one(item, catalog, telegrams, papers, notes)
        except (ValueError, KeyError, TypeError, OSError, ArithmeticError) as err:
            return {"name": item["name"], "status": "incomplete", "lead": item,
                    "coverage_gaps": [type(err).__name__ + ": " + str(err)[:300]]}
    with ThreadPoolExecutor(max_workers=3) as pool:
        data = list(pool.map(one, leads))
    if _UI_DRIVER is not None:
        _UI_DRIVER.quit()
    index = []
    for i, row in enumerate(data):
        file = f"lead-{i + 1:04d}.json"
        (dest / file).write_text(json.dumps(row, indent=2, allow_nan=False) + "\n")
        index.append({"name": row["name"], "status": row["status"], "file": file,
                      "coverage_gaps": row["coverage_gaps"]})
    (dest / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    return index

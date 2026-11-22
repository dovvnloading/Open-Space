"""Targeted searches of object bibliographies, notices and historical spectra."""
import csv
import io
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from tempfile import TemporaryDirectory
from pathlib import Path
import sys
from urllib.parse import quote

from astropy.io.votable import parse
from astropy.table import Table
from bs4 import BeautifulSoup
from common import OUT, PREVIOUS, fetch, status

ALIASES = ["KUG 1256+241", "RX J1258.8+2355", "2MASX J12585149+2355268",
           "SDSS J125851.45+235526.6", "RBS 1195", "2CXO J125851.4+235526",
           "1RXS J125851.4+235532", "XMMSL2 J125851.2+235522"]
TERM = re.compile(r"changing.look|state (?:change|transition)|broad.{0,25}(?:vanish|disappear|emerg)|turn(?:ed)?[ -](?:on|off)", re.I)
CODE = re.compile(r"(?<!\w)(?:19|20)\d\d[A-Za-z&.][A-Za-z&.0-9]{13}[A-Za-z.](?!\w)")
records = []


def clean(s):
    return " ".join(str(s).split())


def ned():
    from astroquery.ipac.ned import Ned
    # Use the installed client to construct the same legacy NED references request.
    args = Ned.get_table(ALIASES[2], table="references", get_query_payload=True)
    endpoints = ["https://ned.ipac.caltech.edu/cgi-bin/nph-refsearch",
                 "https://ned.ipac.caltech.edu/cgi-bin/objsearch"]
    failures = []
    for url in endpoints:
        try:
            raw = fetch("NED.references", url, args, timeout=22)
            vot = parse(io.BytesIO(raw))
            tables = [t.to_table() for t in vot.iter_tables()]
            codes = sorted(set(CODE.findall(raw.decode(errors="replace"))))
            if not tables or not codes:
                raise ValueError("No identifiable bibliography table with bibcodes")
            records.append(dict(search="NED bibliography", references=codes))
            status("NED bibliography", "retrieved", f"{len(codes)} bibcodes in the object bibliography")
            return
        except Exception as exc:
            failures.append(clean(exc)[:200])
    raise ValueError(" / ".join(failures))


def ads():
    old = json.loads((PREVIOUS / "outputs/seyfert_literature.json").read_text())
    papers = old["results"]["simbad_seyfert"]["all_bibliography"]
    refs = {row["bibcode"]: row["title"] for row in papers}
    for r in records:
        if r["search"] == "NED bibliography":
            refs.update({id_:refs.get(id_, "") for id_ in r["references"]})
    # ADS public abstracts can be inspected without credentials even when API search needs a token.
    q = " OR ".join('full:"' + alias + '"' for alias in ALIASES)
    try:
        result = json.loads(fetch("ADS.aliases", "https://api.adsabs.harvard.edu/v1/search/query",
                                  {"q":"("+q+")", "fl":"bibcode,title,abstract", "rows":200},timeout=20))
        for row in result["response"]["docs"]:
            refs[row["bibcode"]] = " ".join(row["title"])
        status("ADS full-text name search", "retrieved", f'{result["response"]["numFound"]} hits for the exact names')
    except Exception as exc:
        status("ADS full-text name search", "unavailable", clean(exc)[:260])
    def abstract(item):
        code, title = item
        url = "https://ui.adsabs.harvard.edu/abs/" + quote(code,safe="") + "/abstract"
        record = dict(search="ADS object bibliography", bibcode=code, title=title, url=url)
        try:
            doc = BeautifulSoup(fetch("ADS."+code, url, timeout=15), "html.parser")
            tag = doc.select_one('meta[name="description"], meta[property="og:description"]')
            abstract_text = clean(tag.get("content","")) if tag else ""
            if len(abstract_text) < 60 or "digital library portal" in abstract_text.lower():
                raise ValueError("No article-specific abstract metadata")
            record["abstract_read"] = True
            record["transition_terms"] = TERM.findall(abstract_text)
            record["names"] = [a for a in ALIASES if a.lower() in abstract_text.lower()]
            if record["transition_terms"]:
                match = TERM.search(abstract_text)
                record["context"] = abstract_text[max(0,match.start()-80):match.end()+80]
        except Exception as exc:
            record["error"] = clean(exc)[:200]
        return record
    with ThreadPoolExecutor(max_workers=4) as pool:
        inspected = list(pool.map(abstract, refs.items()))
    records.extend(inspected)
    n = sum(bool(r.get("abstract_read")) for r in inspected)
    status("ADS object bibliography", "retrieved" if n else "unavailable",
           f"{n}/{len(refs)} linked ADS abstracts readable; title list also deposited")


def tns():
    sys.path.insert(0, str(PREVIOUS.resolve()))
    from tns_html import parse_page
    endpoint = "https://www.wis-tns.org/search"
    def search(key, ra, dec):
        params = dict(ra=ra, decl=dec, radius=35, coords_unit="arcsec", num_page=50)
        rows, next_page, empty, pages = parse_page(fetch("TNS."+key, endpoint, params, timeout=20))
        if next_page or pages - {0} or not (rows or empty):
            raise ValueError("Incomplete cone return")
        return rows
    rows = search("nucleus","12:58:51.49","+23:55:26.8")
    records.append(dict(search="TNS official cone", objects=rows))
    if not rows and not any(r["tns_slug"] == "2018cow" for r in search("control","16:16:00.22","+22:16:04.8")):
        raise ValueError("Empty cone with unsuccessful recovery of AT 2018cow")
    status("TNS official cone", "retrieved" if rows else "no_match", f"{len(rows)} objects in 35 arcsec")


def telegrams():
    def search(name):
        soup = BeautifulSoup(fetch("ATel."+name, "https://www.astronomerstelegram.org/",
                                   {"search":name}, timeout=18), "html.parser")
        notices = []
        for a in soup.find_all("a", href=True):
            m = re.search(r"[?&]read=(\d+)", a["href"])
            if m:
                notices.append(dict(id=m[1], title=clean(a.get_text(" ",strip=True))))
        page = clean(soup.get_text(" ",strip=True))
        return notices, page
    # A known, distinctive telegram guards against a home page silently ignoring the query.
    ctrl, _ = search("XMMSL1 J063045.9-603110")
    if not any(r["id"] == "3811" for r in ctrl):
        raise ValueError("Response did not recover ATel 3811 for its named comparison source")
    for alias in ALIASES[:5]:
        hits, page = search(alias)
        records.append(dict(search="ATel", alias=alias, links=hits[:100]))
        if hits == ctrl or (not hits and not re.search(r"no (?:results|matches|telegrams)|0 results",page,re.I)):
            raise ValueError("Unvalidated search response; see returned links")
    status("ATel search", "searched", "Five aliases; returned links require object-level review")


def papers():
    # Read table contexts from the discovery paper. Keep short excerpts, including headers if near the object.
    urls = ["https://arxiv.org/pdf/astro-ph/0106311",
            "https://www.aanda.org/articles/aa/pdf/2008/20/aa8552-07.pdf"]
    for url in urls:
        try:
            raw = fetch("paper."+url.split("/")[-1], url,timeout=35)
            if not raw.startswith(b"%PDF"):
                raise ValueError("Expected a PDF")
            with TemporaryDirectory() as d:
                pdf = Path(d) / "paper.pdf"
                pdf.write_bytes(raw)
                txt = subprocess.check_output(["pdftotext","-layout",str(pdf),"-"],timeout=40).decode(errors="replace")
            lines = txt.splitlines()
            matches = [dict(line=i+1, text=l.strip(), prior=lines[max(0,i-1)].strip())
                       for i,l in enumerate(lines) if re.search(r"1258\s*\.?\s*8\s*\+\s*2355|1256\s*\+\s*241|125851[.\d]+\+2355|RBS\s+1195",l,re.I)]
            records.append(dict(search="historical paper", url=url, rows=matches,
                                table_headers=[clean(l) for l in lines if re.search(r"FWHM|H.?β.*H.?α",l,re.I) and len(l)<210][:18]))
            status("Historical paper "+url.rsplit("/",1)[1], "retrieved", f"{len(matches)} target-name lines")
        except Exception as exc:
            status("Historical paper "+url.rsplit("/",1)[1], "unavailable", clean(exc)[:260])


def published_tables():
    surveys = ["J/MNRAS/327/673", "J/A+A/484/897", "J/ApJS/267/19", "J/AJ/123/578",
               "J/AJ/119/2540", "J/A+A/416/35"]
    def one(cat):
        try:
            args = {"-source":cat, "-c":"194.714542 +23.924111", "-c.rs":45,
                    "-out.all":1, "-out.max":100, "-out.add":"_r"}
            raw = fetch("VizieR.paper."+cat, "https://vizier.cds.unistra.fr/viz-bin/votable",
                        args, timeout=25)
            vot = parse(io.BytesIO(raw))
            rows = []
            tables = list(vot.iter_tables())
            for table in tables:
                data = table.to_table(use_names_over_ids=True)
                for row in data:
                    rows.append({k:str(v) for k,v in zip(data.colnames,row)})
            records.append(dict(search="published table",catalog=cat,rows=rows,
                                columns={f.name:{"unit":str(f.unit or ""),"description":f.description}
                                         for t in tables for f in t.fields}))
            state = "retrieved" if rows else "no_match" if tables else "unavailable"
            status("VizieR "+cat,state,f"{len(rows)} rows in 45 arcsec; individual IDs and offsets in literature.json")
        except Exception as exc:
            status("VizieR "+cat, "unavailable", str(exc)[:240])
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(one,surveys))


def run():
    def safe(name, fn):
        try:
            fn()
        except Exception as exc:
            status(name, "unavailable", clean(exc)[:400])
    jobs = [("NED bibliography",ned),("TNS official cone",tns),("ATel search",telegrams),("Historical papers",papers), ("Published tables",published_tables)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda pair:safe(*pair),jobs))
    safe("ADS object bibliography",ads)
    (OUT / "literature.json").write_text(json.dumps(records, indent=2, ensure_ascii=False)+"\n")

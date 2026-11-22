"""Audit the prior literature on each independent BAT/MAXI survivor.

Run from the survey's Actions workflow. Inputs are JSONL dictionaries with source,
ra, dec, period_days, uncertainty_days and an optional aliases array. Public query
results are archived losslessly (gzip); the digest is for the HTTP response bytes.
An empty search or a network error can never establish that a period is new.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode

import requests


VIZIERS = ("https://vizier.cds.unistra.fr/viz-bin/asu-xml",
           "https://vizier.cfa.harvard.edu/viz-bin/asu-xml")
SIMBADS = ("https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
           "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync")
CROSSREF = "https://api.crossref.org/works"
ATOMS = {"a": "http://www.w3.org/2005/Atom",
         "o": "http://a9.com/-/spec/opensearch/1.1/"}
# Indexed series names, author checks and static VizieR IDs. Additional editions
# are resolved from the primary papers' A&A volume and article numbers via Crossref.
CATALOGUES = {
    "Ritter-Kolb": ([], [], ["B/cb"]),
    "Fortin HMXB": (["Fortin"], ["high", "mass"], ["J/A+A/671/A149"]),
    "Fortin LMXB": (["Fortin"], ["low", "mass"], []),
    "Neumann HMXB": (["Neumann"], ["high", "mass"], ["J/A+A/677/A134"]),
    "Avakyan-Neumann LMXB": (["Avakyan", "Neumann"], ["low", "mass"], ["J/A+A/675/A199"]),
}


def dump(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


class Archive:
    def __init__(self, output):
        self.output = output
        self.blobs = output / "responses"
        self.blobs.mkdir(parents=True, exist_ok=True)
        self.journal = output / "queries.jsonl"
        self.cache = {}
        self.pdf_text = {}
        self.pdf_program = shutil.which("pdftotext")
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "Open-Space/BAT-MAXI-period-survey (astronomy research)"
        if self.journal.exists():
            for line in self.journal.read_text().splitlines():
                row = json.loads(line)
                if "sha256" in row:
                    self.cache[row["url"]] = row

    def fetch(self, url, params=None, headers=None):
        full = url + ("?" + urlencode(params) if params else "")
        previous = self.cache.get(full)
        if previous:
            raw = gzip.decompress((self.blobs / (previous["sha256"] + ".gz")).read_bytes())
            if hashlib.sha256(raw).hexdigest() != previous["sha256"]:
                raise ValueError("archived query checksum mismatch")
            return raw, previous
        row = dict(url=full, utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        for attempt in range(2):
            try:
                response = self.session.get(full, headers=headers, timeout=(7, 30))
                row["status"] = response.status_code
                response.raise_for_status()
                raw = response.content
                digest = hashlib.sha256(raw).hexdigest()
                row.update(sha256=digest, nbytes=len(raw))
                (self.blobs / (digest + ".gz")).write_bytes(gzip.compress(raw, mtime=0))
                self.cache[full] = row
                with self.journal.open("a") as stream:
                    stream.write(json.dumps(row) + "\n")
                return raw, row
            except requests.RequestException as exc:
                row["error"] = str(exc)[:500]
                if row.get("status") in (400, 401, 403, 404):
                    break
                if attempt == 0:
                    time.sleep(2)
        with self.journal.open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        raise ValueError(row.get("error", "request failed"))


def norm(value):
    # Preserve a coordinate's sign. Ignore separators and the common SIMBAD NAME prefix.
    return re.sub(r"[^a-z0-9+-]", "", re.sub(r"^NAME\s+", "", str(value), flags=re.I)
                  .lower().replace("−", "-").replace("–", "-"))


def split_aliases(values):
    return [w.strip() for v in values for w in re.split(r"[;,|]", str(v)) if w.strip()]


def uniq(values):
    return list(dict.fromkeys(v for v in values if v))


def lit_aliases(values):
    # Survey aliases plus usual XRB designations get priority in full-text queries.
    priorities = ("4u", "3u", "2s", "sax", "igr", "swift", "xte", "rx", "gro", "maxi")
    vs = uniq(values)
    return uniq(vs[:2] + [v for v in vs[2:] if norm(v).startswith(priorities)] + vs[2:])[:9]


def sky_distance(ra, dec, other_ra, other_dec):
    ra, dec, other_ra, other_dec = map(math.radians, [ra, dec, other_ra, other_dec])
    a = math.sin((dec-other_dec)/2)**2 + math.cos(dec)*math.cos(other_dec)*math.sin((ra-other_ra)/2)**2
    return math.degrees(2*math.asin(min(1., math.sqrt(a)))) * 60


def literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def tap(archive, adql, maxrows):
    params = dict(request="doQuery", lang="adql", format="csv", maxrec=str(maxrows), query=adql)
    errors = []
    for url in SIMBADS:
        try:
            raw, meta = archive.fetch(url, params)
            reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
            if not reader.fieldnames or reader.fieldnames[0].startswith("<") or "ERROR" in raw[:250].decode(errors="ignore"):
                raise ValueError("invalid TAP response: " + raw[:300].decode(errors="replace"))
            return list(reader), meta
        except (ValueError, csv.Error) as exc:
            errors.append(str(exc))
    raise ValueError(" | ".join(errors))


def simbad(archive, c):
    aliases = uniq([c["source"]] + list(c.get("aliases") or []))
    coords = c.get("ra") is not None and c.get("dec") is not None
    select = "SELECT DISTINCT b.oid,b.main_id,b.ra,b.dec FROM basic AS b "
    rows, meta = tap(archive, select + " JOIN ident AS i ON i.oidref=b.oid WHERE i.id IN (" +
                     ",".join(map(literal, aliases[:64])) + ")", 128)
    how = "exact identifier"
    if not rows and coords:
        ra, dec = float(c["ra"]), float(c["dec"])
        rows, meta = tap(archive, select + " WHERE 1=CONTAINS(POINT('ICRS', b.ra,b.dec),"
                         + f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.016666667))", 128)
        how = "one-arcminute cone"
    result = dict(lookup=how, matches=rows, query=meta)
    if len(rows) != 1:
        result["status"] = "ambiguous or missing SIMBAD identification"
        return result, aliases
    obj = rows[0]
    if coords and obj.get("ra") and obj.get("dec"):
        result["offset_arcmin"] = sky_distance(float(c["ra"]), float(c["dec"]), float(obj["ra"]), float(obj["dec"]))
        if result["offset_arcmin"] > 3:
            result["status"] = "SIMBAD identifier-coordinate disagreement"
            return result, aliases
    oid = int(obj["oid"])
    ids, m_ids = tap(archive, f"SELECT id FROM ident WHERE oidref={oid}", 2000)
    aliases = uniq([aliases[0], obj["main_id"]] + aliases[1:] + [r["id"] for r in ids])
    result.update(status="identified", aliases_query=m_ids, aliases=aliases)
    try:
        refs, m_refs = tap(archive, "SELECT DISTINCT r.bibcode,r.title FROM ref AS r "
                           f"JOIN has_ref AS h ON h.oidbibref=r.oidbib WHERE h.oidref={oid}", 30000)
        timing = re.compile(r"orbit|superorb|period|ephemer|recurr|precess|outburst|modulat|timing|cycle", re.I)
        result.update(bibliography_query=m_refs, references_total=len(refs), bibliography_truncated=len(refs)>=30000,
                      timing_references=[r for r in refs if timing.search(r.get("title", ""))])
    except ValueError as exc:
        result["bibliography_error"] = str(exc)
    return result, aliases


def asu_tables(raw):
    root = ET.fromstring(raw)
    if root.tag.split("}")[-1] != "VOTABLE":
        raise ValueError("expected a VOTable")
    for node in root.iter():
        node.tag = node.tag.split("}")[-1]
    errors = [x.text for x in root.iter("INFO") if x.get("value", "").upper() == "ERROR"]
    if errors:
        raise ValueError(str(errors))
    tables = []
    for table in root.iter("TABLE"):
        fields = [dict(name=f.get("name", ""), unit=f.get("unit", ""), ucd=f.get("ucd", ""),
                       desc=" ".join((f.findtext("DESCRIPTION") or "").split())) for f in table.findall("FIELD")]
        body = table.find("DATA/TABLEDATA")
        if not fields or body is None:
            continue
        rows = []
        for tr in body.findall("TR"):
            row = {field["name"]: (v.text or "").strip() for field, v in zip(fields, tr.findall("TD"))}
            if any(row.values()):
                rows.append(row)
        tables.append(dict(name=table.get("name", ""), fields=fields, rows=rows))
    if not any(t["rows"] for t in tables):
        raise ValueError("no populated catalogue tables")
    return tables


def xref_rows(archive, params):
    raw, meta = archive.fetch(CROSSREF, params)
    payload = json.loads(raw)["message"]
    return payload["items"], meta


def catalogues(archive):
    inventory, tables = {}, []
    for family, (authors, terms, fallback) in CATALOGUES.items():
        audit = dict(papers=[], downloads=[])
        ids = []
        if authors:
            params = {"query.author": authors[0], "query.bibliographic": "Galactic " + "-".join(terms)
                      + " X-ray binary catalogue", "rows": "12"}
            try:
                works, m = xref_rows(archive, params)
                audit["discovery_query"] = m
                for w in works:
                    title = (w.get("title") or [""])[0]
                    names = [p.get("family", "").lower() for p in w.get("author", [])]
                    if not any(a.lower() in names for a in authors):
                        continue
                    if not all(v in title.lower() for v in terms + ["catalog"]):
                        continue
                    volume, page = w.get("volume", ""), w.get("article-number") or w.get("page", "")
                    aid = re.match(r"A[0-9]+", str(page), re.I)
                    journal = " ".join(w.get("container-title", []))
                    if aid and str(volume).isdigit() and ("astrophys" in journal.lower() or "A&A" in journal):
                        vid = f"J/A+A/{volume}/{aid[0].upper()}"
                        ids.append(vid)
                        audit["papers"].append(dict(title=title, doi=w.get("DOI"), vizier=vid))
            except (ValueError, KeyError) as exc:
                audit["discovery_error"] = str(exc)
        for vid in uniq(ids + fallback):
            errors = []
            for service in VIZIERS:
                try:
                    raw, meta = archive.fetch(service, {"-source": vid, "-out.all": "1", "-out.max": "50000",
                                                        "-out.add": "_RAJ,_DEJ", "-out.form": "VOTable"})
                    parsed = asu_tables(raw)
                    # Require the actual catalogue/resource title to identify the intended
                    # authors and source population, including when a fallback was supplied.
                    root = ET.fromstring(raw)
                    title_lines = []
                    for node in root.iter():
                        if node.tag.split("}")[-1] == "RESOURCE":
                            for child in node:
                                if child.tag.split("}")[-1] == "DESCRIPTION":
                                    title_lines.append(" ".join(" ".join(child.itertext()).split()))
                    identity = " ".join(title_lines)
                    def identifies(title):
                        title = re.sub(r"[^a-z0-9 ]", " ", title.lower())
                        if family == "Ritter-Kolb":
                            return bool(re.search(r"\b(ritter|kolb)\b", title) and
                                        re.search(r"cataclysmic|lmxb|low.mass", title))
                        population = "high" if "HMXB" in family else "low"
                        return (any(re.search(r"\b"+a.lower()+r"\b", title) for a in authors)
                                and bool(re.search(population+r"\s*mass|\b"+population[0]+r"mxbs?\b", title)))
                    readme = None
                    if not identifies(identity):
                        raw_title, readme = archive.fetch("https://cdsarc.cds.unistra.fr/ftp/"+vid+"/ReadMe")
                        header = raw_title.decode(errors="replace").split("Description:")[0].split("File Summary:")[0]
                        identity = " ".join(header.split()[:250])
                    if not identifies(identity):
                        raise ValueError("catalogue author or population did not match its official title")
                    tables.extend((family, vid, t) for t in parsed)
                    audit["downloads"].append(dict(vizier=vid, query=meta, official_title=identity, readme=readme,
                        tables=[dict(name=t["name"], rows=len(t["rows"]), fields=t["fields"]) for t in parsed]))
                    break
                except (ValueError, ET.ParseError) as exc:
                    errors.append(str(exc))
            else:
                audit.setdefault("download_errors", {})[vid] = errors
        audit["status"] = "downloaded" if audit["downloads"] else "unresolved"
        inventory[family] = audit
    return inventory, tables


def angle(value, unit, ra=False):
    try:
        if ":" not in value and " " not in value.strip():
            v = float(value)
            return v * 15 if ra and unit.startswith("h") else v
        vals = list(map(float, re.split(r"[ :]+", value.strip())))
        v = abs(vals[0])+vals[1]/60+(vals[2] if len(vals)>2 else 0)/3600
        return (-v if value.strip().startswith("-") else v) * (15 if ra and not unit.startswith("d") else 1)
    except (ValueError, IndexError):
        return None


def close_catalogues(tables, aliases, candidate):
    alias_keys = {norm(v) for v in aliases if len(norm(v)) >= 4}
    matched = []
    for family, vid, t in tables:
        fields = t["fields"]
        names = [f["name"] for f in fields if re.search("name|alias|designation|^object$|^source$|^id$",
                 f["name"], re.I) or (f["ucd"].startswith("meta.id") and "ref" not in f["name"].lower())]
        ra = next((f for f in fields if f["name"] in ("_RAJ2000", "_RAJ", "RAJ2000", "RA_ICRS")), {})
        dec = next((f for f in fields if f["name"] in ("_DEJ2000", "_DEJ", "DEJ2000", "DE_ICRS")), {})
        pfields = []
        for f in fields:
            desc = f["desc"].lower()
            if not ("period" in desc or "time.period" in f["ucd"] or re.match("p(orb|sup)|orb.per", f["name"], re.I)):
                continue
            if re.search(r"^([erufl]_|[+-])", f["name"]) or re.search("reference|error|uncertainty|flag|limit flag", desc):
                continue
            pfields.append(f)
        for row in t["rows"]:
            ids = split_aliases([row.get(k, "") for k in names])
            common = [v for v in ids if norm(v) in alias_keys]
            d = None
            if ra and dec and candidate.get("ra") is not None and candidate.get("dec") is not None:
                r, b = angle(row.get(ra["name"], ""), ra["unit"], True), angle(row.get(dec["name"], ""), dec["unit"])
                if r is not None and b is not None:
                    d = sky_distance(float(candidate["ra"]), float(candidate["dec"]), r, b)
            if not common and (d is None or d > 1):
                continue
            periods = []
            for f in pfields:
                val = row.get(f["name"], "")
                unit = f["unit"].strip().lower()
                scale = {"d":1., "day":1., "h":1/24., "hr":1/24., "s":1/86400., "min":1/1440.,
                         "yr":365.25, "a":365.25}.get(unit)
                p = dict(field=f["name"], description=f["desc"], value=val, unit=f["unit"])
                if not val:
                    continue
                mt = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?\d+)?)\s*", val)
                if mt and scale and float(mt[1]) > 0:
                    days = float(mt[1])*scale
                    p["period_days"] = days
                    ratio = float(candidate["period_days"])/days
                    p["ratio_measured_to_catalogue"] = ratio
                    tolerance = max(.02, 3*float(candidate.get("uncertainty_days") or 0)/float(candidate["period_days"]))
                    harmonics = (.25, 1/3, .5, 2/3, 1., 1.5, 2., 3., 4.)
                    p["near_ratios"] = [k for k in harmonics if abs(ratio/k-1) <= tolerance]
                periods.append(p)
            matched.append(dict(family=family, vizier=vid, table=t["name"], identifiers=ids, matched_identifiers=common,
                                arcmin=d, periods=periods, row={k:v for k,v in row.items() if v}))
    return matched


PERIOD_TEXT = re.compile(r"(?<![\w.])(\d{1,6}(?:\.\d*)?)(?:\(\d+\)|\s*(?:±|\+/-)\s*\d+(?:\.\d*)?)?\s*"
                         r"(?:-\s*)?(days?\b|d\b|hours?\b|hr?s?\b|years?\b|yrs?\b)", re.I)


def passages(text, target):
    text = " ".join(re.sub(r"<[^>]+>", " ", text or "").split())
    hits = []
    for m in PERIOD_TEXT.finditer(text):
        start, end = m.span()
        unit = m[2].lower()
        days = float(m[1]) * (1/24 if unit.startswith("h") else 365.25 if unit.startswith("y") else 1)
        if days <= 0:
            continue
        ratio = target/days
        if any(abs(ratio/k-1)<.04 for k in (.25,1/3,.5,2/3,1.,1.5,2.,3.,4.)):
            hits.append(dict(text=" ".join(text[:start].split()[-8:]) + " " + text[start:end]
                             + " " + " ".join(text[end:].split()[:8]), period_days=days, ratio=ratio))
    return hits[:8]


def paper_text(archive, identifier):
    """Read a primary arXiv manuscript on the runner; retain digest and snippets."""
    if identifier in archive.pdf_text:
        return archive.pdf_text[identifier]
    url = "https://arxiv.org/pdf/" + identifier
    meta = dict(url=url)
    pages = []
    try:
        if not archive.pdf_program:
            raise ValueError("pdftotext unavailable")
        response = archive.session.get(url, timeout=(7, 35))
        response.raise_for_status()
        raw = response.content
        if not raw.startswith(b"%PDF-"):
            raise ValueError("response is not a PDF")
        meta.update(sha256=hashlib.sha256(raw).hexdigest(), nbytes=len(raw))
        with tempfile.TemporaryDirectory() as temp:
            filename = Path(temp) / "paper.pdf"
            filename.write_bytes(raw)
            proc = subprocess.run([archive.pdf_program, "-layout", str(filename), "-"],
                                  capture_output=True, timeout=25, check=True)
            pages = proc.stdout.decode("utf-8", errors="replace").split("\f")
            if not any(v.strip() for v in pages):
                raise ValueError("no machine-readable paper text")
        meta["pages_with_text"] = sum(bool(v.strip()) for v in pages)
    except (requests.RequestException, ValueError, OSError, subprocess.SubprocessError) as exc:
        meta["error"] = str(exc)[:350]
    archive.pdf_text[identifier] = pages, meta
    return pages, meta


def manuscript_hits(archive, identifiers, aliases, period):
    # Match named sources and scan the same physical page, allowing a period in
    # a table to use units in its column header. Retain page numbers and audit
    # all PDFs, including those which could not be fetched or converted.
    patterns = []
    for a in aliases:
        a = re.sub(r"^NAME\s+", "", a, flags=re.I).strip()
        chunks = re.findall(r"[a-zA-Z]+|[0-9]+|[+-]", a)
        if len(norm(a)) >= 5 and chunks:
            patterns.append(re.compile(r"(?<![a-zA-Z0-9])" + r"[\s_.]*".join(re.escape(v) for v in chunks)
                                       + r"(?![a-zA-Z0-9])", re.I))
    results = []
    for identifier in uniq(identifiers):
        pages, meta = paper_text(archive, identifier)
        page_hits = []
        for n, page in enumerate(pages, 1):
            page = page.replace("−", "-").replace("–", "-")
            matches = [m for pat in patterns for m in pat.finditer(page)]
            if not matches:
                continue
            # A number in another row or on the same page is a lead for human
            # verification, not a published measurement of the requested source.
            snippets = []
            for m in sorted(matches, key=lambda z: z.start()):
                line = page.count("\n", 0, m.start())
                lines = page.splitlines()
                passage = " ".join(lines[line].split()[:22])
                if passage and passage not in snippets:
                    snippets.append(passage)
                if len(snippets) >= 3:
                    break
            page_hits.append(dict(page=n, source_context=snippets,
                                  unit_bearing_periods_on_same_page=passages(page, period)))
        results.append(dict(manuscript=identifier, retrieval=meta, source_pages=page_hits))
    return results


def publications(archive, aliases, target):
    used = lit_aliases(aliases)
    names = {norm(a) for a in used if len(norm(a)) >= 4}
    papers, queries = {}, {}
    identifiers = []
    def add(key, title, abstract, author, url, branch, by_source=False):
        combined = norm(title + " " + re.sub(r"<[^>]+>", " ", abstract or ""))
        matched = sorted(v for v in names if v in combined)
        hits = passages(title + " " + (abstract or ""), target)
        d = papers.setdefault(key, dict(title=title, author=author, url=url, served_by=[], alias_hits=[], passages=[]))
        d["served_by"] = uniq(d["served_by"]+[branch])
        d["alias_hits"] = uniq(d["alias_hits"] + matched)
        if hits:
            d["passages"] = hits
        if by_source:
            d["object_search_hit"] = True
    def query(name, cb):
        try:
            queries[name] = cb()
        except (ValueError, KeyError, ET.ParseError, TypeError) as exc:
            queries[name] = dict(error=str(exc))

    def arxiv(corbet, survey=False):
        aq = " OR ".join('all:"'+a.replace('"', "")+'"' for a in used)
        term = f"({aq})" + (" AND au:Corbet" if corbet else "")
        if survey:
            term = ("au:Corbet AND (ti:orbital OR ti:superorbital OR ti:period "
                    "OR ti:periodic OR ti:recurrence OR ti:modulation) AND "
                    "(all:pulsar OR all:binary OR all:binaries OR all:Swift OR all:RXTE OR all:MAXI)")
        page_size = 200 if survey else 100
        raw, meta = archive.fetch("https://export.arxiv.org/api/query",
                                  dict(search_query=term, start=0, max_results=page_size, sortBy="relevance"))
        root = ET.fromstring(raw)
        if root.tag != "{"+ATOMS["a"]+"}feed":
            raise ValueError("unexpected arXiv response")
        entries = root.findall("a:entry", ATOMS)
        for ent in entries:
            title = " ".join(ent.findtext("a:title", "", ATOMS).split())
            url = ent.findtext("a:id", "", ATOMS)
            if title.lower() == "error" or "api/errors" in url:
                raise ValueError(ent.findtext("a:summary", "arXiv query error", ATOMS))
            authors = [a.findtext("a:name", "", ATOMS) for a in ent.findall("a:author",ATOMS)]
            abstract = ent.findtext("a:summary","",ATOMS)
            add(url, title, abstract, authors, url,
                "arXiv Corbet timing corpus" if survey else "arXiv Corbet" if corbet else "arXiv",
                by_source=not survey)
            identifier = url.split("/abs/")[-1]
            if re.fullmatch(r"(?:[a-z-]+(?:\.[A-Z]+)?/\d{7}|\d{4}\.\d{4,5})(?:v\d+)?", identifier):
                if survey or corbet or (re.search("period|orbit|timing|cycle|recurr|modulat|outburst", title,re.I)
                                        and any(a in norm(title+" "+abstract) for a in names)):
                    identifiers.append(identifier)
        total = int(root.findtext("o:totalResults",str(len(entries)),ATOMS))
        return dict(query=meta, retrieved=len(entries), total=total, truncated=total>len(entries))

    def cross(alias, corbet):
        text = alias + " orbital superorbital period recurrence"
        params = {"query.bibliographic":text, "rows":25}
        if corbet:
            params["query.author"] = "Corbet"
        rows, m = xref_rows(archive, params)
        for work in rows:
            title = (work.get("title") or [""])[0]
            authors = [i.get("family","") for i in work.get("author",[])]
            if corbet and not any(v.lower()=="corbet" for v in authors):
                continue
            doi = str(work.get("DOI") or "")
            add(doi, title, work.get("abstract",""), authors, "https://doi.org/"+doi,
                "Crossref Corbet" if corbet else "Crossref")
        return dict(query=m, retrieved=len(rows), ranked_limit=25)

    for corbet in [False, True]:
        query("arXiv Corbet" if corbet else "arXiv", lambda: arxiv(corbet))
        time.sleep(3.1)
        for alias in used[:3]:
            query(("Crossref Corbet " if corbet else "Crossref ")+alias, lambda: cross(alias, corbet))
    query("arXiv Corbet timing corpus", lambda: arxiv(True, survey=True))
    token = os.getenv("ADS_DEV_KEY")
    if token:
        def ads(corbet):
            aq = " OR ".join('object:"'+a.replace('"',"")+'"' for a in used)
            params = dict(q=f"({aq})" + (' author:"^Corbet"' if corbet else ""),
                          fl="bibcode,title,abstract,author", rows="200", sort="score desc")
            raw, m = archive.fetch("https://api.adsabs.harvard.edu/v1/search/query", params,
                                   {"Authorization":"Bearer "+token})
            data = json.loads(raw)["response"]
            for w in data["docs"]:
                bib = w["bibcode"]
                add(bib, w["title"][0], w.get("abstract",""), w.get("author",[]),
                    "https://ui.adsabs.harvard.edu/abs/"+bib, "ADS Corbet" if corbet else "ADS", True)
            return dict(query=m, retrieved=len(data["docs"]), total=data["numFound"],
                        truncated=data["numFound"]>len(data["docs"]))
        for corbet in [False, True]:
            query("ADS Corbet" if corbet else "ADS", lambda: ads(corbet))
    return dict(aliases_queried=used, queries=queries, papers=papers,
                manuscript_checks=manuscript_hits(archive, identifiers, aliases, target))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    archive = Archive(args.output)
    rows = [json.loads(line) for line in args.candidates.read_text().splitlines() if line.strip()]
    cats, tables = catalogues(archive) if rows else ({}, [])
    dump(args.output / "catalogues.json", cats)
    counts = Counter()
    with (args.output/"triage.jsonl").open("w") as ledger:
        for c in rows:
            aliases = uniq([c["source"]]+list(c.get("aliases") or []))
            try:
                obj, aliases = simbad(archive,c)
            except (ValueError, KeyError, TypeError) as exc:
                obj = dict(status="failed", error=str(exc))
            match = close_catalogues(tables, aliases, c)
            pubs = publications(archive, aliases, float(c["period_days"]))
            exact, harmonics = [], []
            for row in match:
                for p in row["periods"]:
                    near = p.get("near_ratios", [])
                    if near:
                        item = dict(family=row["family"], table=row["table"], matched=row["matched_identifiers"],
                                    arcmin=row["arcmin"], **p)
                        (exact if 1. in near else harmonics).append(item)
            incomplete = [k for k, v in cats.items() if v["status"]!="downloaded"]
            if obj.get("status")!="identified" or "bibliography_error" in obj or obj.get("bibliography_truncated"):
                incomplete.append("SIMBAD identity or bibliography")
            if not any("Corbet" in k and "error" not in v for k,v in pubs["queries"].items()):
                incomplete.append("Corbet search")
            if not any("Corbet" not in k and "error" not in v for k,v in pubs["queries"].items()):
                incomplete.append("literature search")
            state = ("catalogue period match for review" if exact else "possible harmonic for review" if harmonics
                     else "incomplete literature checks" if incomplete else "requires primary-paper review")
            detail = dict(candidate=c, simbad=obj, catalogues=match, bibliography=pubs,
                          similar_periods=exact, possible_harmonics=harmonics, incomplete=incomplete, status=state)
            filename = "candidate-" + hashlib.sha256((c["source"]+"|"+str(c["period_days"])).encode()).hexdigest()[:20]+".json"
            dump(args.output / filename, detail)
            ledger.write(json.dumps(dict(source=c["source"], period_days=c["period_days"], status=state,
                                         incomplete=incomplete, detail=filename,
                                         catalog_matches=len(exact), harmonic_matches=len(harmonics)))+"\n")
            ledger.flush()
            counts[state] += 1
    dump(args.output / "summary.json", dict(survivors=len(rows), triage=dict(counts),
         candidate_input_sha256=hashlib.sha256(args.candidates.read_bytes()).hexdigest(),
         scope="Catalogue matches and ranked abstract searches guide manual review; silence cannot certify novelty."))


if __name__ == "__main__":
    main()

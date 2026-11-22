"""Host and novelty checks for the slew search. Run inside the science job."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import threading
import subprocess
from html.parser import HTMLParser
from urllib.parse import urljoin, urlparse
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

XMATCH = ("https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync",
          "https://xmatch.cds.unistra.fr/xmatch/api/v1/sync")
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
NED = "https://ned.ipac.caltech.edu/tap/sync"
DATALAB = "https://datalab.noirlab.edu/query/query"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
SEARCH_RADIUS = 35.0
# Search terms apply only to papers associated with positional SIMBAD matches.
RELEVANT = re.compile(r"tidal.{0,12}disrupt|transient|changing.{1,4}(look|state)|"
                      r"(x.ray|nuclear).{0,50}(variab|flar|fading|outburst)|"
                      r"(variab|flar|fading|outburst).{0,50}(agn|nuclei|x.ray)|"
                      r"switch(ing|ed).{1,4}off", re.I)
STELLAR = re.compile(r"^(\*|Star|V\*|Fl\*|Er\*|RS\*|BY\*|Ro\*|TT\*|PM\*|"
                     r"YSO|CV\*|EB\*|Mira|PulsV\*|Flare\*|V\*\?|BYDra|RSCVn)$", re.I)
TRANSIENT = re.compile(r"^(TDE|TDe|SN|Nova|XB|XRB|HXB|LXB|CV\*|DN\*)\??$", re.I)
PUBLISHED_SEARCHES = {
    "TDE": "tidal disruption events",
    "changing-look": "changing-look AGN",
    "slew variability": "XMM slew variability",
    "slew transients": "XMM slew transient",
    "Strotjohann": "Strotjohann",
    "Sazonov": "Sazonov tidal disruption",
    "Saxton": "Saxton tidal disruption",
    "ROSAT flares": "ROSAT X-ray flares",
    "eROSITA transients": "eROSITA X-ray transients",
}


PRIOR_PAPERS = {
    "10.1093/mnras/stac714": ("Li2022", "J/MNRAS"),
    "10.1051/0004-6361/201628241": ("Strotjohann", "J/A+A"),
    "10.1111/j.1365-2966.2010.18024.x": ("Starling", "J/MNRAS"),
}


def chunks(values, n):
    for i in range(0, len(values), n):
        yield values[i:i + n]


def scalar(row, *names):
    keys = {str(k).strip().casefold(): v for k, v in row.items() if k is not None}
    for name in names:
        val = keys.get(name.casefold())
        if val is not None and str(val).strip() not in ("", "--", "nan"):
            return str(val).strip()
    return ""


def number(row, *names):
    try:
        value = float(scalar(row, *names))
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def sky(points):
    return SkyCoord([p["ra"] for p in points] * u.deg, [p["dec"] for p in points] * u.deg)


def clipped_failure(exc):
    return re.sub(r"\s+", " ", str(exc))[:420]


class Audit:
    def __init__(self):
        self.events = []
        self.mu = threading.Lock()

    def record(self, payload):
        with self.mu:
            self.events.append(payload)

    def get(self, name, url, *, method="GET", **kwargs):
        params = kwargs.get("data") or kwargs.get("params") or {}
        # Include the submitted coordinates/query or their digest; the original
        # input and code reproduce multipart uploads deterministically.
        signature = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
        try:
            response = requests.request(method, url, timeout=(10, 65), **kwargs)
            raw = response.content
            self.record({"service": name, "url": url, "request_sha256": signature,
                         "status": response.status_code, "bytes": len(raw),
                         "sha256": hashlib.sha256(raw).hexdigest()})
            response.raise_for_status()
            return raw
        except requests.RequestException as exc:
            self.record({"service": name, "url": url, "error": clipped_failure(exc)})
            raise


def csv_rows(raw, expected=()):
    src = raw.decode("utf-8-sig", "replace").lstrip()
    reader = csv.DictReader(io.StringIO(src))
    header = reader.fieldnames or []
    if src.startswith("<") or len(header) < 2 or any(k.lower() not in
            {h.lower() for h in header} for k in expected):
        raise ValueError("Unrecognized response: " + src[:220])
    return list(reader)


def tap(audit, name, url, query):
    return csv_rows(audit.get(name, url, method="POST", data={
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query,
        "MAXREC": "100001"}))


def select_points(rows, controls):
    points = []
    for i, row in enumerate(rows):
        ra, dec = float(row["ra"]), float(row["dec"])
        sig = number(row, "poserr_arcsec")
        gate = min(SEARCH_RADIUS, max(5, 2.45 * (sig if sig and sig > 0 else 8.0)))
        points.append(dict(k=i * 5, i=i, shift=0, ra=ra, dec=dec, gate=gate))
        if controls:
            c = SkyCoord(ra * u.deg, dec * u.deg)
            for j, pa in enumerate((0, 90, 180, 270), 1):
                s = c.directional_offset_by(pa * u.deg, 240 * u.arcsec)
                points.append(dict(k=i * 5 + j, i=i, shift=j, ra=s.ra.deg,
                                   dec=s.dec.deg, gate=gate))
    return points


def match_fields(input_points, rows, names, *, radius=SEARCH_RADIUS):
    """Return all in-cone rows with a great-circle separation in arcseconds."""
    results = defaultdict(list)
    fields, coords = [], []
    for row in rows:
        ra, dec = number(row, *names[0]), number(row, *names[1])
        if ra is None or dec is None or abs(dec) > 90:
            continue
        coords.append({"ra": ra, "dec": dec})
        fields.append(row)
    if rows and not coords:
        raise ValueError("Catalogue rows have no recognized decimal coordinates")
    if coords and input_points:
        ip, ic, d, _ = search_around_sky(sky(input_points), sky(coords), radius * u.arcsec)
        for i, j, sep in zip(ip, ic, d.arcsec):
            results[input_points[i]["k"]].append(dict(fields[j], sep_arcsec=round(sep, 3)))
    return results


def xmatch(audit, cat, points, fields):
    result, status = defaultdict(list), {}
    failures,reason = 0,""
    for part in chunks(points, 900):
        if failures >= 2:
            status.update({p["k"]:"service unavailable after repeated failures: "+reason for p in part})
            continue
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(("qkey", "input_ra", "input_dec"))
        for p in part:
            writer.writerow((p["k"], p["ra"], p["dec"]))
        try:
            data = {"request": "xmatch", "distMaxArcsec": str(SEARCH_RADIUS), "RESPONSEFORMAT": "csv",
                    "cat2": "vizier:" + cat, "colRA1": "input_ra", "colDec1": "input_dec"}
            data["RESPONSEFORMAT"] = "csv"
            for mirror in XMATCH:
                try:
                    raw = audit.get(cat, mirror, method="POST", data=data,
                                    files={"cat1": ("positions.csv", buf.getvalue(), "text/csv")})
                    rows = csv_rows(raw, ("qkey","angDist"))
                    break
                except Exception:
                    if mirror == XMATCH[-1]:
                        raise
            failures = 0
            for row in rows:
                k = int(row["qkey"])
                item = {field: scalar(row, *aliases) for field, aliases in fields.items()}
                item["sep_arcsec"] = float(row["angDist"])
                result[k].append(item)
            status.update({p["k"]: "ok" for p in part})
        except Exception as exc:
            reason = clipped_failure(exc)
            failures += 1
            status.update({p["k"]: reason for p in part})
    return result, status


def region_query(audit, name, points, make_query, names):
    result, status = defaultdict(list), {}
    failures,reason = 0,""
    for part in chunks(points, 35):
        if failures >= 2:
            status.update({p["k"]:"service unavailable after repeated failures: "+reason for p in part})
            continue
        try:
            rows = make_query(audit, name, part)
            if len(rows) >= 100000:
                raise ValueError("service row limit reached")
            matched = match_fields(part, rows, names)
            failures = 0
            for k, found in matched.items():
                result[k].extend(found)
            status.update({p["k"]: "ok" for p in part})
        except Exception as exc:
            reason = clipped_failure(exc)
            failures += 1
            status.update({p["k"]: reason for p in part})
    return result, status


def cones(part):
    return " OR ".join(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{p['ra']:.9f},"
                       f"{p['dec']:.9f},{SEARCH_RADIUS / 3600:.9f}))" for p in part)


def simbad_query(audit, name, part):
    return tap(audit, name, SIMBAD, "SELECT oid,main_id,otype,ra,dec FROM basic WHERE " + cones(part))


def ned_query(audit, name, part):
    rows = tap(audit, name, NED, "SELECT * FROM NEDTAP.objdir WHERE " + cones(part))
    return [{"ra": scalar(r, "ra"), "dec": scalar(r, "dec"),
             "preferred_name": scalar(r, "prefname", "objname", "name"),
             "type": scalar(r, "prefphytype", "pretype", "objtype", "type"), "z": scalar(r, "z", "redshift")}
            for r in rows]


def legacy_query(audit, name, part):
    # DR10 includes earlier northern imaging. Deduplication uses stable IDs.
    where = " OR ".join(f"q3c_radial_query(ra,dec,{p['ra']:.9f},{p['dec']:.9f},"
                         f"{SEARCH_RADIUS / 3600:.9f})" for p in part)
    query = "SELECT ls_id,type,ra,dec,flux_r,flux_z FROM ls_dr10.tractor WHERE " + where
    raw = audit.get(name, DATALAB, method="GET", params={"sql": query, "ofmt": "csv",
                    "async": "False"}, headers={"X-DL-AuthToken": "anonymous.0.0.anon_access"})
    return csv_rows(raw, ("ls_id", "ra", "dec"))


GAIA = {
    "id": ["Source", "source_id"], "ra": ["RAdeg", "RA_ICRS", "RAJ2000", "_RAJ2000"],
    "dec": ["DEdeg", "DE_ICRS", "DEJ2000", "_DEJ2000"], "plx": ["Plx"], "e_plx": ["e_Plx"],
    "pmra": ["pmRA"], "e_pmra": ["e_pmRA"], "pmdec": ["pmDE"], "e_pmdec": ["e_pmDE"],
    "ruwe": ["RUWE"], "gmag": ["Gmag"], "variable": ["VarFlag", "VARIABLE"]}

WISE = {"id": ["AllWISE", "WISE", "designation"], "ra": ["RAJ2000", "RAdeg"],
        "dec": ["DEJ2000", "DEdeg"], "w1": ["W1mag"], "w2": ["W2mag"], "ccf": ["ccf"]}

VSX = {"id": ["Name", "OID"], "type": ["Type"], "ra": ["RAJ2000"], "dec": ["DEJ2000"]}


def stellar_astrometry(row):
    # Gaia proper motions are the two projected components in mas/yr. Requiring
    # one significant component avoids treating a noisy vector length as Gaussian.
    plx, ep = number(row, "plx"), number(row, "e_plx")
    ruwe = number(row, "ruwe")
    if ruwe and ruwe > 1.4:
        return False
    if plx and ep and ep > 0 and plx / ep >= 5:
        return True
    for field in ("pmra", "pmdec"):
        v, e = number(row, field), number(row, "e_" + field)
        if v and e and e > 0 and abs(v) >= 5 and abs(v) / e >= 5:
            return True
    return False


def associate_refs(audit, simbad, points):
    ids = sorted({int(r["oid"]) for p in points for r in simbad.get(p["k"], [])
                  if r.get("oid") and float(r["sep_arcsec"]) <= p["gate"]})
    refs, checks = defaultdict(list), {}
    failures = 0
    for part in chunks(ids, 60):
        if failures >= 2:
            checks.update({k:"service unavailable after repeated bibliography failures" for k in part})
            continue
        keys = ",".join(map(str, part))
        try:
            query = ("SELECT h.oidref,r.bibcode,r.title FROM has_ref AS h JOIN ref AS r "
                     f"ON h.oidbibref=r.oidbib WHERE h.oidref IN ({keys})")
            rows = tap(audit, "SIMBAD-bibliography", SIMBAD, query)
            if len(rows) >= 100000:
                raise ValueError("bibliography truncated")
            for r in rows:
                if RELEVANT.search(r.get("title") or ""):
                    refs[int(r["oidref"])].append({"bibcode": r["bibcode"], "title": r["title"]})
            failures = 0
            checks.update({k: "ok" for k in part})
        except Exception as exc:
            failures += 1
            checks.update({k: clipped_failure(exc) for k in part})
    return refs, checks


def parse_coord(ra, dec):
    try:
        rd, dd = float(ra), float(dec)
        if 0 <= rd < 360 and abs(dd) <= 90:
            return rd, dd
    except (TypeError, ValueError):
        pass
    pos = SkyCoord(str(ra).replace(" ", ":"), str(dec).replace(" ", ":"),
                   unit=(u.hourangle, u.deg))
    return pos.ra.deg, pos.dec.deg


def tns(audit, points, previous):
    try:
        raw = audit.get("TNS", "https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objects.csv.zip")
        archive = zipfile.ZipFile(io.BytesIO(raw))
        file = next(s for s in archive.namelist() if s.lower().endswith(".csv"))
        # First line of the official export may contain a standalone timestamp.
        stream = io.TextIOWrapper(archive.open(file), encoding="utf-8-sig")
        header = next(stream)
        if "," not in header:
            header = next(stream)
        columns = next(csv.reader([header]))
        if "ra" not in {s.lower() for s in columns}:
            raise ValueError("TNS export has unrecognized columns")
        objs = []
        missed = 0
        for r in csv.DictReader(stream, fieldnames=columns):
            try:
                ra, dec = parse_coord(scalar(r, "RA"), scalar(r, "DEC", "declination"))
                objs.append({"ra": ra, "dec": dec, "name": scalar(r, "Name Prefix") + scalar(r, "Name"),
                             "type": scalar(r, "Obj. Type", "Type"), "date": scalar(r,"Discovery Date (UT)")})
            except (TypeError, ValueError):
                missed += 1
        if not objs:
            raise ValueError("TNS export returned zero usable objects")
        matched = match_fields(points, objs, (("ra",), ("dec",)))
        audit.record({"service":"TNS", "parsed":len(objs), "unparsed":missed})
        return matched, {p["k"]: "ok" if not missed else f"partial ({missed} unreadable)" for p in points}
    except Exception as exc:
        archive_error = clipped_failure(exc)
    # AT 2018cow validates the TNS coordinate convention, filters and returned
    # coordinate fields before any empty cone is counted as checked.
    def cone(ra, dec):
        c = SkyCoord(ra * u.deg, dec * u.deg)
        raw = audit.get("TNS-cone", "https://www.wis-tns.org/search", params={
            "ra": c.ra.to_string(unit=u.hourangle, sep=":", precision=3),
            "decl": c.dec.to_string(sep=":", alwayssign=True, precision=2),
            "radius": SEARCH_RADIUS, "coords_unit": "arcsec", "format": "csv", "num_page": 100})
        rows = csv_rows(raw)
        header = {s.strip().lower() for s in next(csv.reader([raw.decode("utf-8-sig").splitlines()[0]]))}
        if not {"ra", "dec"} <= header or "name" not in header or len(rows) >= 100:
            raise ValueError("TNS cone schema/truncation error")
        return rows

    result, status = defaultdict(list), {}
    try:
        control = cone(244.000925, 22.268008)
        if not any("2018cow" in scalar(r, "Name").replace(" ", "") for r in control):
            raise ValueError("TNS official search missed positive control AT 2018cow")
    except Exception as exc:
        return result, {p["k"]: archive_error + " | " + clipped_failure(exc) for p in points}

    def lookup(p):
        hits = []
        try:
            rows = cone(p["ra"], p["dec"])
            for r in rows:
                ra, dec = parse_coord(scalar(r, "RA"), scalar(r, "DEC", "declination"))
                center = SkyCoord(p["ra"] * u.deg, p["dec"] * u.deg)
                sep = center.separation(SkyCoord(ra*u.deg, dec*u.deg)).arcsec
                if sep > SEARCH_RADIUS + 0.2:
                    raise ValueError("TNS search returned an object outside requested cone")
                hits.append({"name": scalar(r,"Name"), "type":scalar(r,"Type","Obj. Type"),
                             "ra":ra,"dec":dec,"sep_arcsec":round(sep,3)})
            return hits, "ok", False
        except Exception as exc:
            blocked = isinstance(exc, requests.HTTPError) and exc.response.status_code in (401,403,429)
            return hits, archive_error + " | " + clipped_failure(exc), blocked

    coords = defaultdict(list)
    for p in points:
        # A successful bulk download above covers even the disqualified inputs.
        # On the per-target fallback, avoid hammering TNS about known X-ray binaries.
        if previous.get(p["k"],False):
            status[p["k"]] = "previously classified; TNS bulk unavailable"
        else:
            coords[(p["ra"],p["dec"])].append(p)
    with ThreadPoolExecutor(max_workers=2) as pool:
        blocked = ""
        for group in chunks(list(coords.values()), 2):
            if blocked:
                for members in group:
                    for q in members:
                        status[q["k"]] = "not attempted after service block: "+blocked
                continue
            for members, (hits, state, stop) in zip(group,pool.map(lookup,[arr[0] for arr in group])):
                for q in members:
                    result[q["k"]].extend(hits)
                    status[q["k"]] = state
                if stop:
                    blocked = state
    return result, status


def exact_papers(audit):
    records = {}
    for doi,(label,journal) in PRIOR_PAPERS.items():
        try:
            msg = json.loads(audit.get("Crossref-DOI", "https://api.crossref.org/works/"+doi))["message"]
            page = str(msg.get("page") or msg.get("article-number") or "").split("-")[0].strip()
            volume = str(msg.get("volume","")).strip()
            resource = journal+"/"+volume+"/"+page if page and volume else ""
            records[doi] = (msg,resource,label)
            audit.record({"service":"DOI-to-CDS", "doi":doi, "title":msg.get("title"),
                          "volume":volume,"first_page":page,"cds_lookup":resource})
        except Exception as exc:
            audit.record({"service":"DOI-to-CDS","doi":doi,"error":clipped_failure(exc)})
    return records


def publication_catalogs(audit,papers):
    from astroquery.vizier import VizierClass
    viz = VizierClass(columns=["*", "_RAJ2000", "_DEJ2000"], row_limit=20001)
    viz.TIMEOUT = 28
    catalogues = {}
    coverage = {}
    mirror = VizierClass(columns=["*", "_RAJ2000", "_DEJ2000"],row_limit=20001,
                         vizier_server="vizier.cfa.harvard.edu")
    mirror.TIMEOUT = 28
    for label, terms in PUBLISHED_SEARCHES.items():
        try:
            try:
                found = viz.find_catalogs(terms)
            except Exception:
                found = mirror.find_catalogs(terms)
            for name, resource in found.items():
                desc = str(resource.description or "")
                if label == "Strotjohann" and not re.search(r"slew|tidal|AGN|X-ray", desc, re.I):
                    continue
                catalogues.setdefault(name, set()).add(label)
            coverage[label] = {"catalogue_count":len(found), "status":"ok"}
        except Exception as exc:
            coverage[label] = {"status":clipped_failure(exc)}
    # Independent published sample with explicit archive identifier. The returned
    # table shape and coordinates still have to pass the same parsing checks.
    for _,resource,label in papers.values():
        if resource:
            catalogues.setdefault(resource,set()).add(label)
    audit.record({"service":"VizieR-publication-discovery", "search_terms":PUBLISHED_SEARCHES,
                  "coverage":coverage, "catalogues":sorted(catalogues)})
    return (viz,mirror), catalogues


def positional_columns(tab):
    keys = set(tab.colnames)
    for pair in (("_RAJ2000", "_DEJ2000"), ("RAJ2000", "DEJ2000"), ("RA_ICRS", "DE_ICRS"),
                 ("RAdeg", "DEdeg"), ("RA", "DE"), ("RA", "Dec")):
        if set(pair) <= keys:
            return pair
    coords = {}
    for name in keys:
        ucd = str(tab[name].meta.get("ucd", "")).lower()
        if ucd in ("pos.eq.ra;meta.main", "pos.eq.dec;meta.main"):
            coords[ucd.split(";")[0]] = name
    return coords.get("pos.eq.ra"), coords.get("pos.eq.dec")


def hash_table(tab):
    digest = hashlib.sha256()
    digest.update(json.dumps(tab.colnames).encode() + b"\n")
    for row in tab:
        digest.update(json.dumps([str(x) for x in row], ensure_ascii=True).encode() + b"\n")
    return digest.hexdigest()


class ArticleHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.parts = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "a" and a.get("href"):
            self.links.append(a["href"])
        if tag == "meta" and a.get("name") in ("citation_fulltext_html_url","citation_pdf_url"):
            self.links.append(a.get("content",""))

    def handle_data(self, data):
        self.parts.append(data)


IAU_NAME = re.compile(r"\b(XMMSL[12]?|[12]RXS|eRASS[0-9]?|[234]XMM|SRGe|SDSS|RX)"
                      r"\s*J?\s*(\d{6}(?:\.\d{1,3})?)\s*([+-])\s*(\d{6}(?:\.\d{1,2})?)",re.I)


def article_objects(audit, points, papers):
    """Find coordinate-coded objects in open primary papers; a hit flags review."""
    searches = [
        ("slew","The evolution of highly variable AGN in the XMM-Newton slew survey Strotjohann"),
        ("slew","XMM-Newton slew survey ROSAT X-ray transient sources tidal disruption Saxton"),
        ("slew","highly variable X-ray sources XMM-Newton Slew Survey Li"),
        ("TDE","First tidal disruption events discovered by SRG eROSITA Sazonov"),
    ]
    matches, searched, seendois = defaultdict(list),set(),set()
    allsearch = [("exact",doi) for doi in papers] + searches
    for kind, phrase in allsearch:
        try:
            if kind=="exact":
                records = [papers[phrase][0]]
            else:
                raw = audit.get("Crossref", "https://api.crossref.org/works",
                                params={"query.bibliographic":phrase,"rows":5,"select":"DOI,title,link,resource"})
                records = json.loads(raw)["message"]["items"]
            for rec in records:
                doi, title = rec["DOI"], rec.get("title",[""])[0]
                if doi in seendois:
                    continue
                if kind == "slew" and not ("slew" in title.lower() and re.search(
                        r"varia[bnt]|transient|disrupt|flare|outburst",title,re.I)):
                    continue
                if kind == "TDE" and not (re.search(r"tidal.{0,10}disrupt",title,re.I) and
                                           re.search(r"eROSITA|SRG",title,re.I)):
                    continue
                seendois.add(doi)
                urls = [rec.get("resource",{}).get("primary",{}).get("URL","")]
                urls.extend(x.get("URL","") for x in rec.get("link",[]))
                urls = [x for x in dict.fromkeys(urls) if x.startswith("https://")]
                visited, seenids = set(), {}
                for link in urls:
                    if link in visited:
                        continue
                    visited.add(link)
                    try:
                        page = audit.get(doi,link)
                        if page.startswith(b"%PDF"):
                            run = subprocess.run(["pdftotext","-layout","-","-"],input=page,
                                                 capture_output=True,timeout=45,check=True)
                            body = run.stdout.decode("utf-8","replace")
                        else:
                            doc = ArticleHTML()
                            doc.feed(page.decode("utf-8","replace"))
                            body = " ".join(doc.parts)
                            # Full-text and stand-alone table URLs from the
                            # publisher's own article, excluding citation links.
                            for href in doc.links:
                                stem = re.search(r"(aa[0-9]+[-_][0-9]+)",link,re.I)
                                if (stem and stem.group().lower() in href.lower() and
                                    re.search(r"(_T\d+\.html|/full_html/|/pdf/)",href,re.I)) :
                                    dest = urljoin(link,href).split("#")[0]
                                    if urlparse(dest).hostname==urlparse(link).hostname and dest not in urls:
                                        urls.append(dest)
                        body = body.replace("−","-").replace("–","-")
                        for m in IAU_NAME.finditer(body):
                            prefix, a, sign, d = m.groups()
                            try:
                                ra,dc = parse_coord(a[:2]+":"+a[2:4]+":"+a[4:],
                                                    sign+d[:2]+":"+d[2:4]+":"+d[4:])
                            except ValueError:
                                continue
                            alias = f"{prefix} J{a}{sign}{d}"
                            seenids[alias] = {"name":alias,"ra":ra,"dec":dc,"doi":doi,"title":title}
                    except Exception as exc:
                        audit.record({"service":"paper-full-text","doi":doi,"url":link,"error":clipped_failure(exc)})
                if seenids:
                    searched.add(kind)
                    mapped = match_fields(points,list(seenids.values()),(("ra",),("dec",)))
                    for k,hits in mapped.items():
                        matches[k].extend(hits)
                audit.record({"service":"paper-full-text","doi":doi,"title":title,
                              "distinct_coordinate_names":len(seenids),"urls_checked":sorted(visited)})
        except Exception as exc:
            audit.record({"service":"Crossref","query":phrase,"error":clipped_failure(exc)})
    return matches, searched


def publications(audit, points):
    result, inventory = defaultdict(list), []
    exact = exact_papers(audit)
    (viz,mirror), catalogues = publication_catalogs(audit,exact)
    for cat, terms in sorted(catalogues.items()):
        try:
            try:
                table_list = viz.get_catalogs(cat)
                if not len(table_list):
                    raise ValueError("empty primary resource")
            except Exception:
                table_list = mirror.get_catalogs(cat)
            if not len(table_list):
                raise ValueError("empty resource or service failure")
            for table_name in table_list.keys():
                tab = table_list[table_name]
                item = {"catalog":cat, "table":table_name, "terms":sorted(terms), "rows":len(tab)}
                inventory.append(item)
                if len(tab) >= 20001:
                    item["status"] = "truncated; excluded"
                    continue
                ra, dec = positional_columns(tab)
                if not ra or not dec:
                    item["status"] = "table lacks paired sky coordinates"
                    continue
                # Hash the exact table projection used for the crossmatch.
                item["rows_sha256"] = hash_table(tab)
                objects, failures = [], 0
                idfields = [x for x in tab.colnames if x.lower() in {
                    "name","objid","id","source","sourceid","tns","tnsname","simbad","iau",
                    "xmmsl1","xmmsl2","rosat","erosita","object","atname","galaxy","event"}]
                for idx, r in enumerate(tab):
                    try:
                        x,y = parse_coord(r[ra],r[dec])
                    except (ValueError,TypeError):
                        failures += 1
                        continue
                    ids = " | ".join(str(r[k]) for k in idfields if str(r[k]) != "--")
                    objects.append({"ra":x, "dec":y, "catalog":cat, "table":table_name,
                                    "record":idx+1, "name":ids, "search_terms":sorted(terms)})
                matched = match_fields(points, objects, (("ra",),("dec",)))
                for k,v in matched.items():
                    result[k].extend(v)
                item.update(status="ok", usable=len(objects), invalid_coordinates=failures)
        except Exception as exc:
            inventory.append({"catalog":cat, "terms":sorted(terms), "status":clipped_failure(exc)})
    audit.record({"service":"VizieR-publication-tables", "tables":inventory})
    # A TDE sample and a slew variability sample must both have searchable
    # positions. Other omissions are reported in the full coverage inventory.
    complete = {t for i in inventory if i.get("status")=="ok" and i.get("usable",0)
                and not i.get("invalid_coordinates") for t in i["terms"]}
    fulltext,article_kinds = article_objects(audit,points,exact)
    for k,hits in fulltext.items():
        for h in hits:
            result[k].append(dict(h,table="doi:"+h["doi"]))
    covered = ("TDE" in complete or "TDE" in article_kinds) and (
        bool(complete & {"Strotjohann", "slew variability"}) or "slew" in article_kinds)
    absent = []
    if not covered:
        absent.append("TDE/slew samples")
    if "Li2022" not in complete:
        absent.append("Li et al. 2022 source table (10.1093/mnras/stac714)")
    return result, {p["k"]:"ok" if not absent else "incomplete "+", ".join(absent) for p in points}



def compressed_json(path, records):
    with path.open("wb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as gz:
        for record in records:
            gz.write((json.dumps(record, sort_keys=True, default=str) + "\n").encode())


def write_table(path, rows, columns):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def is_galaxy(name):
    return bool(re.match(r"^(G|Galaxy|AGN|QSO|Sy1|Sy2|Seyfert|LINER|Bla|BLLac|EmG|GiP|GiG|rG|RadioG|SBG)\??$",name,re.I))


def extended(source):
    return str(source.get("type","")).strip().upper() in ("EXP","DEV","SER","REX","COMP")


def distance(first, second):
    try:
        a=SkyCoord(float(first["ra"])*u.deg,float(first["dec"])*u.deg)
        b=SkyCoord(float(second["ra"])*u.deg,float(second["dec"])*u.deg)
        return a.separation(b).arcsec
    except (ValueError,TypeError,KeyError):
        return math.inf


def host_summary(found):
    candidates = []
    for key, ident, typ in (("SIMBAD","main_id","otype"),("NED","preferred_name","type")):
        for r in found[key]:
            if is_galaxy(scalar(r,typ)):
                candidates.append((scalar(r,ident), scalar(r,typ), key, r))
    for r in found["Legacy"]:
        if extended(r):
            candidates.append((scalar(r,"ls_id"),scalar(r,"type"),"Legacy",r))
    if not candidates:
        for key in ("SIMBAD", "NED", "Gaia", "Legacy", "WISE"):
            fields = {"SIMBAD":("main_id","otype"),"NED":("preferred_name","type"),
                      "Gaia":("id",""),"Legacy":("ls_id","type"),"WISE":("id","")}
            if found[key]:
                r = found[key][0]
                f,t = fields[key]
                return (scalar(r,f), scalar(r,t), key, r["sep_arcsec"], False, "unresolved host class")
        return ("", "", "", "", False, "no optical/IR counterpart identified")
    name,typ,key,first = candidates[0]
    rivals = [r for _,_,_,r in candidates[1:] if distance(first,r)>2.0]
    # Unclassified/pointlike optical alternatives also matter.
    rivals.extend(r for r in found["Legacy"] + found["Gaia"] if distance(first,r)>2.0)
    # A secure host needs an independently identified nucleus/generic galaxy,
    # or an extended Legacy source independently coincident with AllWISE.
    independent = key in ("NED","SIMBAD") or any(distance(first,r)<=2 for r in found["WISE"])
    reason = "competing optical counterparts" if rivals else (
        "one extended source without independent identification" if not independent else "")
    return (name,typ,key,first["sep_arcsec"],not reason,reason)

MPE_BASE = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"
MPE_LISTS = {
    "HamStar-likely":"FreundS_DR1/HamStar_eRASS1_Main_Likely_Identifications_v1.1.fits.tgz",
    "HamStar-possible":"FreundS_DR1/HamStar_eRASS1_Main_Possible_Identifications_v1.1.fits.tgz",
    "BlazEr1":"HaemmerichS_DR1/BlazEr1_catalog.fits",
    "BLAZE":"HaemmerichS_DR1/BLAZE_catalog.fits",
    "HMXB":"KaltenbrunnerD_DR1/eRASS1_HMXB_LMC_v1.0.fits.tgz",
    "ULX":"WeberP_DR1/ulx_erass1_main.fits",
    "CV":"SchwopeA_DR1/allnewdr1cvs_fin_2cds.fits",
}


def ancillary(audit, points):
    import tarfile
    from astropy.io import fits
    result, states = defaultdict(list), {}
    if not points:
        return result, states
    centers = sky(points)
    for cat, path in MPE_LISTS.items():
        try:
            raw = audit.get(cat, MPE_BASE + path)
            if path.endswith(".tgz"):
                with tarfile.open(fileobj=io.BytesIO(raw),mode="r:*") as tf:
                    entry = next(m for m in tf if m.isfile() and m.name.lower().endswith((".fits", ".fit")))
                    raw = tf.extractfile(entry).read()
            with fits.open(io.BytesIO(raw),memmap=False) as hdus:
                tab = next(h.data for h in hdus if getattr(h,"columns",None) and h.data is not None)
                cols = {n.lower():n for n in tab.names}
                pairs = [("ra","dec"),("ra_deg","dec_deg"),("ra_x","dec_x"),("ra_xray","dec_xray"),
                         ("ra_corr","dec_corr"),("ra_erosita","dec_erosita"),("ra_erass1","dec_erass1"),
                         ("raj2000","dej2000"), ("ra_opt","dec_opt"), ("ra_icrs","de_icrs")]
                ra,dec = next((cols[a],cols[b]) for a,b in pairs if a in cols and b in cols)
                xx, yy = np.asarray(tab[ra], dtype=float).ravel(), np.asarray(tab[dec],dtype=float).ravel()
                mask = np.isfinite(xx) & np.isfinite(yy) & (np.abs(yy)<=90) & (xx>=0) & (xx<360)
                if not np.any(mask):
                    raise ValueError("zero usable positions")
                flat = np.flatnonzero(mask)
                catpos = SkyCoord(xx[mask]*u.deg,yy[mask]*u.deg)
                a,b,sep,_ = search_around_sky(centers,catpos,SEARCH_RADIUS*u.arcsec)
                names = [cols[n] for n in ("iauname","name","detuid","erass1","designation","source_id",
                                             "srcid","object","gaia_dr3","gaia_id","xid") if n in cols]
                for i,j,d in zip(a,b,sep.arcsec):
                    idx = flat[j]
                    ids = [str(np.asarray(tab[n]).ravel()[idx]).strip() for n in names]
                    result[points[i]["k"]].append({"catalog":cat, "id":" | ".join(ids), "row":int(idx+1),
                                                   "ra":xx[idx], "dec":yy[idx], "sep_arcsec":round(d,3)})
                states[cat] = "ok"
                audit.record({"service":cat,"ra_col":ra,"dec_col":dec,
                              "valid_rows":int(mask.sum()), "total_rows":len(xx)})
        except Exception as exc:
            states[cat] = clipped_failure(exc) or type(exc).__name__
    audit.record({"service":"MPE-ancillary-coverage", "coverage":states})
    return result, {p["k"]:"ok" if all(s=="ok" for s in states.values()) else
                             "unsearched: "+",".join(c for c,s in states.items() if s!="ok")
                    for p in points}


def vet_rows(rows, output_dir):
    """Return one verdict per input row. Write all positional alternatives and service outcomes."""
    rows = list(rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    columns = ["id","ra","dec","match_radius_arcsec","vet_status","host_id","host_type","host_catalog",
               "host_offset_arcsec","rationale","known_associations","prior_bibcodes",
               "shifted_host_controls","service_outcomes"]
    if not rows:
        write_table(out/"vetting.csv",[],columns)
        return []
    points = select_points(rows, controls=True)
    real = [p for p in points if not p["shift"]]
    audit = Audit()
    tasks = {
        "Gaia":(xmatch,(audit,"I/355/gaiadr3",points,GAIA)),
        "WISE":(xmatch,(audit,"II/328/allwise",points,WISE)),
        "VSX":(xmatch,(audit,"B/vsx/vsx",real,VSX)),
        "Legacy":(region_query,(audit,"Legacy",points,legacy_query,(("ra",),("dec",)))),
        "NED":(region_query,(audit,"NED",real,ned_query,(("ra",),("dec",)))),
        "SIMBAD":(region_query,(audit,"SIMBAD",real,simbad_query,(("ra",),("dec",)))),
        "papers":(publications,(audit,real)),
        "MPE":(ancillary,(audit,real)),
    }
    found, status = {}, {}
    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        futures = {k:executor.submit(fn,*args) for k,(fn,args) in tasks.items()}
        for k,fut in futures.items():
            try:
                found[k],status[k] = fut.result()
            except Exception as exc:
                found[k] = {}
                status[k] = {p["k"]:clipped_failure(exc) or type(exc).__name__ for p in points}
    previously_classified = {}
    for p in real:
        k = p["k"]
        previous_xray = any(r["sep_arcsec"] <= p["gate"] and TRANSIENT.fullmatch(r.get("otype",""))
                            for r in found["SIMBAD"].get(k,[]))
        previous_mpe = any(r["sep_arcsec"] <= p["gate"] and r["catalog"] in ("HamStar-likely","CV","HMXB","ULX")
                           for r in found["MPE"].get(k,[]))
        previously_classified[k] = previous_xray or previous_mpe
    found["TNS"],status["TNS"] = tns(audit,real,previously_classified)
    refs, ref_status = associate_refs(audit, found["SIMBAD"], real)
    output, evidence, controls = [], [], []
    for p in real:
        k = p["k"]
        source = rows[p["i"]]
        sid = str(source["id"])
        in_gate = {name:sorted((h for h in hits.get(k,[]) if h["sep_arcsec"]<=p["gate"]),
                               key=lambda h:h["sep_arcsec"]) for name,hits in found.items()}
        checks = {name:states.get(k,"not searched") for name,states in status.items()}
        sim_ids = [int(h["oid"]) for h in in_gate["SIMBAD"] if h.get("oid")]
        checks["SIMBAD_bibliography"] = ("ok" if checks["SIMBAD"]=="ok" and
                all(ref_status.get(n)=="ok" for n in sim_ids) else "incomplete")
        citations = [ref for oid in sim_ids for ref in refs.get(oid,[])]
        host,typ,cat,sep,secure,host_issue = host_summary(in_gate)
        known, uncertain = [], []
        xtype = scalar(source,"xmm_xtyp")
        xident = scalar(source,"xmm_counterpart","xmm_altident")
        xidsep = number(source,"xmm_iddist")
        if xident and xtype and xidsep is not None and 0 <= xidsep <= p["gate"]:
            if (re.search(r"\bstar\b|dwarf|stellar|cataclysmic|x.ray binary|flare star",xtype,re.I)
                    or STELLAR.fullmatch(xtype) or TRANSIENT.fullmatch(xtype)):
                if "?" in xtype:
                    uncertain.append("tentative slew counterpart "+xident+" "+xtype)
                else:
                    known.append("slew counterpart "+xident+" "+xtype)
        for h in in_gate["VSX"]:
            # VSX includes unconfirmed variables and miscellaneous sources.
            uncertain.append("VSX "+h["id"]+" "+h["type"])
        for h in in_gate["SIMBAD"]:
            if TRANSIENT.fullmatch(h["otype"]):
                known.append("SIMBAD "+h["main_id"]+" "+h["otype"])
            elif STELLAR.fullmatch(h["otype"]) or "*" in h["otype"]:
                uncertain.append("stellar SIMBAD "+h["main_id"]+" "+h["otype"])
        for h in in_gate["MPE"]:
            tag = h["catalog"]+" "+(h["id"] or f"row {h['row']}")
            if h["catalog"] in ("HamStar-likely","CV","HMXB","ULX"):
                known.append(tag)
            else:
                uncertain.append(tag)
        for h in in_gate["TNS"]:
            uncertain.append("TNS "+h["name"]+" "+h["type"])
        for h in in_gate["papers"]:
            uncertain.append(h["table"]+" "+(h["name"] or f"row {h['record']}"))
        for h in in_gate["Gaia"]:
            if stellar_astrometry(h):
                uncertain.append("significant Gaia stellar astrometry "+h["id"])
                if any(distance(h,s)<=2 for s in in_gate["VSX"]):
                    known.append("astrometric star coincident with VSX "+h["id"])
            if scalar(h, "variable").upper() == "VARIABLE":
                uncertain.append("Gaia variable "+h["id"])
        for ref in citations:
            uncertain.append("related published variability paper "+ref["bibcode"])

        ctl = {}
        for label,catalog,qualifies in (("Gaia","Gaia",bool),("Legacy","Legacy",bool),
                    ("WISE","WISE",bool),("Gaia_significant_astrometry","Gaia",stellar_astrometry),
                    ("Legacy_extended","Legacy",extended)):
            flags = [any(x["sep_arcsec"]<=p["gate"] and qualifies(x) for x in found[catalog].get(k+j,[]))
                     for j in range(1,5) if status[catalog].get(k+j)=="ok"]
            ctl[label] = [sum(flags),len(flags)]
            controls.append({"id":sid,"catalog":label,"shifted_with_counterpart":sum(flags),
                             "shifted_checked":len(flags),"real_with_counterpart":int(any(qualifies(x) for x in in_gate[catalog])),
                             "real_checked":int(checks[catalog]=="ok"),"match_radius_arcsec":round(p["gate"],2)})
        incomplete = [n for n,s in checks.items() if s!="ok" and n!="MPE"]
        if known:
            outcome,why = "known_or_stellar","positional association with a previously classified source"
        elif uncertain or incomplete or not secure:
            outcome = "inconclusive"
            why = " | ".join(filter(None,[host_issue, "prior association or stellar/variable alternative" if uncertain else "",
                                          "incomplete "+",".join(incomplete) if incomplete else ""]))
        else:
            outcome,why = "candidate","unique nonstellar host; no prior event in completed positional and bibliography checks"
        output.append(dict(id=sid,ra=p["ra"],dec=p["dec"],match_radius_arcsec=round(p["gate"],2),
                           vet_status=outcome,host_id=host,host_type=typ,host_catalog=cat,host_offset_arcsec=sep,
                           rationale=why,known_associations=" | ".join(dict.fromkeys(known+uncertain)),
                           prior_bibcodes=" | ".join(sorted({r["bibcode"] for r in citations})),
                           shifted_host_controls=json.dumps(ctl,sort_keys=True),
                           service_outcomes=json.dumps(checks,sort_keys=True)))
        evidence.append({"id":sid,"ra":p["ra"],"dec":p["dec"],"gate_arcsec":p["gate"],
                         "alternatives_within_35arcsec":{n:h.get(k,[]) for n,h in found.items()},
                         "slew_identification":{"id":xident,"class":xtype,"separation_arcsec":xidsep},
                         "bibliography_hits":citations,"checks":checks})
    write_table(out/"vetting.csv",output,columns)
    write_table(out/"host_shifted_controls.csv",controls,list(controls[0]))
    compressed_json(out/"counterpart_evidence.jsonl.gz",evidence)
    compressed_json(out/"vetting_provenance.jsonl.gz",audit.events)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    with open(args.input, newline="", encoding="utf-8") as inp:
        vet_rows(csv.DictReader(inp), args.output_dir)

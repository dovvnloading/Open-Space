"""Recover and audit the Neumann et al. Galactic HMXB catalogue on an Actions runner."""
import argparse
import csv
import gzip
import io
import json
import re
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

from astropy.table import Table
from bs4 import BeautifulSoup

from catalogue_archive import binary_tables, display
from literature import Archive, asu_tables, close_catalogues

VID = "J/A+A/677/A134"
DOI = "10.1051/0004-6361/202245728"
PAPER = "https://www.aanda.org/articles/aa/full_html/2023/09/aa45728-22/aa45728-22.html"
H = "https://vizier.cfa.harvard.edu/viz-bin/"
CDS = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
HEA = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"


def first_lines(path):
    return [json.loads(s) for s in path.read_text().splitlines() if s.strip()]


def source_columns(t):
    def physical(f):
        name = f["name"].lower()
        text = name + " " + f["desc"].lower()
        return (re.search(r"p[_\s]*orb|orb(?:ital)?[._\s]*per|orbital period", text) and
                not re.search(r"^([erufl]_|err|ref)|uncertainty|reference|error", name+" "+f["desc"].lower()))
    return [f for f in t["fields"] if physical(f)]


def accepted(tables):
    keep = [t for t in tables if len(t["rows"]) >= 100 and source_columns(t)]
    if not keep:
        raise ValueError("response lacks a Galactic-source-scale table with a documented orbital-period column")
    return keep


def tsv_tables(raw):
    text = raw.decode("utf-8", errors="replace")
    if VID not in text or "Neumann" not in text:
        raise ValueError("TSV lacks the exact catalogue/author identity")
    tables = []
    for section in re.split(r"(?m)^#Table\s+", text):
        data = [s for s in section.splitlines() if s.strip() and not s.startswith("#")]
        if len(data)<4 or "\t" not in data[0]:
            continue
        rows = list(csv.reader(data, delimiter="\t"))
        names = [v.strip() for v in rows[0]]
        units = [v.strip() for v in rows[1]]
        fields = []
        for i, k in enumerate(names):
            m = re.search(r"(?m)^#Column\s+" + re.escape(k) + r"\s+\([^)]*\)\s*(.*)$", section)
            desc = m[1].strip() if m else ""
            fields.append(dict(name=k, unit=units[i] if i<len(units) else "", desc=desc,
                               ucd=desc.split("[ucd=")[-1].rstrip("]") if "[ucd=" in desc else ""))
        good = [{n:v.strip() for n,v in zip(names, row)} for row in rows[3:] if len(row)==len(names)]
        name = re.search(r"(?m)^#Name:\s*(.+)", section)
        tables.append(dict(name=name[1].strip() if name else VID+"/hmxbcat", fields=fields, rows=good))
    return tables


def tap_request(arc, endpoint, query, fmt):
    return arc.fetch(endpoint, {"REQUEST":"doQuery","LANG":"ADQL","MAXREC":"20000","FORMAT":fmt,"QUERY":query})


def registry(arc, log):
    for host in (H, "https://vizier.iucaa.in/viz-bin/", "https://vizier.cds.unistra.fr/viz-bin/"):
        for fmt, name in (("asu-xml",VID+"/hmxbcat"), ("asu-tsv",VID)):
            try:
                raw, receipt = arc.fetch(host+fmt, {"-out.max":"20000","-source":name,"-out.add":"_RAJ,_DEJ","-out.all":"1"})
                if VID.encode() not in raw[:200000] or b"Neumann" not in raw[:200000]:
                    raise ValueError("response lacks the exact Neumann catalogue identity")
                tables = accepted(asu_tables(raw) if fmt=="asu-xml" else tsv_tables(raw))
                return tables, {"provider":"VizieR", "query":receipt, "catalogue":VID}
            except Exception as e:
                log.append({"route":host+fmt, "table":name,"error":str(e)[:500]})
    raise ValueError("no ASU mirror supplied the published catalogue")


def tap_catalogues(arc, log):
    for endpoint in (CDS, HEA):
        if endpoint==CDS:
            where = "table_name LIKE '%677/A134%'"
        else:
            where = ("table_name LIKE '%hmxb%' OR table_name LIKE '%xrbcat%' "
                     "OR description LIKE '%Neumann%' OR description LIKE '%XRBcats%'")
        try:
            raw, search = tap_request(arc,endpoint,"SELECT table_name,description FROM TAP_SCHEMA.tables WHERE "+where,"csv")
            tables = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
            if not tables:
                raise ValueError("no catalogue identity in registry response")
            for record in tables:
                name = record["table_name"].strip('"')
                if not re.fullmatch(r"[A-Za-z0-9_+/ .\-]+",name):
                    continue
                identity = record.get("description","")
                source = {"registry":search,"table":name,"description":identity}
                if endpoint == CDS and VID.lower() not in name.lower():
                    continue
                if endpoint == HEA:
                    if not re.search("hmxb|xrbcat|high.mass",name+" "+identity,re.I):
                        continue
                    if not (DOI in identity or ("Neumann" in identity and re.search("high.mass",identity,re.I))):
                        try:
                            slug = name.split(".")[-1]
                            page, meta = arc.fetch("https://heasarc.gsfc.nasa.gov/W3Browse/all/"+slug.lower()+".html")
                            doc = BeautifulSoup(page,"html.parser").get_text(" ",strip=True)
                            if not ((DOI in doc or "2023A&A...677A.134N" in doc or
                                     ("Neumann" in doc and re.search(r"677.{0,30}A134",doc))) and re.search("high.mass",doc,re.I)):
                                continue
                            source["identity_page"] = meta
                        except Exception as e:
                            log.append({"route":name, "error":str(e)[:300]})
                            continue
                try:
                    adql_name = '"'+name+'"' if endpoint == CDS else name
                    blob, receipt = tap_request(arc,endpoint,"SELECT * FROM "+adql_name,"votable")
                    parsed = accepted(binary_tables(blob))
                    for t in parsed:
                        t["name"] = name
                    return parsed, {"provider":"VizieR TAP" if endpoint==CDS else "HEASARC TAP",
                                    "query":receipt, **source}
                except Exception as e:
                    log.append({"route":name, "error":str(e)[:350]})
        except Exception as e:
            log.append({"route":endpoint, "error":str(e)[:400]})
    raise ValueError("no identity-verified full TAP source table")


def readable_table(raw, path):
    if raw.startswith(b"\x1f\x8b"):
        raw = gzip.decompress(raw)
        path = path.removesuffix(".gz")
    if raw.startswith(b"SIMPLE"):
        t = Table.read(io.BytesIO(raw), format="fits")
    elif path.lower().endswith((".vot",".votable",".xml")):
        return accepted(binary_tables(raw))
    else:
        txt = raw.decode("utf-8-sig")
        delimiter = "\t" if path.lower().endswith(".tsv") else None
        opts = {"delimiter":delimiter} if delimiter else {}
        t = Table.read(txt.splitlines(), format="ascii.ecsv" if path.lower().endswith(".ecsv") else "ascii", **opts)
    fields = []
    for k in t.colnames:
        unit = str(t[k].unit or "")
        if not unit:
            unit = next((u for pat,u in ((r"\b(?:days?|d)\b","d"), (r"\b(?:hours?|h)\b","h"), (r"\b(?:seconds?|s)\b","s"))
                         if re.search(pat,k.lower().replace("_"," "))), "")
        fields.append(dict(name=k, unit=unit, desc=str(t[k].description or k), ucd=""))
    return accepted([dict(name=path, fields=fields,rows=[{k:display(r[k]) for k in t.colnames} for r in t])])


def author_repo(arc, owner, repo, origin, errors):
    base = "https://api.github.com/repos/"+owner+"/"+repo
    metadata, repo_receipt = arc.fetch(base)
    branch = json.loads(metadata)["default_branch"]
    raw, commit_receipt = arc.fetch(base+"/commits/"+quote(branch,safe=""))
    commit = json.loads(raw)["sha"]
    content, tree_receipt = arc.fetch(base+"/git/trees/"+commit, {"recursive":"1"})
    files = [s for s in json.loads(content)["tree"] if s["type"]=="blob" and s.get("size",0)<20000000]
    immutable = f"https://raw.githubusercontent.com/{owner}/{repo}/{commit}/"
    readme = next((s for s in files if s["path"].lower() in ("readme.md","readme.rst","readme.txt")), None)
    if not readme:
        raise ValueError("repository lacks identity documentation")
    rb, rm = arc.fetch(immutable+readme["path"])
    readme_text = rb.decode("utf8",errors="replace")
    if not (DOI in readme_text or (re.search("Neumann",readme_text,re.I) and re.search("high.mass|HMXB",readme_text,re.I))):
        raise ValueError("repository README does not attribute this HMXB catalogue to Neumann")
    data = [s["path"] for s in files if re.search(r"\.(fits|csv|tsv|ecsv|vot|xml|dat)(\.gz)?$",s["path"],re.I)
            and "lmxb" not in s["path"].lower()]
    ordered = sorted(data,key=lambda x:(not bool(re.search("hmxb|high.mass",x,re.I)),not bool(re.search("catalog|source",x,re.I)),x))
    for f in ordered[:12]:
        try:
            b, meta = arc.fetch(immutable+quote(f))
            t = readable_table(b,f)
            return t,dict(provider="author-maintained repository, pinned revision",query=meta,commit=commit,
                          linked_from=origin,repository=repo_receipt,branch_query=commit_receipt,tree=tree_receipt,readme=rm)
        except Exception as e:
            errors.append({"route":f,"error":str(e)[:320]})
    raise ValueError("no source-scale identifiable machine-readable table in linked repository")


def author(arc, errors):
    # Only follow repositories linked from the actual paper or its full-text
    # preprint. Identity must additionally match the repository documentation.
    pending = [PAPER,"https://doi.org/"+DOI]
    linked_sites = set()
    try:
        query = 'ti:"XRBcats" AND au:Neumann'
        import xml.etree.ElementTree as ET
        z,_ = arc.fetch("https://export.arxiv.org/api/query",{"search_query":query,"max_results":"5"})
        for e in ET.fromstring(z).findall("{*}entry"):
            title = e.findtext("{*}title", "")
            if re.search("high.mass",title,re.I):
                ident=e.findtext("{*}id","").split("/abs/")[-1]
                if re.fullmatch(r"\d{4}\.\d+(?:v\d+)?",ident):
                    pending += ["https://arxiv.org/html/"+ident, "https://ar5iv.labs.arxiv.org/html/"+ident.removesuffix("v1")]
    except Exception as e:
        errors.append({"route":"arXiv paper lookup","error":str(e)[:320]})
    done=set()
    while pending and len(done)<12:
        url=pending.pop(0)
        if url in done:
            continue
        done.add(url)
        try:
            raw, receipt=arc.fetch(url)
            doc=BeautifulSoup(raw,"html.parser")
            body = doc.get_text(" ",strip=True)
            identity = ("Neumann" in body or DOI in raw.decode("utf8",errors="replace")) and re.search("high.mass|hmxb",body,re.I)
            if not identity:
                continue
            for a in doc.select("a[href]"):
                dest=urljoin(url,a["href"])
                m=re.match(r"https?://github\.com/([\w.-]+)/([\w.-]+)",dest)
                if m and m[2].lower() not in ("astropy","numpy","pandas","pyvo"):
                    try:
                        return author_repo(arc,m[1],m[2].removesuffix(".git"),receipt,errors)
                    except Exception as e:
                        errors.append({"route":dest,"error":str(e)[:320]})
                if "aanda.org/" in dest and "full_html/" in dest and not urlparse(dest).fragment and "45728" in dest:
                    pending.append(dest)
                if url not in linked_sites and re.search(r"(xrbcats|xray-binaries|xrbcat|hmxbcat)",dest,re.I) and not re.search(r"\.(?:gz|dat|csv|fits|pdf|vot|xml|txt)$",dest,re.I):
                    host = urlparse(dest).hostname or ""
                    if host not in {"github.com", "ui.adsabs.harvard.edu", "vizier.cds.unistra.fr"} and dest not in done:
                        pending.append(dest)
                        linked_sites.add(dest)
                if (url in linked_sites or re.search(r"hmxb|high.mass|xrbcats",dest,re.I)) and "lmxb" not in dest.lower() and re.search(r"\.(fits|csv|tsv|ecsv|vot)(\.gz)?$",dest,re.I):
                    try:
                        b,meta=arc.fetch(dest)
                        return readable_table(b,dest),dict(provider="published article data link",query=meta,linked_from=receipt)
                    except Exception as e:
                        errors.append({"route":dest,"error":str(e)[:320]})
        except Exception as e:
            errors.append({"route":url,"error":str(e)[:320]})
    raise ValueError("no identity-verified author table located via the paper")


def canonical(table):
    # Preserve every original cell. Supply semantic aliases for provider-specific
    # RA/DEC and source-ID labels so the same survey matcher handles all mirrors.
    for f in table["fields"].copy():
        old = f["name"]
        key = re.sub("[^a-z0-9]","",old.lower())
        new = ({"ra":"RAJ2000","raj2000":"RAJ2000","raicrs":"RA_ICRS","ra2000":"RAJ2000",
                "dec":"DEJ2000","decj2000":"DEJ2000","de2000":"DEJ2000","deicrs":"DE_ICRS",
                "dej2000":"DEJ2000","mainid":"Source_name","popularid":"Other_name",
                "simbadid":"SIMBAD_name","srcname":"Source_name","otherid":"Other_name"}.get(key))
        if new and new not in {x["name"] for x in table["fields"]}:
            table["fields"].append(dict(f, name=new))
            for r in table["rows"]:
                r[new] = r.get(old,"")


def run(evidence,out):
    arc=Archive(out/"receipts")
    unresolved=[]
    tables=[]
    access={}
    for route in (registry,tap_catalogues,author):
        try:
            tables,access=route(arc,unresolved)
            break
        except Exception as exc:
            unresolved.append({"route":route.__name__,"error":str(exc)[:500]})
    for t in tables:
        canonical(t)
    cands=first_lines(evidence/"survivors.jsonl")
    simids={}
    triage=evidence/"literature"/"triage.jsonl"
    if triage.is_file():
        for s in first_lines(triage):
            d=json.loads((triage.parent/s["detail"]).read_text())
            simids[d["candidate"]["source_id"]]=d.get("simbad",{}).get("aliases",[])
    rows=[]
    for c in cands:
        aliases = list(dict.fromkeys([c["source"],*c.get("aliases",[]),*simids.get(c["source_id"],[])]))
        found=close_catalogues([("Neumann HMXB",VID,t) for t in tables],aliases,c)
        rows.append(dict(source_id=c["source_id"],source=c["source"],candidate_days=c["period_days"],
                         ra_deg=c.get("ra"),dec_deg=c.get("dec"),identifiers=aliases,rows=found,
                         status=("catalogue_unavailable" if not tables else "no_match" if not found else
                                 "one_match" if len(found)==1 else "multiple_row_matches")))
    with (out/"candidates.jsonl").open("w") as f:
        for r in rows:
            f.write(json.dumps(r,allow_nan=False)+"\n")
    meta=dict(citation="Neumann et al. 2023, A&A 677, A134",doi=DOI,published_scope="Galactic HMXBs",
              access=access,status="retrieved" if tables else "unavailable",attempts=unresolved,
              n_published_sources=169,
              exact_published_row_count_present=any(len(t["rows"])==169 for t in tables),
              tables=[dict(name=t["name"],source_rows=len(t["rows"]),orbital_fields=source_columns(t)) for t in tables],
              candidates=len(rows),with_matches=sum(bool(r["rows"]) for r in rows))
    (out/"summary.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    print("Neumann HMXB",meta["status"],"source rows",meta["tables"],"matches",meta["with_matches"],"/",len(rows),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    run(a.evidence,a.output)

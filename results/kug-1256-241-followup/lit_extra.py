"""Resolve primary optical papers and check published changing-look source lists by position."""
import csv
import io
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urljoin
from astropy.io.votable import parse
from bs4 import BeautifulSoup
from pypdf import PdfReader
from common import OUT, HOST_RA, HOST_DEC, fetch, status

notes = []


def lookup_catalogs():
    url = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
    query = ("SELECT table_name,description FROM TAP_SCHEMA.tables WHERE "
             "description LIKE '%hanging-look%' OR description LIKE '%hanging look%' "
             "OR description LIKE '%hanging-state%'")
    try:
        raw = fetch("CL catalog directory",url,dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query),timeout=(6,22))
    except Exception:
        # The alias-based search and published-table test remain usable without the global schema index.
        raw = b"table_name,description\n"
    header=raw.decode(errors="replace").lower().splitlines()[0]
    if "table_name" not in header or "description" not in header:
        raise ValueError("Unrecognized TAP metadata response")
    items = list(csv.DictReader(io.StringIO(raw.decode())))
    catalogs = sorted({r["table_name"].strip('"').rsplit("/",1)[0] for r in items
                       if "/" in r["table_name"] and not re.search("white dwarf|binary|pulsar", r["description"],re.I)})
    # These two discovery catalogues give an independent route when the schema service is down.
    catalogs = sorted(set(catalogs) | {"J/ApJ/862/109","J/MNRAS/457/389"})
    def lookup(cat):
        try:
            args={"-source":cat,"-c":f"{HOST_RA} {HOST_DEC}","-c.rs":15,"-out.all":1,"-out.max":50,"-out.add":"_r"}
            raw=fetch("CL list "+cat,"https://vizier.cfa.harvard.edu/viz-bin/votable",args,timeout=20)
            if re.search(rb'<INFO[^>]*(?:name="Error"|value="ERROR")',raw,re.I):
                raise ValueError("VizieR service error")
            v=parse(io.BytesIO(raw))
            tables=list(v.iter_tables())
            n=sum(len(t.array) for t in tables)
            rows=[{k:str(a) for k,a in zip(r.colnames,r)} for t in tables for r in t.to_table(use_names_over_ids=True)]
            notes.append(dict(kind="changing-look catalogue",catalog=cat,entries=n,rows=rows,
                              state="rows" if n else "empty" if tables or b"No astronomical object found" in raw else "no recognized table"))
        except Exception as exc:
            notes.append(dict(kind="changing-look catalogue",catalog=cat,error=str(exc)[:240]))
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lookup,catalogs))
    status("Published changing-look catalogues","searched" if catalogs else "unavailable",
           f"{len(catalogs)} VizieR catalogues identified from their table descriptions and queried at the nucleus; scope and responses in literature-extra.json.")


def primary():
    papers = [("2001MNRAS.327..673V", "10.1046/j.1365-8711.2001.04763.x"),
              ("2008A&A...484..897K","10.1051/0004-6361:20078552"),
              ("2023ApJS..267...19I","10.3847/1538-4365/acd783")]
    def getpaper(pair):
        code,doi=pair
        urls=["https://articles.adsabs.harvard.edu/pdf/"+code]
        if code.startswith("2001"):
            urls.insert(0,"https://arxiv.org/pdf/astro-ph/0106311")
        try:
            raw=fetch("Crossref-links-"+code,"https://api.crossref.org/works/"+doi,timeout=15)
            meta=json.loads(raw)["message"]
            for v in meta.get("link",[]):
                if "pdf" in (v.get("content-type","")+v["URL"]).lower():
                    urls.insert(0,v["URL"])
            landing=meta.get("resource",{}).get("primary",{}).get("URL","")
            if landing:
                raw=fetch("publisher-"+code,landing,timeout=25)
                page=BeautifulSoup(raw,"html.parser")
                tag=page.select_one('meta[name="citation_pdf_url"]')
                if tag and tag.get("content"):
                    urls.insert(0,urljoin(landing,tag["content"]))
        except Exception as exc:
            notes.append(dict(kind="publisher metadata",bibcode=code,error=str(exc)[:180]))
        for url in dict.fromkeys(urls):
            try:
                raw=fetch("full text "+code,url,timeout=25)
                if not raw.startswith(b"%PDF"):
                    continue
                with TemporaryDirectory() as d:
                    p=Path(d)/"article.pdf";p.write_bytes(raw)
                    try:
                        s=subprocess.check_output(["pdftotext","-layout",str(p),"-"],timeout=30).decode(errors="replace")
                    except FileNotFoundError:
                        s="\n".join(page.extract_text(extraction_mode="layout") or "" for page in PdfReader(str(p)).pages)
                if len(s.strip())<1000:
                    continue
                hit=re.compile(r'1258\s*\.?\s*8\s*\+\s*2355|1256\s*\+\s*241|125851[.\d]+\s*\+\s*2355|RBS\s+1195',re.I)
                lines=[l.strip() for l in s.splitlines() if hit.search(l)]
                notes.append(dict(kind="historical paper",bibcode=code,doi=doi,url=url,target_rows=lines[:30]))
                status("Full text "+code,"retrieved", f"{len(lines)} exact target-name lines in machine-readable text")
                return
            except Exception as exc:
                notes.append(dict(kind="paper access",bibcode=code,error=str(exc)[:180]))
        status("Full text "+code,"unavailable","No readable full-text PDF from publisher or ADS")
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(getpaper,papers))


def literature_tables():
    # Published catalogues at the CfA mirror. Save each row with the table name and column definitions.
    catalogs=["J/MNRAS/327/673","J/A+A/484/897","J/ApJS/267/19","J/AJ/123/578","J/AJ/119/2540","J/A+A/416/35"]
    def one(c):
        try:
            data=fetch("published table CfA "+c, "https://vizier.cfa.harvard.edu/viz-bin/votable",
                       {"-source":c,"-c":f"{HOST_RA} {HOST_DEC}","-c.rs":30,"-out.max":50,"-out.all":"","-out.add":"_r"},timeout=18)
            if re.search(rb'<INFO[^>]*(?:name="Error"|value="ERROR")',data,re.I):
                raise ValueError("VizieR query error")
            entries=[]
            for t in parse(io.BytesIO(data)).iter_tables():
                rows=[{k:str(a) for k,a in zip(r.colnames,r)} for r in t.to_table(use_names_over_ids=True)]
                entries.append(dict(table=t.name,columns={f.name:[str(f.unit or ""),f.description] for f in t.fields},rows=rows))
            notes.append(dict(kind="targeted published data",catalog=c,tables=entries))
            status("CfA VizieR "+c,"retrieved" if any(t['rows'] for t in entries) else "empty_or_unrecognized",
                   f"{sum(len(t['rows']) for t in entries)} rows within 30 arcsec")
        except Exception as exc:
            status("CfA VizieR "+c,"unavailable",str(exc)[:240])
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(one,catalogs))


def notices():
    # Object-linked references, including ATels, are available from the public SIMBAD bibcode index.
    # A same-interface positive control distinguishes a working index from missing metadata.
    def linked(key,oid):
        query=("SELECT r.bibcode,r.title,r.abstract FROM has_ref AS h JOIN ref AS r "
               f"ON h.oidbibref=r.oidbib WHERE h.oidref={oid}")
        raw=fetch(key,"https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
                  dict(REQUEST="doQuery",LANG="ADQL",FORMAT="votable",QUERY=query),timeout=22)
        return list(parse(io.BytesIO(raw)).get_first_table().to_table(use_names_over_ids=True))
    control=linked("ATel bibliography control",8399363)
    if not any("ATel.3811" in str(r["bibcode"]) for r in control):
        raise ValueError("SIMBAD comparison source did not recover ATel 3811")
    rows=linked("nuclear object bibliography",2047449)
    from literature import ALIASES, TERM
    res=[]
    for r in rows:
        title,ab=str(r["title"]),str(r["abstract"])
        res.append(dict(bibcode=str(r["bibcode"]),title=title,
                        names=[a for a in ALIASES if a.casefold() in (ab+" "+title).casefold()],
                        transition_terms=TERM.findall(ab+" "+title)))
    notes.append(dict(kind="object-linked notices and abstracts",rows=res))
    status("Object-linked ATels and abstracts","retrieved",
           f"{len(rows)} linked references, {sum('ATel' in r['bibcode'] for r in res)} ATel bibcodes; known comparison ATel recovered.")


def ned_official():
    from urllib.parse import quote
    codes=set()
    for obj in ("KUG 1256+241","2MASX J12585149+2355268"):
        for host in ("ned.ipac.caltech.edu","nedwww.ipac.caltech.edu"):
            try:
                params={"objname":obj,"search_type":"Ref_id","extend":"no","ref_extend":"no",
                        "of":"xml_main","begin_year":"1800","end_year":"2026"}
                raw=fetch("NED legacy references "+obj+" "+host,
                          "https://"+host+"/cgi-bin/nph-refsearch",params,timeout=(6,15))
                soup=BeautifulSoup(raw,"html.parser")
                cells=[t.get_text(" ",strip=True) for t in soup.find_all(["td","a"])]
                codes.update(v for v in cells if re.fullmatch(r"(?:19|20)\d{2}[A-Za-z&.][A-Za-z&.0-9]{13}[A-Za-z.]",v))
                if codes:
                    break
            except Exception:
                continue
        if codes:
            break
    notes.append(dict(kind="NED returned official bibliographic identifiers",codes=sorted(codes)))
    status("NED reference lookup by historical names","retrieved" if codes else "unavailable",
           f"{len(codes)} validated bibcodes in the returned NED reference tables.")


def run():
    for name,func in [("changing-look source lists",lookup_catalogs),("target optical papers",primary),("published tables",literature_tables),("object-linked ATels",notices),("NED by alias",ned_official)]:
        try:
            func()
        except Exception as exc:
            status(name,"unavailable",str(exc)[:300])
    (OUT/"literature-extra.json").write_text(json.dumps(notes,ensure_ascii=False,indent=2)+"\n")

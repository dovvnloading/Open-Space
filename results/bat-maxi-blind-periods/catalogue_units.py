"""Validate the period columns used in the LMXB literature audit, on Actions."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import parse_qs, urlparse, urljoin, urlencode

from bs4 import BeautifulSoup
import numpy as np
from literature import Archive, asu_tables, norm, paper_text, sky_distance

FORTIN = "J/A+A/684/A124"
AVAKYAN = "J/A+A/675/A199"
FIELDS = {FORTIN:"Period", AVAKYAN:"Porb"}
DOIS = {FORTIN:"10.1051/0004-6361/202347908", AVAKYAN:"10.1051/0004-6361/202346522"}
FROZEN_DIGESTS = {FORTIN:"c56a1d87a45b291cd7475663d38b89575e1c74b186c9e1946337fa65d0f25b2b",
                  AVAKYAN:"6055f0b07f772a9c35b7c071a8aa7b94bd67ca6dd39a9a08d12c0e75877eb78d"}
RADIUS = 60.
WIDTH = .05
ATOM = "{http://www.w3.org/2005/Atom}"


def pretty(path, val):
    path.write_text(json.dumps(val, indent=2, ensure_ascii=False, allow_nan=False)+"\n")


def tables_for(raw, ident, ids, req):
    root=ET.fromstring(raw)
    for k in root.iter():
        k.tag=k.tag.split("}")[-1]
    for xt in root.iter("TABLE"):
        if xt.get("name","").strip('"')!=ident+"/lmxbcat":
            continue
        fs={f.get("name"):f for f in xt.findall("FIELD")}
        field=FIELDS[ident]
        if not set(ids+[field])<=fs.keys():
            continue
        col=fs[field]
        desc=col.findtext("DESCRIPTION") or ""
        if "time.period" not in col.get("ucd","") or "orbital period" not in desc.lower():
            raise ValueError("field does not describe an orbital period")
        block=xt.find("DATA/CSV")
        if block is not None and block.text:
            # Deserialize only this validated table. These original receipts
            # use VizieR's three-line-header, tab-delimited Astrores extension.
            lines=block.text.lstrip("\r\n").splitlines()
            headings=[v.strip() for v in lines[0].split("\t")]
            if not set(ids+[field])<=set(headings):
                raise ValueError("archive CSV lacks the required named fields")
            start=int(block.get("headlines","3"))
            rows=[dict(zip(headings,(v.strip() for v in line.split("\t")))) for line in lines[start:] if line.strip()]
        else:
            rows=next(t["rows"] for t in asu_tables(raw) if t["name"]==xt.get("name"))
        if not rows:
            raise ValueError("empty source table")
        return rows,dict(table=xt.get("name"),doi=DOIS[ident],title=xt.findtext("DESCRIPTION"),
                         query=req,field_attributes=col.attrib,description=desc,
                         identity_fields={k:dict(attributes=fs[k].attrib,
                         description=fs[k].findtext("DESCRIPTION")) for k in ids},n_catalogue_rows=len(rows))
    raise ValueError("expected table or named fields absent")


def catalogue(directory, ident, ids, arc):
    errors=[]
    seen=set()
    params={"-source":ident,"-out.all":"1","-out.max":"50000","-out.add":"_RAJ,_DEJ","-out.form":"VOTable"}
    host="https://vizier.cfa.harvard.edu/viz-bin/asu-xml"
    # The first-run SHA keys also work if a resumed job wrote a fresh query log.
    pin={"sha256":FROZEN_DIGESTS[ident],"url":host+"?"+urlencode(params),"routing":"verified first-run digest"}
    path=directory/"queries.jsonl"
    records=[json.loads(s) for s in path.read_text().splitlines() if s.strip()] if path.exists() else []
    # Preserve the full recorded request/date for the pinned response where available.
    pin=next((r for r in records if r.get("sha256")==pin["sha256"]),pin)
    for req in [pin]+records:
        digest=req.get("sha256")
        if not digest or digest in seen or ident not in parse_qs(urlparse(req["url"]).query).get("-source",[]):
            continue
        seen.add(digest)
        try:
            raw=gzip.decompress((directory/"responses"/(digest+".gz")).read_bytes())
            if digest!=hashlib.sha256(raw).hexdigest():
                raise ValueError("checksum mismatch")
            return tables_for(raw,ident,ids,req)
        except (ValueError,KeyError,StopIteration,ET.ParseError,OSError) as exc:
            errors.append(str(exc))
    # Read the verified public route again on the runner if receipts are missing.
    # A newer response gets its own hash and must pass the same schema checks.
    for endpoint in (host,"https://vizier.cds.unistra.fr/viz-bin/asu-xml","https://vizier.nao.ac.jp/viz-bin/asu-xml"):
        try:
            raw,req=arc.fetch(endpoint,params)
            return tables_for(raw,ident,ids,req)
        except (ValueError,KeyError,StopIteration,ET.ParseError,OSError) as exc:
            errors.append(str(exc))
    raise ValueError("Cannot verify "+ident+" catalogue: "+" | ".join(errors))


def names(row, fields):
    return {v for k in fields if len(v:=norm(row.get(k,"")))>3 and v not in {"unknown","none","unidentified"}}


def coord(row):
    try:
        ra,dec=float(row["_RAJ2000"]),float(row["_DEJ2000"])
        return (ra,dec) if 0<=ra<360 and -90<=dec<=90 else None
    except (ValueError,KeyError):
        return None


def posvalue(v):
    try:
        n=float(v)
        return n if math.isfinite(n) and n>0 else None
    except (ValueError,TypeError):
        return None


def matches(left,right):
    aliases=defaultdict(set)
    for j,b in enumerate(right):
        for name in names(b,("Name","AName")):
            aliases[name].add(j)
    options,reverse=defaultdict(set),defaultdict(set)
    distance={}
    for i,a in enumerate(left):
        xy=coord(a)
        if xy is None:
            continue
        near=set().union(*(aliases[n] for n in names(a,("MainID","PopularID"))))
        for j in near:
            pos=coord(right[j])
            if pos:
                arcsec=60*sky_distance(*xy,*pos)
                if arcsec<=RADIUS:
                    options[i].add(j)
                    reverse[j].add(i)
                    distance[i,j]=arcsec
    result=[]
    for i,js in options.items():
        if len(js)!=1:
            continue
        j=next(iter(js))
        if len(reverse[j])!=1:
            continue
        a,b=left[i],right[j]
        p,q=posvalue(a.get("Period")),posvalue(b.get("Porb"))
        ratio=p/q if p and q else None
        bin_name=("unavailable" if ratio is None else "near 24" if abs(ratio/24-1)<=WIDTH
                  else "near 1" if abs(ratio-1)<=WIDTH else "other")
        result.append(dict(fortin_main_id=a["MainID"],fortin_popular_id=a["PopularID"],avakyan_name=b["Name"],
            avakyan_alternate=b["AName"],fortin_ra_deg=a["_RAJ2000"],fortin_dec_deg=a["_DEJ2000"],
            avakyan_ra_deg=b["_RAJ2000"],avakyan_dec_deg=b["_DEJ2000"],separation_arcsec=distance[i,j],
            fortin_value=a.get("Period"),fortin_error=a.get("e_Period"),fortin_flag=a.get("f_Period"),
            fortin_bibcode=a.get("r_Period"),avakyan_value=b.get("Porb"),avakyan_bibcode=b.get("OrbRef"),
            ratio_fortin_over_avakyan=ratio,ratio_bin=bin_name))
    return sorted(result,key=lambda r:r["avakyan_name"].casefold())


def cds_doc(archive,ident):
    field=FIELDS[ident]
    urls=["https://cdsarc.cds.unistra.fr/ftp/"+ident+"/ReadMe",
          "https://cdsarc.u-strasbg.fr/ftp/"+ident+"/ReadMe",
          "https://cdsarc.cds.unistra.fr/ftp/cats/"+ident+"/ReadMe",
          "https://vizier.cfa.harvard.edu/viz-bin/ReadMe/"+ident+"?format=html&tex=false",
          "https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/"+ident+"?format=html&tex=false"]
    errors=[]
    for url in urls:
        try:
            raw,meta=archive.fetch(url)
            txt=raw.decode(errors="replace")
            if "<html" in txt[:2000].lower():
                txt=BeautifulSoup(txt,"html.parser").get_text("\n")
            for line in txt.splitlines():
                m=re.search(r"^\s*\d+\s*(?:-\s*\d+)?\s+\S+\s+(\S+)\s+"+field+r"\s+(.+)",line)
                if m and "orbit" in m[2].lower():
                    return dict(status="verified",unit=m[1],column_definition=line.strip(),response=meta)
            errors.append(url+" lacks a readable orbital-period line")
        except ValueError as exc:
            errors.append(str(exc))
    return dict(status="unresolved",attempts=errors)


def statements(lines):
    clock=r"(?:orbital periods?|\bP\s*[_-]?\s*orb\b)"
    unit=re.compile(clock+r"\s*(?:of the binary\s*)?(?:expressed\s+|provided\s+|given\s+)?(?:in\s+|\(\s*|\[\s*)(hours?\b|hrs?\b|h\b|days?\b|d\b)",re.I)
    found=[]
    for i,line in enumerate(lines):
        if not re.search(clock,line,re.I):
            continue
        txt=" ".join(lines[max(0,i-1):i+2])
        m=unit.search(txt)
        if m:
            found.append(" ".join(txt[max(0,m.start()-35):m.end()+45].split()))
    return list(dict.fromkeys(found))[:8]


def article(archive):
    rec=[]
    queries=[]
    def html(url,is_supplement=False):
        raw,meta=archive.fetch(url)
        soup=BeautifulSoup(raw,"html.parser")
        title=soup.title.get_text(" ",strip=True) if soup.title else ""
        paras=[n.get_text(" ",strip=True) for n in soup.select("p,tr,caption")]
        recognized="catalog" in title.lower() and bool(re.search("low.mass",title,re.I))
        if not recognized and not (is_supplement and soup.select_one("table")):
            return dict(status="unrecognized page",title=title,url=url,response=meta)
        links=[]
        for a in soup.select("a[href]"):
            href=a.get("href","")
            label=a.get_text(" ",strip=True)
            if not href.startswith("#") and re.search("T[0-9]+[.]html|/full_html/|github[.]com|zenodo[.]org",href,re.I):
                links.append([label,href])
        table_rows=[]
        for tr in soup.select("tr"):
            cells=[" ".join(td.stripped_strings) for td in tr.find_all(["td","th"],recursive=False)]
            joined=" ".join(cells)
            if re.search(r"orbital period|\bp\s*orb\b",joined,re.I):
                standalone=[c for c in cells if re.fullmatch(r"[(\[]?\s*(hours?|hrs?|h|days?|d)\s*[)\]]?",c,re.I)]
                if standalone:
                    table_rows.append(dict(cells=[x[:150] for x in cells[:12]],unit_cells=standalone))
        return dict(status="read",url=url,title=title,response=meta,
                    unit_statements=statements(paras)+table_rows,table_or_data_links=links[:30])
    base="https://www.aanda.org/articles/aa/full_html/2024/04/aa47908-23/aa47908-23.html"
    for u in [base,"https://doi.org/"+DOIS[FORTIN]]:
        try:
            rec.append(html(u))
        except ValueError as exc:
            rec.append(dict(url=u,error=str(exc)))
    visited={r["url"] for r in rec}
    for r in rec:
        for _,link in r.get("table_or_data_links",[]):
            url=urljoin(r["url"],link).split("#")[0]
            if url in visited or not url.startswith("https://www.aanda.org/"):
                continue
            # Follow the publisher's own complete article and tables only. The
            # article suffix guards against unrelated full-text citations.
            if not re.search(r"aa47908-23(?:[-_.]?T[0-9]+)?[.]html$",url,re.I):
                continue
            visited.add(url)
            try:
                rec.append(html(url, bool(re.search("T[0-9]+[.]html",url,re.I))))
            except ValueError as exc:
                rec.append(dict(url=url,error=str(exc)))
    try:
        params={"search_query":'au:Fortin AND ti:"catalogue of low-mass X-ray binaries"',"start":0,"max_results":10}
        raw,meta=archive.fetch("https://export.arxiv.org/api/query",params)
        queries.append(meta)
        for ent in ET.fromstring(raw).findall(ATOM+"entry"):
            title=" ".join(ent.findtext(ATOM+"title","").split())
            auth=[n.findtext(ATOM+"name","") for n in ent.findall(ATOM+"author")]
            ident=ent.findtext(ATOM+"id","").split("/abs/")[-1]
            if not re.search("low.mass",title,re.I) or not any("fortin" in a.lower() for a in auth):
                continue
            if not re.fullmatch(r"\d{4}\.\d{4,5}(?:v\d+)?",ident):
                continue
            pages,pdf=paper_text(archive,ident)
            rec.append(dict(status="read" if pages else "unresolved",title=title,authors=auth,
                            url="https://arxiv.org/abs/"+ident,response=pdf,
                            unit_statements=[dict(page=i,passages=hits) for i,p in enumerate(pages,1)
                                             if (hits:=statements(p.splitlines()))]))
    except (ValueError,ET.ParseError) as exc:
        queries.append(dict(error=str(exc)))
    updates=[]
    try:
        raw,receipt=archive.fetch("https://api.crossref.org/works/"+DOIS[FORTIN])
        w=json.loads(raw)["message"]
        updates.append(dict(response=receipt,doi=w.get("DOI"),relations=w.get("relation",{}),
                            update_to=w.get("update-to",[]),updated_by=w.get("updated-by",[])))
        raw,receipt=archive.fetch("https://api.crossref.org/works",{"query.title":
            "Corrigendum A catalogue of low-mass X-ray binaries in the Galaxy","query.author":"Fortin","rows":20})
        works=json.loads(raw)["message"]["items"]
        related=[dict(title=w.get("title"),doi=w.get("DOI")) for w in works
                 if re.search("corrigendum|erratum|correction"," ".join(w.get("title",[])),re.I)]
        updates.append(dict(response=receipt,ranked_limit=20,possible_corrections=related))
    except (ValueError,KeyError) as exc:
        updates.append(dict(error=str(exc)))
    return dict(records=rec,queries=queries,publisher_and_metadata_update_queries=updates)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    dest=args.output
    dest.mkdir(parents=True,exist_ok=True)
    lit=args.evidence/"literature"
    arc=Archive(dest/"documents")
    try:
        a,ma=catalogue(lit,FORTIN,["Period","MainID","PopularID","_RAJ2000","_DEJ2000"],arc)
        b,mb=catalogue(lit,AVAKYAN,["Porb","Name","AName","_RAJ2000","_DEJ2000"],arc)
        rows=matches(a,b)
        if not rows:
            raise ValueError("No reciprocal name-and-position matches")
    except (ValueError,TypeError,KeyError,IndexError,StopIteration,OSError,ET.ParseError) as exc:
        pretty(dest/"summary.json",{"status":"incomplete","reason":str(exc),"expected_first_run_hashes":FROZEN_DIGESTS})
        (dest/"README.md").write_text("# LMXB catalogue-unit audit\n\nThe source-table cross-check was inconclusive. See [retrieval diagnostics](summary.json).\n")
        print("Ancillary catalogue-unit audit incomplete:",exc)
        return
    with (dest/"source-comparison.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    docs={}
    for cat in (FORTIN,AVAKYAN):
        try:
            docs[cat]=cds_doc(arc,cat)
        except (ValueError,TypeError,KeyError,IndexError,OSError) as exc:
            docs[cat]=dict(status="unresolved",error=str(exc))
    try:
        primary=article(arc)
    except (ValueError,TypeError,KeyError,IndexError,OSError,ET.ParseError) as exc:
        primary=dict(records=[],queries=[],error=str(exc))
    pretty(dest/"primary-publication.json",primary)
    good=[r for r in rows if r["ratio_fortin_over_avakyan"] is not None]
    cnt=Counter(r["ratio_bin"] for r in rows)
    med=float(np.median([r["ratio_fortin_over_avakyan"] for r in good])) if good else None
    same_reference=[r for r in good if r.get("fortin_bibcode") and
                    r["fortin_bibcode"].casefold()==(r.get("avakyan_bibcode") or "").casefold()]
    sensitivity=[dict(relative_half_width=z,near_24=sum(abs(r["ratio_fortin_over_avakyan"]/24-1)<=z for r in good),
                 near_one=sum(abs(r["ratio_fortin_over_avakyan"]-1)<=z for r in good)) for z in (.01,.05,.10)]
    summary=dict(status="measured",fortin=ma,avakyan=mb,official_readmes=docs,reciprocal_matches=len(rows),both_positive_numeric=len(good),
                 counts=cnt,median_ratio_all=med,max_separation_arcsec=RADIUS,relative_bin_half_width=WIDTH,
                 common_primary_bibcode=len(same_reference),
                 common_primary_bibcode_counts=Counter(r["ratio_bin"] for r in same_reference),
                 tolerance_sensitivity=sensitivity)
    pretty(dest/"summary.json",summary)
    want={"Her X-1","Sco X-1","Cen X-4","4U 1820-30","Aql X-1","Cyg X-2"}
    selected=[r for r in good if want.intersection({r['fortin_popular_id'],r['avakyan_name'],r['avakyan_alternate']})]
    txt=["# Orbital-period units in two LMXB catalogues","",
         f"The Fortin Period and Avakyan Porb VOTable fields declare {ma['field_attributes']['unit']} and {mb['field_attributes']['unit']}. Among {len(good)} common objects with positive numerical periods, {cnt['near 24']} Fortin-to-Avakyan ratios fell within 5% of 24, {cnt['near 1']} within 5% of one, and {cnt['other']} outside both intervals.","",
         f"The match requires an exact normalized primary or alternative name, positions within {RADIUS:g} arcseconds, and a unique association in both directions. Coordinate signs are preserved. There were {len(rows)} reciprocal name-and-position matches before requiring positive periods. The [full comparison](source-comparison.csv) includes the flags and original period-reference bibcodes.","",
         "| Source | Fortin Period | Avakyan Porb | Ratio |","| --- | ---: | ---: | ---: |"]
    for r in selected:
        txt.append(f"| {r['avakyan_name']} | {r['fortin_value']} | {r['avakyan_value']} | {r['ratio_fortin_over_avakyan']:.4g} |")
    txt+=["","A numerical ratio near 24 is compatible with an hours-to-days conversion. The counts apply to these published table versions and matched records. Different ephemerides, ambiguous periods and aliases can affect individual ratios. The survey relies on independent published ephemerides for source-specific decisions.",""]
    txt.extend([f"The catalogues cite an identical nonempty period-reference bibcode for {len(same_reference)} numerical pairs. Of those pairs, {sum(r['ratio_bin']=='near 24' for r in same_reference)} have a ratio within 5% of 24. Separate catalogue entries can draw on the same original measurement.",""])
    for ident,m in ((FORTIN,ma),(AVAKYAN,mb)):
        d=docs[ident]
        info="The CDS ReadMe gives "+repr(d["unit"])+"." if d["status"]=="verified" else "The CDS ReadMe units line remains unresolved."
        txt.append(f"- {ident}, field {m['field_attributes']['name']}, UCD {m['field_attributes']['ucd']}, description “{m['description']}”, VOTable unit {m['field_attributes']['unit']}. {info}")
    statements_found=[p for p in primary["records"] if p.get("unit_statements")]
    txt+=["",f"The paper and supplement search found an explicit orbital-period unit statement in {len(statements_found)} retrieved records. [Publication receipts](primary-publication.json) and [catalogue metadata](summary.json).",""]
    (dest/"README.md").write_text("\n".join(txt))
    print(json.dumps(dict(overlap=len(rows),numeric=len(good),bins=cnt,cds={k:v["status"] for k,v in docs.items()},article_unit_records=len(statements_found))))


if __name__=="__main__":
    main()

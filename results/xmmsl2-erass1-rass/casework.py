"""Source dossiers from archival notices, object aliases and bright-star astrometry."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import html
import io
import json
from pathlib import Path
import re
import tarfile
import threading
from urllib.parse import quote, urljoin
from xml.etree import ElementTree

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from dateutil import parser as dateparser

PRIORITY = ("090822.4-643751","063045.7-603110","051627.4-650426","105848.7+082026","162701.0-694411","052443.6-613005","125851.2+235522","080026.9+101303")
HOST = re.compile(r"galaxy|AGN|Sy[12]|QSO|Seyfert|XRB|ULX|BL.?Lac|LINER|^G$",re.I)
VARIABLE = re.compile(r"\b(?:nova|supernova|outburst|tidal disruption|transient|flares?|flaring|quiescen\w*|turn(?:ed)?[- ](?:off|on)|changing[- ]look|variab\w*|optical loading|dimm\w*|brighten\w*)\b",re.I)
CANON = str.maketrans({c:"-" for c in "−–—‐﹣"})
ASTRO = re.compile(r"(?<!\d)J?\s*(\d{2})([0-5]\d)([0-5]\d(?:\.\d+|\d{1,2})?)([+\-])([0-8]\d)([0-5]\d)([0-5]\d(?:\.\d+|\d)?)")
GCN = "https://gcn.nasa.gov/circulars/"
ATEL = "https://www.astronomerstelegram.org/"
HEA = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
LOCK = threading.Lock()
RECEIPTS = []
AGENT = {"User-Agent":"Open-Space astronomical archival case review/1.0"}


def get(url, params=None, timeout=23, cap=30_000_000):
    rec = {"url":url,"params":params}
    try:
        with requests.get(url,params=params,headers=AGENT,timeout=(9,timeout),stream=True) as r:
            rec.update(status=r.status_code,final=r.url)
            r.raise_for_status()
            with io.BytesIO() as fp:
                for chunk in r.iter_content(131072):
                    fp.write(chunk)
                    if fp.tell()>cap:
                        raise ValueError("response exceeds size limit")
                b = fp.getvalue()
        rec.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
        return b
    except Exception as exc:
        rec["error"] = str(exc)[:260]
        raise
    finally:
        with LOCK:
            RECEIPTS.append(rec)


def clean(x):
    return " ".join(str(x or "").replace("\x00","").translate(CANON).split()).strip()


def norm(x):
    return re.sub("[^a-z0-9+-]","",clean(x).lower())


def coord_token(s):
    m=ASTRO.search(clean(s))
    if not m:
        return None
    h,mi,ss,sgn,d,dm,ds=m.groups()
    def second(x):
        return float(x if "." in x or len(x)==2 else x[:2]+"."+x[2:])
    if int(h)>23:
        return None
    return (int(h)+int(mi)/60+second(ss)/3600)*15, (1 if sgn=="+" else -1)*(int(d)+int(dm)/60+second(ds)/3600)


def degdist(a,b):
    x1,y1,x2,y2=np.deg2rad([a[0],a[1],b[0],b[1]])
    v=np.sin((y2-y1)/2)**2+np.cos(y1)*np.cos(y2)*np.sin((x2-x1)/2)**2
    return float(2*np.rad2deg(np.arcsin(np.sqrt(np.clip(v,0,1))))*3600)


def raw_text(s):
    return clean(BeautifulSoup(str(s),"html.parser").get_text(" "))


def search_aliases(row, ev):
    aliases=set()
    for k in ("name","ero_name","rass_name","xmm_counterpart","xmm_oldrass","xmm_altident","host_id"):
        for name in re.split(r"\s*\|\s*|;\s*",clean(row.get(k))):
            if 4<len(name)<110 and re.search("[A-Za-z]",name):
                aliases.add(name)
    for k,keys in (("NED",("preferred_name",)),("SIMBAD",("main_id",)),("WISE",("id",))):
        for item in ev.get("alternatives_within_35arcsec",{}).get(k,[]):
            if float(item.get("sep_arcsec",999))<=float(ev.get("gate_arcsec",15)):
                for key in keys:
                    a=clean(item.get(key))
                    if a:
                        if k=="WISE" and a[0].isdigit():
                            a="WISE J"+a
                        aliases.add(a)
    if "090822.4-643751" in row["name"]:
        # Supplied optical/IR identification is a query alias, not an inferred match.
        aliases.add("WISE J090822.49-643751.7")
    return sorted(aliases)


def bootstrap(out):
    leads=pd.read_csv(out/"input.csv",low_memory=False).fillna("")
    vt=pd.read_csv(out/"vetting/vetting.csv",low_memory=False).fillna("")
    joined=leads.merge(vt.drop(columns=["ra","dec"],errors="ignore"),on="id",validate="one_to_one")
    with gzip.open(out/"vetting/counterpart_evidence.jsonl.gz","rt") as f:
        evidence={r["id"]:r for r in map(json.loads,f)}
    cases=[]
    for r in joined.to_dict("records"):
        special=any(k in r["name"] for k in PRIORITY)
        promising=HOST.search(str(r.get("host_type",""))) or HOST.search(str(r.get("xmm_xtyp","")))
        unknown=clean(r.get("xmm_xtyp","")).lower() in ("x","unknown","unidentified","","ir","radio","candidate_agn")
        if not (special or promising or (r["vet_status"]!="known_or_stellar" and unknown)):
            continue
        ev=evidence.get(r["id"],{})
        case=dict(id=r["id"],name=r["name"],ra=r["ra"],dec=r["dec"],xmm_ra=r.get("xmm_ra",r["ra"]),
                  xmm_dec=r.get("xmm_dec",r["dec"]),xmm_position_error=r.get("xmm_poserr"),xmm_epoch=r["xmm_time"],
                  gate=ev.get("gate_arcsec",15),guarded_ratio=r["guarded_ratio"],high=r["high"],low=r["low"],
                  host=r.get("host_id"),host_type=r.get("host_type"),slew_class=r.get("xmm_xtyp"),
                  optical_loading_catalog_flag=r.get("xmm_veropt"),alias=search_aliases(r,ev),
                  priority=special,source_status=r["vet_status"],evidence=[],checks={})
        for cit in ev.get("bibliography_hits",[]):
            if cit.get("bibcode"):
                case["evidence"].append(dict(kind="bibliography of positional SIMBAD object",
                    title=cit.get("title"),url="https://ui.adsabs.harvard.edu/abs/"+quote(cit["bibcode"]),
                    bibcode=cit["bibcode"],interpretation="object-associated paper; epoch and association require inspection"))
        cases.append(case)
    return cases


def tap(url,sql):
    raw=get(url,{"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"votable","QUERY":sql,"MAXREC":130000},timeout=65)
    vt=parse(io.BytesIO(raw))
    for res in vt.resources:
        for info in res.infos:
            if info.name=="QUERY_STATUS" and str(info.value).upper() in ("OVERFLOW","ERROR"):
                raise ValueError(clean(info.content or info.value))
    return vt.get_first_table().to_table(use_names_over_ids=True).to_pandas()


def bright_stars(cases):
    # The Hipparcos Main Catalogue provides epoch-1991.25 positions and V.
    # Propagate angular offsets to each slew epoch. A bright nearby star is a
    # competing hypothesis; this proximity alone cannot diagnose optical loading.
    table=None
    failure=[]
    for url in [HEA]:
        try:
            names=tap(url,"SELECT column_name,description,unit FROM TAP_SCHEMA.columns WHERE table_name='hipparcos'")
            lookup={str(x).lower():str(x) for x in names.column_name}
            def field(*options):
                return next(lookup[x] for x in options if x in lookup)
            fields=dict(ra=field("ra"),dec=field("dec"),name=field("hip_number","hip","name"),
                        vmag=field("vmag","v_mag","v_magnitude"),pmra=field("pm_ra","pmra","proper_motion_ra"),
                        pmdec=field("pm_dec","pmdec","proper_motion_dec"))
            sql=f"SELECT {','.join(fields.values())} FROM hipparcos WHERE {fields['vmag']} <= 11"
            table=tap(url,sql).rename(columns={v:k for k,v in fields.items()})
            if len(table)<10000:
                raise ValueError("Hipparcos query returned unexpectedly few rows")
            # HEASARC public coordinates are commonly standardized to the J2000
            # epoch. Use the original Hipparcos data whenever its epoch is explicit.
            epoch_note = str(names.loc[names.column_name==fields["ra"],"description"].iloc[0])
            epoch = 1991.25 if "1991.25" in epoch_note else None
            if epoch is None:
                raise ValueError("HEASARC position epoch unspecified; using original release")
            break
        except Exception as e:
            failure.append(str(e)[:140])
            table=None
    if table is None:
        try:
            raw=get("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",{"-source":"I/239/hip_main","-out":"HIP,Vmag,RAICRS,DEICRS,pmRA,pmDE","-out.max":"unlimited"},timeout=80)
            table=pd.read_csv(io.BytesIO(raw),sep="\t",comment="#",dtype=str,skip_blank_lines=True)
            table.columns=[x.strip() for x in table.columns]
            table=table.rename(columns={"HIP":"name","Vmag":"vmag","RAICRS":"ra","DEICRS":"dec","pmRA":"pmra","pmDE":"pmdec"})
            if len(table)<118218 or not {"ra","dec","pmra","pmdec","name","vmag"}<=set(table):
                raise ValueError("primary Hipparcos survey schema/size failed")
            epoch=1991.25
        except Exception as exc:
            table=None
            failure.append(str(exc)[:160])
    if table is None:
        urls=["https://cdsarc.cds.unistra.fr/ftp/I/239/hip_main.dat.gz",
              "https://vizier.cfa.harvard.edu/ftp/cats/I/239/hip_main.dat.gz",
              "https://cdsarc.u-strasbg.fr/ftp/I/239/hip_main.dat.gz",
              "https://heasarc.gsfc.nasa.gov/FTP/heasarc/dbase/tdat_files/heasarc_hipparcos.tdat.gz"]
        # The final HEASARC TDAT file is not in the original CDS column format.
        for url in urls[:-1]:
            try:
                raw=gzip.decompress(get(url,timeout=75,cap=90_000_000))
                table=pd.read_csv(io.BytesIO(raw),sep="|",header=None,usecols=[1,5,8,9,12,13],
                                  names=None,low_memory=False)
                table.columns=["name","vmag","ra","dec","pmra","pmdec"]
                if len(table)<117000:
                    raise ValueError("Hipparcos primary file incomplete")
                epoch=1991.25
                break
            except Exception as e:
                table=None
                failure.append(str(e)[:140])
    if table is None:
        for c in cases:
            c["checks"]["Hipparcos"]="unavailable: "+"; ".join(failure)
        return
    for key in ("vmag","ra","dec","pmra","pmdec"):
        table[key]=pd.to_numeric(table[key],errors="coerce")
    table=table.loc[table.ra.between(0,360)&table.dec.between(-90,90)&table.vmag.le(11)].copy()
    check=table.loc[pd.to_numeric(table.name,errors="coerce").eq(91262)]
    if len(table)<100000 or len(check)!=1 or degdist((float(check.iloc[0].ra),float(check.iloc[0].dec)),(279.2347,38.7837))>60:
        for c in cases:
            c["checks"]["Hipparcos"]="incomplete catalogue or failed Vega position control"
        return
    position=SkyCoord(table.ra.to_numpy(),table.dec.to_numpy(),unit="deg")
    for c in cases:
        centre=SkyCoord(c["xmm_ra"],c["xmm_dec"],unit="deg")
        pool=table.loc[position.separation(centre).arcsec<600]
        t=pd.to_datetime(c["xmm_epoch"],errors="coerce")
        if pd.isna(t):
            c["checks"]["Hipparcos"]="missing slew epoch"
            continue
        years=t.year+(t.dayofyear-1)/365.25-epoch
        for r in pool.itertuples(index=False):
            start=SkyCoord(r.ra,r.dec,unit="deg")
            motion_known=np.isfinite(r.pmra) and np.isfinite(r.pmdec)
            p=start.spherical_offsets_by(r.pmra*years*u.mas,r.pmdec*years*u.mas) if motion_known else start
            separation=p.separation(centre).arcsec
            if separation<=65:
                c["evidence"].append(dict(kind="bright Hipparcos star near slew position",
                    title="HIP "+str(r.name),v_mag=float(r.vmag),separation_arcsec=round(float(separation),3),
                    pmra_mas_yr=float(r.pmra) if motion_known else None,pmdec_mas_yr=float(r.pmdec) if motion_known else None,
                    epoch=epoch,url="https://cdsarc.cds.unistra.fr/viz-bin/cat/I/239",
                    interpretation="inspect stellar flare/optical-loading alternative; brightness and separation alone do not establish detector loading"))
        c["checks"]["Hipparcos"]="V≤11; 65 arcsec at slew epoch, including PSF wings; 600 arcsec motion preselection"


def body_hit(c,text):
    """Match a precise name or angular coincidence with a coordinate-coded name."""
    text=clean(text)
    compact=norm(text)
    for alias in c["alias"]:
        if len(norm(alias))>=6 and re.search(re.escape(norm(alias)) + (r"(?![0-9])" if alias[-1].isdigit() else ""),compact):
            return alias,"exact supplied or catalogue alias"
    for m in ASTRO.finditer(text):
        p=coord_token(m.group())
        if p is not None and degdist(p,(c["ra"],c["dec"]))<=min(35,c["gate"]):
            return m.group(),"position-coded designation in search circle"
    return None


def add_publication(c,source,title,date,url,body):
    title,body=clean(title),clean(body)
    hit=body_hit(c,title+" "+body)
    if not hit:
        return
    token,how=hit
    located=clean(body).lower().find(clean(token).lower())
    if located<0:
        alt=ASTRO.search(token)
        if alt:
            located=body.find(alt.group()[0:6])
    around=body[max(0,located-115):located+215] if located>=0 else ""
    flags=sorted(set(m.group().lower() for m in VARIABLE.finditer(title+" "+around)))
    statement=("source-named event publication; verify which epoch and instrument" if body_hit(c,title) and flags
               else "event terms in local source context; review full publication" if flags
               else "source mentioned; reported epoch variation unestablished")
    c["evidence"].append(dict(kind=source,title=title[:650],date=date,url=url,matched=token,
                              match=how,event_terms=flags,context=around[:270],interpretation=statement))


def words_for_alias(c):
    words=[]
    for alias in c["alias"]:
        # Shared coarse designations are search keys only. A publication
        # without a precise alias or in-gate coordinates is never accepted.
        m=ASTRO.search(alias)
        if m:
            token=m.group(0).replace(" ","").lstrip("J")
            if re.search(r"XMMSL(?:1|2)?|1XMMSL",alias,re.I):
                words += [alias,"J"+token]
            words += ["J"+token[:6], "J"+m[1]+m[2]+m[4]+m[5]+m[6]]
        if (not m or not re.match(r"WISE|\dRXS|Gaia|AllWISE|1eRASS|XMMSL",alias,re.I)) and len(alias)<80:
            words.append(alias)
    return list(dict.fromkeys(words))


def query_primary(c):
    # arXiv records expose author-supplied titles and abstracts. Search ADS
    # links are supplied independently; lack of an arXiv hit is not clearance.
    unit="https://export.arxiv.org/api/query"
    terms=words_for_alias(c)
    if not terms:
        return
    sample=terms[:8]
    q=" OR ".join('all:"'+t.replace('"','')+'"' for t in sample)
    try:
        doc=ElementTree.fromstring(get(unit,{"search_query":q,"start":0,"max_results":40},timeout=18))
        if doc.tag.split("}")[-1]!="feed":
            raise ValueError("unexpected arXiv result format")
        entries=doc.findall("{*}entry")
        if any("error" in (e.findtext("{*}id") or "") for e in entries):
            raise ValueError("arXiv rejected the search")
        for e in entries:
            add_publication(c,"arXiv author abstract",e.findtext("{*}title",""),
                            e.findtext("{*}published",""),e.findtext("{*}id",""),e.findtext("{*}summary",""))
        c["checks"]["arXiv"]=dict(terms=sample,returned=len(entries),scope="source aliases; title and abstract")
    except Exception as e:
        c["checks"]["arXiv"]={"error":str(e)[:200],"terms":sample}
    # ADS indexes ATels/GCNs, which are sparsely represented by DOI registries.
    # ADS search requires authenticated access. Supply reproducible deep links
    # even when no token is available, without claiming a completed search.
    c["ADS_search_url"]="https://ui.adsabs.harvard.edu/search/q="+quote(" OR ".join('full:"'+s+'"' for s in sample),safe="")
    # Crossref searches here include coarse aliases for papers and publisher
    # records not selected by the exact-host search in history.py.
    if c["priority"] or HOST.search(c["host_type"] or c["slew_class"] or ""):
        counts=[]
        for term in sample[:3]:
            try:
                r=json.loads(get("https://api.crossref.org/works",{"query.bibliographic":term,"rows":25,
                                 "select":"DOI,title,abstract,published"},timeout=24))
                works=r["message"]["items"]
                for w in works:
                    title=" ".join(w.get("title",[]))
                    date="-".join(str(x) for x in (w.get("published",{}).get("date-parts") or [[]])[0])
                    add_publication(c,"publisher metadata via Crossref",title,date,
                                    "https://doi.org/"+w.get("DOI",""),raw_text(w.get("abstract","")))
                counts.append({"term":term,"returned":len(works)})
            except Exception as e:
                counts.append({"term":term,"error":str(e)[:200]})
        c["checks"]["Crossref case queries"]=counts


def gcn_notices(cases):
    missing=[]
    raw=None
    for url in (GCN+"all.tar.gz","https://gcn.gsfc.nasa.gov/gcn3/all_gcn_circulars.tar.gz",GCN+"all.jsonl"):
        try:
            b=get(url,timeout=80,cap=240_000_000)
            if not (b.startswith(b"\x1f\x8b") or b.lstrip().startswith(b"{")):
                raise ValueError("not a GCN notice archive")
            raw=b
            break
        except Exception as e:
            missing.append(str(e)[:145])
    if raw is None:
        for c in cases:
            c["checks"]["GCN"]="bulk archive unavailable: "+"; ".join(missing)
        return
    def entries():
        if raw.startswith(b"\x1f\x8b"):
            with tarfile.open(fileobj=io.BytesIO(raw),mode="r|gz") as tf:
                for member in tf:
                    if member.isfile() and 20<member.size<2_000_000:
                        yield member.name,tf.extractfile(member).read().decode("utf-8","replace")
        else:
            for line in raw.splitlines():
                try:
                    r=json.loads(line)
                    yield str(r.get("circularId",r.get("id",""))),json.dumps(r,ensure_ascii=False)
                except ValueError:
                    continue
    n=0
    ids=[]
    # Quickly select a small subset of notices by precise RA seconds, short
    # coordinate designation or a distinct galaxy name before full checks.
    shortlist=defaultdict(list)
    for c in cases:
        for t in words_for_alias(c):
            if len(t)>=5:
                shortlist[norm(t).lstrip("j")].append(c)
        for a in c["alias"]:
            m=ASTRO.search(a)
            if m:
                for offset in range(-2,3):
                    s=int(float(m[3][:2]))+offset
                    if 0<=s<60:
                        shortlist[m[1]+m[2]+f"{s:02d}"].append(c)
    test=re.compile("|".join(re.escape(k) for k in sorted(shortlist,key=len,reverse=True))) if shortlist else None
    try:
        for filename,body in entries():
            code=re.search(r"(\d{1,6})(?:\.gcn3|\.txt|\.md|\.json)?$",filename)
            if not code:
                code=re.search(r'(?:^NUMBER\s*:\s*|"circularId"\s*:\s*)(\d{1,6})',body,re.M|re.I)
            if not code:
                continue
            ident=int(code[1])
            n+=1
            ids.append(ident)
            flat=norm(body)
            if not test:
                continue
            use={c["id"]:c for m in test.finditer(flat) for c in shortlist[m.group()]}
            if not use:
                continue
            title=re.search(r"^\s*SUBJECT\s*:\s*(.+)$",body,re.I|re.M) or re.search(r"^\s*TITLE\s*:\s*(.+)$",body,re.I|re.M)
            date=re.search(r"^\s*DATE\s*:\s*(.+)$",body,re.I|re.M)
            if body.lstrip().startswith("{"):
                try:
                    j=json.loads(body)
                    title=j.get("subject",j.get("title","GCN "+str(ident)))
                    date=str(j.get("createdOn",j.get("date","")))
                    body=str(j.get("body",j.get("text",body)))
                except ValueError:
                    pass
            for c in use.values():
                add_publication(c,"NASA GCN Circular",title[1] if isinstance(title,re.Match) else title or "GCN "+str(ident),
                                date[1] if isinstance(date,re.Match) else date or "",GCN+str(ident),body)
        if n<1000:
            raise ValueError(f"only {n} identifiable circulars")
        check=dict(source=url,parsed=n,lowest_id=min(ids),highest_id=max(ids),scope="designations in notice text")
    except Exception as e:
        check=dict(source=url,parsed=n,error=str(e)[:240])
    for c in cases:
        c["checks"]["GCN"]=check


def simbad_aliases(cases):
    mirrors=("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
             "https://simbad.cds.unistra.fr/simbad/sim-tap/sync")
    active=None
    for base in mirrors:
        try:
            probe=tap(base,"SELECT TOP 1 basic.main_id,ids.ids FROM basic JOIN ids ON basic.oid=ids.oidref WHERE basic.main_id='3C 273'")
            if len(probe) and "3c273" in norm(str(probe.iloc[0].to_dict())):
                active=base
                break
        except Exception:
            continue
    if active is None:
        for c in cases:
            c["checks"]["SIMBAD historical aliases"]="identifier table unavailable or positive control failed"
        return
    try:
        meta=tap(active,"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='ref'")
        refcols={str(c).lower() for c in meta.column_name}
    except Exception:
        refcols=set()
    for start in range(0,len(cases),15):
        part=cases[start:start+15]
        circles=" OR ".join(f"1=CONTAINS(POINT('ICRS',basic.ra,basic.dec),CIRCLE('ICRS',{c['ra']:.7f},{c['dec']:.7f},{min(24,c['gate'])/3600:.8f}))" for c in part)
        try:
            rows=tap(active,"SELECT basic.oid,basic.main_id,basic.otype,basic.ra,basic.dec,ids.ids FROM basic JOIN ids ON basic.oid=ids.oidref WHERE "+circles)
            for c in part:
                choices=[]
                for r in rows.itertuples(index=False):
                    separation=degdist((r.ra,r.dec),(c["ra"],c["dec"]))
                    if separation<=min(24,c["gate"]):
                        names=[clean(n) for n in str(r.ids).split("|") if re.search(r"XMMSL|^1?RXS|^2RX|^Swift|^MAXI|^AT |^SN |^GALEX|^IRAS|^PGC |^NGC|^2MASX|^Mrk|^LEDA|^ESO |^HD |^HIP |^V\*|^RX |^IGR |^GRB|^ASASSN",n.strip(),re.I)]
                        choices.append((separation,clean(r.main_id),names,clean(r.otype),int(r.oid)))
                for d,m,aliases,typ,oid in choices:
                    c["evidence"].append(dict(kind="SIMBAD historical identifiers",title=m,otype=typ,
                        aliases=aliases,separation_arcsec=round(d,3),url="https://simbad.cds.unistra.fr/simbad/sim-id?Ident="+quote(m),
                        interpretation="positional object; use its aliases only after inspecting the association"))
                if len(choices)==1:
                    c["alias"]=sorted(set(c["alias"]+choices[0][2]+[choices[0][1]]))
                    c["simbad_oid"]=choices[0][4]
                c["checks"]["SIMBAD historical aliases"]=dict(objects_in_gate=len(choices),endpoint=active)
            oids=sorted({c["simbad_oid"] for c in part if c.get("simbad_oid")})
            if oids:
                optional="".join(",ref."+field for field in ("abstract","doi","journal") if field in refcols)
                sql=("SELECT has_ref.oidref,ref.bibcode,ref.title"+optional+
                     " FROM has_ref JOIN ref ON has_ref.oidbibref=ref.oidbib WHERE has_ref.oidref IN ("+
                     ",".join(map(str,oids))+")")
                try:
                    hits=tap(active,sql)
                    for c in part:
                        if not c.get("simbad_oid"):
                            continue
                        subset=hits[hits.oidref==c["simbad_oid"]]
                        c["checks"]["SIMBAD primary references"]={"records":len(subset),"metadata":sorted(refcols & {"abstract","doi","title"})}
                        for r in subset.fillna("").to_dict("records"):
                            title,body,code=clean(r["title"]),clean(r.get("abstract","")),clean(r["bibcode"])
                            if not (VARIABLE.search(title+" "+body) or "ATel" in code or "GCN" in code):
                                continue
                            link="https://doi.org/"+clean(r["doi"]) if clean(r.get("doi")) else "https://ui.adsabs.harvard.edu/abs/"+quote(code)
                            c["evidence"].append(dict(kind="SIMBAD object reference",title=title,bibcode=code,date=code[:4],
                                url=link,interpretation="published event subject associated with this object's bibliography; inspect scope and epoch"))
                            add_publication(c,"source-named SIMBAD paper",title,code[:4],link,body)
                except Exception as exc:
                    for c in part:
                        c["checks"]["SIMBAD primary references"]={"error":str(exc)[:170]}
        except Exception as exc:
            for c in part:
                c["checks"]["SIMBAD historical aliases"]={"error":str(exc)[:200]}


def atel_notices(cases):
    import functools
    from datetime import timedelta
    months=r"(?:Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sept?(?:ember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)"
    dates=re.compile(r"\b(\d{1,2}\s+"+months+r"\.?\s+\d{4})\b",re.I)

    @functools.lru_cache(maxsize=1600)
    def read(n):
        url=ATEL+"?read="+str(n)
        page=BeautifulSoup(get(url,timeout=17,cap=2_000_000),"html.parser")
        # Retain the individual telegram; remove menus, related-notice lists
        # and site-wide links, which can name an unrelated object.
        heading=clean(page.title.get_text(" ") if page.title else "")
        heading=re.sub(r"^ATel\s*#?\s*"+str(n)+r"\s*[:\-]?\s*","",heading,flags=re.I)
        main=page.select_one("#telegram, #telegramcontent, #telegramContent, .telegram, article")
        if main is None:
            main=page.select_one("main") or page
        for tag in main.select("header, footer, aside, nav, script, style, .related, #related"):
            tag.decompose()
        body=clean(main.get_text(" "))
        date=None
        for m in dates.finditer(body):
            context=body[max(0,m.start()-8):m.start()].lower()
            if "on " in context or date is None:
                date=dateparser.parse(m[1],dayfirst=True).date()
                if "on " in context:
                    break
        return heading,date,body,url

    try:
        front=BeautifulSoup(get(ATEL+"?browse",timeout=18,cap=8_000_000),"html.parser")
        links={}
        for a in front.find_all("a",href=True):
            m=re.search(r"(?:\?|&)read=(\d+)",a["href"])
            if m:
                links[int(m[1])]=clean(a.get_text(" "))
        if not links:
            raise ValueError("no numbered telegrams in index")
        for n,title in links.items():
            for c in cases:
                if body_hit(c,title):
                    try:
                        title,date,body,url=read(n)
                        add_publication(c,"The Astronomer's Telegram",title,str(date or ""),url,body)
                    except Exception:
                        continue
        hi=max(links)+1
        for c in cases:
            c["checks"]["ATel title index"]=dict(parsed_titles=len(links),url=ATEL+"?browse")
    except Exception as exc:
        for c in cases:
            c["checks"]["ATel title index"]={"error":str(exc)[:200]}
        return

    def near(day):
        low,high=1,hi
        while low<high:
            mid=(low+high)//2
            _,date,_,_=read(mid)
            if date is None:
                raise ValueError(f"no publication date on ATel {mid}")
            if date<day:
                low=mid+1
            else:
                high=mid
        return low
    for c in cases:
        # Also inspect ATels cited in the object's own SIMBAD bibliography,
        # irrespective of the year or date-window restriction.
        cited=set()
        for entry in c["evidence"]:
            m=re.search(r"ATel\.+(\d+)",entry.get("bibcode",""))
            if m:
                cited.add(int(m[1]))
        if c["priority"]:
            window={}
            try:
                epoch=dateparser.parse(str(c["xmm_epoch"])).date()
                start,end=epoch-timedelta(days=5),epoch+timedelta(days=31)
                lower,upper=near(start),near(end+timedelta(days=1))
                window.update(start=str(start),end=str(end),first_id=lower,last_id=upper-1)
                if not 0<=upper-lower<=400:
                    raise ValueError("unexpected number of notices in date window")
                cited.update(range(lower,upper))
            except Exception as exc:
                window["error"]=str(exc)[:200]
            c["checks"]["ATel slew-date window"]=window
        good,bad=0,0
        for n in sorted(cited):
            try:
                title,date,body,url=read(n)
                add_publication(c,"The Astronomer's Telegram",title,str(date or ""),url,body)
                good+=1
            except Exception:
                bad+=1
        c["checks"]["ATel individual telegrams"]=dict(inspected=good,unavailable=bad)


def main(out):
    dest=out/"casework"
    dest.mkdir(exist_ok=True)
    cases=bootstrap(out)
    # Expand aliases before searching notices. Independent providers have
    # separate failure records and are never used to clear another provider.
    simbad_aliases(cases)
    with ThreadPoolExecutor(max_workers=4) as pool:
        parallel={name:pool.submit(fn,cases) for name,fn in [
            ("Hipparcos",bright_stars),("GCN archive",gcn_notices),("ATel archive",atel_notices)]}
        priors=sorted(cases,key=lambda c:(not c["priority"],-c["guarded_ratio"]))
        for c in priors:
            # Keep the primary case sample and all galaxies; the bulk GCN,
            # ATel title, host and positional checks still cover every case.
            if c["priority"] or HOST.search(c["host_type"] or "") or HOST.search(c["slew_class"] or ""):
                try:
                    query_primary(c)
                except Exception as e:
                    c["checks"]["primary search"]={"error":str(e)[:220]}
        for name,job in parallel.items():
            try:
                job.result()
            except Exception as e:
                for c in cases:
                    c["checks"][name]={"error":str(e)[:240]}
    rows=[]
    with (dest/"casebook.jsonl").open("w") as fp:
        for c in sorted(cases,key=lambda c:(not c["priority"],-c["guarded_ratio"])):
            unique={}
            for e in c["evidence"]:
                key=(e["kind"],e.get("url"),e.get("title"))
                unique[key]=e
            c["evidence"]=list(unique.values())
            events=[e for e in c["evidence"] if e.get("event_terms")]
            if events:
                status="possible published event; inspect source, instrument, and epoch"
            elif any(e["kind"].startswith("bright Hipparcos") for e in c["evidence"]):
                status="bright-star alternative; inspect astrometry and slew spectrum"
            else:
                status="no source-specific prior event established by these checks"
            c["interpretation"]=status
            fp.write(json.dumps(c,allow_nan=False)+"\n")
            rows.append(dict(detection_uid=c["id"],name=c["name"],priority=c["priority"],xmm_epoch=c["xmm_epoch"],
                             guarded_ratio=c["guarded_ratio"],interpretation=status,prior_event_leads=len(events),
                             publications=" | ".join(sorted({e["url"] for e in events}))))
    pd.DataFrame(rows).to_csv(dest/"case_index.csv",index=False)
    with gzip.open(dest/"retrievals.jsonl.gz","wt") as fp:
        for rec in RECEIPTS:
            fp.write(json.dumps(rec)+"\n")
    print(f"Case dossiers: {len(cases)} sources, {len(RECEIPTS)} archival requests",flush=True)


def historical_slew(cases, leads):
    """Resolve historical XMMSL1 names at the original slew epoch."""
    from astropy.io import fits
    master=leads.set_index("id")
    for release in ("D6","D5","D4"):
        url=f"https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/xmmsl1{release}_total.fits.gz"
        try:
            raw=get(url,timeout=75,cap=90_000_000)
            with fits.open(io.BytesIO(gzip.decompress(raw)),memmap=False) as h:
                d=h[1].data
                cols={re.sub("[^a-z0-9]","",str(k).lower()):k for k in d.names}
                def col(*keys):
                    return next(cols[k] for k in keys if k in cols)
                rc,dc=col("ra","raj2000"),col("dec","dej2000")
                oc=col("obsid","observationid","observationnumber","observation")
                names=[cols[k] for k in ["srcname","sourcename","uniquesrcname","iauname","name"] if k in cols]
                if not names:
                    raise ValueError("missing historical source designation columns")
                ids=np.array([clean(x.decode() if isinstance(x,bytes) else x).lstrip("0") for x in d[oc]])
                coords=SkyCoord(np.array(d[rc],dtype=float)*u.deg,np.array(d[dc],dtype=float)*u.deg)
                for c in cases:
                    obs=str(master.loc[c["id"],"xmm_obs"]).split(".")[0].strip().lstrip("0")
                    indices=np.flatnonzero(ids==obs)
                    center=SkyCoord(c["xmm_ra"],c["xmm_dec"],unit="deg")
                    hits=[]
                    for k in indices:
                        sep=coords[k].separation(center).arcsec
                        if sep>max(12,min(35,c["gate"])):
                            continue
                        aliases=set()
                        for n in names:
                            rawname=d[n][k]
                            a=clean(rawname.decode() if isinstance(rawname,bytes) else rawname)
                            if "XMMSL1" in a.upper():
                                aliases.add(a)
                        if aliases:
                            hits.append((float(sep),sorted(aliases)))
                    if len(hits)==1:
                        sep,aliases=hits[0]
                        new=sorted(set(aliases)-set(c["alias"]))
                        c["alias"]=sorted(set(aliases+c["alias"]))
                        c["evidence"].append(dict(kind=f"XMMSL1{release} same-observation association",title=" | ".join(aliases),
                             url=url,observation_id=obs,separation_arcsec=round(sep,3),additional_query_aliases=new,
                             interpretation="published catalogue measurement at same slew epoch; does not establish a prior variability report"))
                    c["checks"][f"XMMSL1{release}"]=dict(same_observation_detections=len(indices),positional_matches=len(hits))
        except Exception as e:
            for c in cases:
                c["checks"][f"XMMSL1{release}"]={"error":str(e)[:220]}


def literature_documents(record):
    """Route a cited primary paper to open article text, keeping genuine identifiers."""
    import difflib
    import shutil
    import subprocess
    title=record["title"]
    canon=lambda t: re.sub(r"[^a-z0-9 ]"," ",clean(t).lower()).split()
    status=[]
    if record.get("url","").startswith("https://doi.org/"):
        doi=record["url"].split("doi.org/",1)[1]
    else:
        q=json.loads(get("https://api.crossref.org/works",{"query.title":title,"rows":3,"select":"DOI,title"}))
        doi=""
        for r in q["message"]["items"]:
            sim=difflib.SequenceMatcher(None,canon(title),canon(" ".join(r.get("title",[])))).ratio()
            if sim>=.86:
                doi=clean(r.get("DOI"))
                break
    if not doi:
        return [],["no verified DOI resolution"]
    urls=["https://doi.org/"+doi]
    try:
        meta=json.loads(get("https://api.openalex.org/works/https://doi.org/"+doi))
        for item in [meta.get("best_oa_location")]+meta.get("locations",[]):
            if item and item.get("is_oa"):
                for ukey in ("pdf_url","landing_page_url"):
                    addr=clean(item.get(ukey))
                    if addr.startswith("https://") and addr not in urls:
                        urls.append(addr)
        status.append("metadata via OpenAlex for publisher DOI")
    except Exception as e:
        status.append("OpenAlex: "+str(e)[:150])
    # Prefer the primary publisher and arXiv copies when present.
    urls=sorted(urls,key=lambda x:("doi.org" not in x and "arxiv.org" not in x,not x.endswith(".pdf")))
    docs=[]
    for url in list(dict.fromkeys(urls))[:5]:
        try:
            b=get(url,timeout=30,cap=45_000_000)
            if b.startswith(b"%PDF") and shutil.which("pdftotext"):
                text=subprocess.run(["pdftotext","-layout","-","-"],input=b,stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL,check=True,timeout=45).stdout.decode("utf8","replace")
            elif b[:3000].lstrip().startswith(b"<"):
                dom=BeautifulSoup(b,"html.parser")
                for node in dom.select("nav,footer,script,style,header,aside"):
                    node.decompose()
                main=dom.select_one("article,main") or dom
                text=main.get_text(" ",strip=True)
                for a in dom.find_all("a",href=True):
                    ref=urljoin(url,a["href"])
                    if ref.startswith("https://") and "/full_html/" in ref and ref not in urls:
                        try:
                            raw=get(ref,timeout=30,cap=12_000_000)
                            htmltext=BeautifulSoup(raw,"html.parser").get_text(" ",strip=True)
                            if len(htmltext)>len(text):
                                docs.append((ref,htmltext))
                        except Exception:
                            pass
                        break
            else:
                status.append("unsupported content: "+url)
                continue
            if len(text)>400 and not re.search(r"^.{0,120}(?:captcha|access denied|checking your browser)",text,re.I):
                docs.append((url,text))
            if len(text)>18000:
                break
        except Exception as e:
            status.append(str(e)[:150])
    return docs,status


def historical(out):
    dest=out/"casework"
    with (dest/"casebook.jsonl").open() as fp:
        cases=[json.loads(r) for r in fp]
    leads=pd.read_csv(out/"input.csv",low_memory=False)
    for c in cases:
        c["priority"]=any(t in c["name"] for t in PRIORITY)
    historical_slew(cases,leads)
    with gzip.open(out/"vetting/counterpart_evidence.jsonl.gz","rt") as fp:
        snapshots={c["id"]:c for c in map(json.loads,fp)}
    for c in cases:
        ev=snapshots.get(c["id"],{})
        c["catalog_counterpart_observables"]={k:v for k,v in ev.get("alternatives_within_35arcsec",{}).items()
                                             if k in ("Gaia","WISE","Legacy","NED","SIMBAD")}
    recovered=[]
    for c in cases:
        aliases=sorted({a for e in c["evidence"] for a in e.get("additional_query_aliases",[])})
        if aliases and (c["priority"] or HOST.search(c["host_type"] or "") or HOST.search(c["slew_class"] or "")):
            c["alias"]=aliases+[a for a in c["alias"] if a not in aliases]
            recovered.append(c)
    priority={c["id"]:c for c in cases if c["priority"]}
    for c in recovered:
        priority[c["id"]]=c
    if priority:
        focused=list(priority.values())
        named_object_references([c for c in focused if c["priority"]])
        gcn_notices(focused)
        atel_publication_index(focused)
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(query_primary,focused))
    targets=defaultdict(list)
    for c in cases:
        if not (c["priority"] or HOST.search(c["host_type"] or "") or HOST.search(c["slew_class"] or "")):
            continue
        refs=[r for r in c["evidence"] if r.get("kind")=="SIMBAD object reference" and
              r.get("title") and not re.search(r"ATel|GCN",r.get("bibcode",""))]
        refs.sort(key=lambda r:(not bool(re.search("slew|eROSITA|RASS|ROSAT",r["title"],re.I)),
                                not bool(re.search("tidal|outburst|nova|X-ray transient",r["title"],re.I))))
        # The all-sky variability studies may list this source in the primary
        # body or table even when an object bibliography omitted the link.
        if c["priority"]:
            refs=[{"title":"Populations of highly variable X-ray sources in the XMM–Newton slew survey","url":"https://doi.org/10.1093/mnras/stac714","date":"2022"},
                  {"title":"Highly variable AGN from the XMM-Newton slew survey","url":"https://doi.org/10.1051/0004-6361/201628241","date":"2016"}]+refs
        for r in refs[:10]:
            targets[r["title"]].append((c,r))
    def fetch_one(title):
        try:
            return title,*literature_documents(targets[title][0][1])
        except Exception as e:
            return title,[],[str(e)[:180]]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for title,docs,status in pool.map(fetch_one,targets):
            for c,r in targets[title]:
                c["checks"].setdefault("article body searches",[]).append(dict(title=title,
                    readable_locations=[url for url,_ in docs],status=status))
                for url,body in docs:
                    add_publication(c,"cited primary paper body",title,r.get("date",""),url,body)
    with (dest/"historical_casebook.jsonl").open("w") as fp:
        for c in cases:
            fp.write(json.dumps(c,allow_nan=False)+"\n")
            if c["priority"]:
                name=re.sub("[^A-Za-z0-9+-]","_",c["name"])
                (dest/(name+".json")).write_text(json.dumps(c,indent=2,allow_nan=False)+"\n")
    with gzip.open(dest/"historical_retrievals.jsonl.gz","wt") as fp:
        for rec in RECEIPTS:
            fp.write(json.dumps(rec)+"\n")
    print(f"Historical source checks: {len(cases)} dossiers; {len(targets)} primary titles",flush=True)


def atel_publication_index(cases):
    """Read indexed original telegram titles/abstracts, with a named-event control."""
    root="https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"
    def year(y):
        sql=f"SELECT bibcode,title,abstract FROM ref WHERE bibcode LIKE '{y}ATel%'"
        try:
            table=tap(root,sql)
            return y,table.fillna("").to_dict("records"),None
        except Exception as exc:
            return y,[],str(exc)[:180]
    covered,missing,docs={}, {},{}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for y,records,error in pool.map(year,range(1997,datetime.now(timezone.utc).year+1)):
            if error:
                missing[y]=error
            else:
                covered[y]=len(records)
            for r in records:
                docs[clean(r["bibcode"])]=r
    known=any("3811." in k and "063045" in clean(v["title"]) for k,v in docs.items())
    for c in cases:
        c["checks"]["SIMBAD ATel index"]={"control_3811":known,"records":len(docs),"year_counts":covered,
                                           "failed_years":missing,"scope":"only telegrams indexed by SIMBAD"}
    for code,r in docs.items():
        title,body=clean(r["title"]),clean(r["abstract"])
        n=re.search(r"ATel\.+(\d+)",code)
        if not n:
            continue
        url=ATEL+"?read="+str(int(n[1]))
        for c in cases:
            before=len(c["evidence"])
            add_publication(c,"original ATel indexed by SIMBAD",title,code[:4],url,body)
            if len(c["evidence"])>before:
                c["evidence"][-1]["bibcode"]=code
                c["evidence"][-1]["ads_url"]="https://ui.adsabs.harvard.edu/abs/"+quote(code)
                c["evidence"][-1]["abstract_available"]=bool(body)


def named_object_references(cases):
    """Keep neighbouring objects separate when a slew error circle is crowded."""
    from urllib.parse import urlencode
    root="https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"
    for c in cases:
        names=set(c["alias"])
        for finding in c["evidence"]:
            if finding["kind"]=="SIMBAD historical identifiers":
                names.add(finding["title"])
                names.update(finding.get("aliases",[]))
        quoted=",".join("'"+n.replace("'","''")+"'" for n in sorted(names) if n)
        try:
            objects=tap(root,"SELECT DISTINCT b.oid,b.main_id,b.otype,b.ra,b.dec FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE i.id IN ("+quoted+")")
            ids=",".join(str(int(x)) for x in objects.oid)
            if not ids:
                c["checks"]["named SIMBAD alternatives"]={"aliases":len(names),"objects":0}
                continue
            refs=tap(root,"SELECT h.oidref,r.bibcode,r.title,r.doi,r.abstract FROM has_ref AS h JOIN ref AS r ON h.oidbibref=r.oidbib WHERE h.oidref IN ("+ids+")").fillna("")
            status=[]
            for obj in objects.itertuples(index=False):
                name=clean(obj.main_id)
                d=degdist((obj.ra,obj.dec),(c["ra"],c["dec"]))
                subset=refs.loc[refs.oidref==obj.oid]
                status.append(dict(name=name,otype=clean(obj.otype),separation_arcsec=round(d,3),records=len(subset)))
                for r in subset.to_dict("records"):
                    title,abstract,code=clean(r["title"]),clean(r.get("abstract","")),clean(r["bibcode"])
                    if not (VARIABLE.search(title+" "+abstract) or "ATel" in code or "GCN" in code):
                        continue
                    url="https://doi.org/"+clean(r["doi"]) if clean(r.get("doi")) else "https://ui.adsabs.harvard.edu/abs/"+quote(code)
                    c["evidence"].append(dict(kind="SIMBAD object reference",title=title,bibcode=code,
                        date=code[:4],url=url,indexed_object=name,otype=clean(obj.otype),object_separation_arcsec=round(d,3),
                        interpretation="reference belongs to this named object; a neighbouring object's paper alone does not identify the X-ray emitter"))
                    if d<=c["gate"]:
                        add_publication(c,"paper for named positional alternative",title,code[:4],url,abstract)
            c["checks"]["named SIMBAD alternatives"]={"aliases":len(names),"objects":status}
        except Exception as exc:
            c["checks"]["named SIMBAD alternatives"]={"error":str(exc)[:200]}
        # Follow the reference lists for named NED galaxies, including sources
        # where SIMBAD has no central object. Do not transfer a positive entry
        # to a different galaxy in the error circle.
    for c in cases:
        hosts=[c.get("host","")]+[n for n in c["alias"] if n.startswith(("IRAS F","KUG ","2MASX ","WISEA J"))][:3]
        all_results=[]
        for name in dict.fromkeys(clean(n) for n in hosts if n):
            if name.isdigit():
                continue
            record={"name":name}
            try:
                params={"objname":name,"search_type":"Search","ref_extend":"no","start_year":1000,
                        "end_year":datetime.now(timezone.utc).year,"of":"xml_main"}
                url="https://ned.ipac.caltech.edu/cgi-bin/nph-refsearch"
                data=get(url,params,timeout=25)
                xml=parse(io.BytesIO(data))
                table=xml.get_first_table().to_table(use_names_over_ids=True).to_pandas().fillna("")
                fields={re.sub("[^a-z0-9]","",k.lower()):k for k in table}
                rc=next(k for n,k in fields.items() if n in ("refcode","bibcode","referencecode"))
                tc=next(k for n,k in fields.items() if n in ("title","reftitle","papertitle"))
                record.update(count=len(table),query=url+"?"+urlencode(params))
                for r in table.to_dict("records"):
                    title,code=clean(r[tc]),clean(r[rc])
                    if VARIABLE.search(title) or "ATel" in code or "GCN" in code:
                        link="https://ui.adsabs.harvard.edu/abs/"+quote(code)
                        c["evidence"].append(dict(kind="NED named object bibliography",title=title,url=link,
                            date=code[:4],bibcode=code,indexed_object=name,
                            interpretation="reference for named galaxy; exact flare correspondence requires event evidence"))
            except Exception as e:
                record["error"]=str(e)[:200] or type(e).__name__
            all_results.append(record)
        c["checks"]["named NED object bibliographies"]=all_results


if __name__=="__main__":
    cmd=argparse.ArgumentParser()
    cmd.add_argument("out",type=Path)
    cmd.add_argument("--historical",action="store_true",help="deepen an existing casebook with original XMMSL1 detections and primary paper text")
    args=cmd.parse_args()
    historical(args.out) if args.historical else main(args.out)

"""Find and audit public TNS catalogue replicas. No credentials or alert-only negatives."""
import gzip
import io
import json
import re
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urljoin, urlparse

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io.votable import parse
from bs4 import BeautifulSoup

# A southern and northern sample across seven years, with different event types.
CONTROLS = {"2016gkg":(23.5603,-29.4403), "2017gfo":(197.4504,-23.3815),
            "2018cow":(244.000925,22.268008), "2023ixf":(210.9107,54.31165)}
CONTROL_DATES={"2016gkg":"2016-09-20","2017gfo":"2017-08-17","2018cow":"2018-06-16","2023ixf":"2023-05-19"}
TAP = {
    "NOIRLab": "https://datalab.noirlab.edu/tap/sync",
    "HEASARC": "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync",
    "GAVO": "https://dc.g-vo.org/tap/sync",
    "IRSA": "https://irsa.ipac.caltech.edu/TAP/sync",
}
PUBLIC_DOCS = {
    "Lasair TNS schema": "https://lasair-ztf.lsst.ac.uk/schema/tns/",
    "ANTARES external catalogues": "https://antares.noirlab.edu/catalogs/",
    "ANTARES catalogue inventory": "https://api.antares.noirlab.edu/v1/catalogs",
    "SNAD crossmatch": "https://xmatch.snad.space/",
}

def key(s):
    return re.sub(r"[^a-z0-9]","",str(s).lower())

def info(audit, **data):
    audit.record({"service":"TNS-replica-audit",**data})

def short(err):
    return re.sub(r"\s+"," ",str(err))[:370]

def tap(audit,endpoint,query,maxrec):
    # Always use VOTable, which carries overflow and failure statuses.
    buf=audit.get("TNS-public-TAP",endpoint,params={"REQUEST":"doQuery","LANG":"ADQL",
                  "FORMAT":"votable","QUERY":query,"MAXREC":maxrec})
    v=parse(io.BytesIO(buf))
    flags=list(v.infos)+[i for resource in v.resources for i in resource.infos]
    bad=[i for i in flags if i.name=="QUERY_STATUS" and i.value in ("OVERFLOW","ERROR")]
    if bad:
        raise ValueError(f"{bad[0].value}: {bad[0].content}")
    return v.get_first_table().to_table(use_names_over_ids=True).to_pandas()

def completeness_claim(description):
    text=" ".join(description.lower().split())
    if re.search(r"classified (objects|sources|transients) only|classified subset|spectroscopically classified|training sample|ztf footprint|only matches|subset of",text):
        return False
    target=r"(?:\btns\b|transient name server)"
    census=r"(?:export|snapshot|mirror|copy|list|catalog(?:ue)?|database|census|archive|distribution)"
    return any(re.search(p,text) for p in (
        r"\b(?:full|complete|entire)\s+(?:public\s+|daily\s+|nightly\s+)?"+target+r"\s+(?:public\s+)?(?:objects|"+census+r")",
        r"\b(?:full|complete|entire)\s+(?:public\s+|daily\s+|nightly\s+)?"+census+r"\s+(?:of|from)\s+(?:the\s+)?"+target,
        r"\ball\s+(?:(?:public|publicly reported|reported|registered|named)\s+)?(?:objects|transients|entries|sources|records)\s+(?:registered |reported |listed |named )?(?:in|by|on|at|from)\s+(?:the\s+)?"+target,
        target+r"\s+(?:public\s+)?"+census+r"\s+in its entirety"))


def snapshot_time(description):
    # A most-recent event, HTTP header, DOI publication date, and local download
    # timestamp each mean something different from the coverage of a catalogue.
    match=re.search(r"(?:as of|complete through|snapshot (?:date|taken|of)|data (?:last )?updated(?: at| on)?)"
                    r"[^\d]{0,16}(20\d\d-[01]\d-[0-3]\d(?:[T ][0-2]\d:[0-5]\d(?::[0-5]\d)?)?)",
                    description,re.I)
    return match[1] if match else None

def inventory_tap(audit,item):
    service,url=item
    candidates=[]
    try:
        q=("SELECT table_name,description FROM TAP_SCHEMA.tables WHERE "
           "table_name LIKE '%tns%' OR table_name LIKE '%TNS%' OR "
           "description LIKE '%Transient Name Server%' OR description LIKE '%Transient Name Service%'")
        table=tap(audit,url,q,4000)
        table.columns=table.columns.str.lower()
        for _,r in table.iterrows():
            name=str(r.table_name).strip('"')
            desc=str(r.description)
            if name.lower().startswith("mydb.") or re.match(r"\s*select\s+(?:top|count)",desc,re.I):
                info(audit,archive=service,excluded_saved_query=name,reason="user query result / subset rather than archive distribution")
                continue
            if re.search(r"(?:^|[./_])tns(?:$|[./_])|tns_public_objects|transient.name.server",name+" "+desc,re.I):
                candidates.append({"kind":"tap","label":service+":"+name,"url":url,
                                   "table":name,"description":desc})
        info(audit,archive=service,tables=[[r["table"],r["description"][:380]] for r in candidates[:25]])
    except Exception as exc:
        info(audit,archive=service,error=short(exc))
    return candidates

def inventory_registry(audit):
    # The IVOA registry indexes publicly advertised data resources at other
    # archives as well as the individual TAP services checked above.
    endpoint="https://dc.g-vo.org/tap/sync"
    q="""SELECT DISTINCT t.ivoid,t.table_name,t.table_description,i.access_url
         FROM rr.res_table AS t JOIN rr.capability AS c ON t.ivoid=c.ivoid
         JOIN rr.interface AS i ON i.ivoid=c.ivoid AND i.cap_index=c.cap_index
         WHERE c.standard_id='ivo://ivoa.net/std/tap' AND
         (1=ivo_nocasematch(t.table_name,'%tns_public_objects%') OR
          1=ivo_nocasematch(t.table_description,'%transient name server%'))"""
    try:
        data=tap(audit,endpoint,q,1500)
        found=[]
        for _,r in data.iterrows():
            url=str(r.access_url).rstrip("/")
            host=urlparse(url).hostname or ""
            if urlparse(url).scheme!="https" or not (
                host.endswith((".edu",".org",".ac.uk",".esa.int",".aip.de",".unistra.fr",
                               ".in2p3.fr",".cam.ac.uk",".mpe.mpg.de")) and not re.search(r"^\d",host)):
                continue
            if not url.endswith("/sync"):
                url+="/sync"
            found.append(dict(kind="tap",label="IVOA "+str(r.ivoid),url=url,
                              table=str(r.table_name).strip('"'),description=str(r.table_description)))
        info(audit,archive="IVOA resource registry",distributions=[[c["table"],c["url"]] for c in found[:35]])
        return found
    except Exception as exc:
        info(audit,archive="IVOA resource registry",error=short(exc))
        return []

def public_json_distributions(doc,url,service):
    """Read explicit catalogue-file links in public provider metadata records."""
    results=[]
    def visit(obj):
        if isinstance(obj,list):
            for child in obj:
                visit(child)
        elif isinstance(obj,dict):
            parts={k:v for k,v in obj.items() if isinstance(v,(str,int,float))}
            label=" ".join(str(v) for v in parts.values())
            if re.search(r"tns_public_objects|transient name server",label,re.I):
                for k,v in parts.items():
                    if isinstance(v,str) and v.startswith("https://") and re.search(r"\.(?:parquet|csv(?:\.gz|\.zip)?)$",v,re.I):
                        results.append(dict(kind="download",label=service,url=v,description=label))
            for child in obj.values():
                if isinstance(child,(list,dict)):
                    visit(child)
    visit(doc)
    return results

def inventory_docs(audit,item):
    service,url=item
    try:
        raw=audit.get("TNS-provider-documentation",url)
        if raw.lstrip().startswith((b"{",b"[")):
            doc=json.loads(raw)
            direct=public_json_distributions(doc,url,service)
            # Use the actual self relation returned by the ANTARES inventory.
            # The relation may be advertised as HTTP; upgrade it to the same public TLS host.
            for item in doc.get("data",[]) if isinstance(doc,dict) else []:
                if not isinstance(item,dict) or item.get("attributes",{}).get("name")!="tns_public_objects":
                    continue
                ref=item.get("links",{}).get("self","")
                route=urlparse(ref)
                if route.hostname!="api.antares.noirlab.edu" or not re.fullmatch(r"/v1/catalogs/\d+",route.path):
                    continue
                catalog_url="https://api.antares.noirlab.edu"+route.path
                try:
                    payload=json.loads(audit.get("TNS-public-catalogue-record",catalog_url))
                    direct.extend(public_json_distributions(payload,catalog_url,service))
                    info(audit,provider=service,archive_record_url=catalog_url,
                         archive_record_metadata=json.dumps(payload)[:6000])
                except Exception as exc:
                    info(audit,provider=service,archive_record_url=catalog_url,error=short(exc))
            content=json.dumps(doc)
            match=re.search(r"tns_public_objects|transient name server",content,re.I)
            excerpt=content[max(match.start()-200,0):match.end()+800] if match else ""
            info(audit,provider=service,url=url,structured_metadata_excerpt=excerpt,
                 explicitly_linked_snapshots=[d["url"] for d in direct])
            return direct
        soup=BeautifulSoup(raw,"html.parser")
        plain=" ".join(soup.stripped_strings)
        # Follow only direct public downloads explicitly linked by the publisher.
        links=[]
        candidates=[]
        for a in soup.select("a[href]"):
            href=urljoin(url,a.get("href",""))
            host=urlparse(href).hostname or ""
            if "tns" not in (str(a)+" "+href).lower():
                continue
            links.append(href)
            if not (host==urlparse(url).hostname or host.endswith((".noirlab.edu",".snad.space",".lsst.ac.uk"))):
                continue
            if re.search(r"\.csv(?:\.zip|\.gz)?$|\.parquet$",href,re.I):
                scope_text=a.parent.get_text(" ",strip=True)
                candidates.append({"kind":"download","label":service,"url":href,
                                   "description":scope_text})
        m=re.search(r"transient name server|tns_public_objects|\bTNS\b",plain,re.I)
        snippet=plain[max(0,m.start()-120):m.end()+650] if m else plain[:360]
        info(audit,provider=service,url=url,documentation_excerpt=snippet,tns_links=links[:20])
        return candidates
    except Exception as exc:
        info(audit,provider=service,url=url,error=short(exc))
        return []

def inventory_zenodo(audit):
    # Only publisher-deposited public snapshots with TNS-specific distribution files.
    url="https://zenodo.org/api/records"
    try:
        data=json.loads(audit.get("TNS-archive-inventory",url,params={
                   "q":'"Transient Name Server"',"size":30}))
        candidates=[]
        for r in data["hits"]["hits"]:
            meta=r["metadata"]
            description=BeautifulSoup(meta.get("description",""),"html.parser").get_text(" ",strip=True)
            title=meta.get("title","")
            for f in r.get("files",[]):
                fn=f.get("key","")
                if not re.search(r"tns.*\.(csv|csv\.gz|csv\.zip|parquet)$",fn,re.I) or f.get("size",0)>350_000_000:
                    continue
                link=f.get("links",{}).get("self","")
                if (urlparse(link).hostname or "").lower()!="zenodo.org":
                    continue
                candidates.append({"kind":"download","label":"Zenodo "+str(r.get("doi",r.get("id"))),
                    "url":link,"description":title+" "+description,"declared_checksum":f.get("checksum","")})
        info(audit,archive="Zenodo",public_candidate_files=[c["label"] for c in candidates])
        return candidates
    except Exception as exc:
        info(audit,archive="Zenodo",error=short(exc))
        return []

def columns(names):
    d={key(v):v for v in names}
    wanted={"name":["tnsname","fullname","objname","name","objidname"],"ra":["radeg","ra","raj2000"],
            "dec":["decdeg","dec","declination","decl","dej2000"],"prefix":["nameprefix","prefix"],
            "type":["objtype","tnstype","classification","type"],
            "date":["discoverydateut","discoverydate","discoverydatetime","discoverytimestamp","discovered"],
            "id":["objid","tnsid","id"]}
    pick={k:next((d[s] for s in choices if s in d),None) for k,choices in wanted.items()}
    if not all(pick[k] for k in ("name","ra","dec")):
        raise ValueError("missing TNS name and coordinates in " + str(list(names))[:380])
    return {k:v for k,v in pick.items() if v is not None}

def sql_name(s):
    # VizieR-style embedded slashes are a single quoted table identifier.
    fields=[s] if "/" in s else s.split(".")
    return ".".join('"'+v.strip('"').replace('"','""')+'"' for v in fields)

def fetch_table(audit,source):
    cutoff=snapshot_time(source["description"])
    expected=None
    if source["kind"]=="tap":
        name=source["table"]
        esc=name.replace("'","''")
        meta=tap(audit,source["url"],f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{esc}' OR table_name='\"{esc}\"'",3000)
        use=columns(meta.iloc[:,0].str.strip('"'))
        table=sql_name(name)
        expected=int(tap(audit,source["url"],"SELECT COUNT(*) FROM "+table,2).iloc[0,0])
        if not 1<=expected<=2_000_000:
            raise ValueError(f"unsupported table row count {expected}")
        fields=",".join(sql_name(f) for f in use.values())
        df=tap(audit,source["url"],"SELECT "+fields+" FROM "+table,expected+1)
    else:
        data=audit.get("TNS-public-archive",source["url"])
        if source.get("declared_checksum","").startswith("md5:"):
            import hashlib
            if hashlib.md5(data).hexdigest()!=source["declared_checksum"].split(":")[1]:
                raise ValueError("depositor checksum mismatch")
        if data.startswith(b"PK"):
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                files=[n for n in z.namelist() if n.lower().endswith(".csv")]
                if len(files)!=1 or z.getinfo(files[0]).file_size>650_000_000:
                    raise ValueError("ambiguous/oversize CSV archive")
                data=z.read(files[0])
        if data.startswith(b"\x1f\x8b"):
            data=gzip.decompress(data)
        if data.startswith(b"PAR1"):
            df=pd.read_parquet(io.BytesIO(data))
        else:
            first,_,rest=data.partition(b"\n")
            timestamp=first.decode("utf8","replace").lstrip("\ufeff").strip('"\r ')
            if re.fullmatch(r"20\d\d-\d\d-\d\d[ T]\d\d:\d\d:\d\d(?:\.\d+)?(?: UTC|Z)?",timestamp):
                cutoff=timestamp
                data=rest
            df=pd.read_csv(io.BytesIO(data),dtype=str,low_memory=False)
        use=columns(df.columns)
    if expected is not None and len(df)!=expected:
        raise ValueError(f"query truncated: {len(df)} of {expected}")
    return df.rename(columns={v:k for k,v in use.items()})[list(use)],cutoff,expected

def coordinates(df):
    ra=pd.to_numeric(df.ra,errors="coerce").to_numpy(float)
    dec=pd.to_numeric(df.dec,errors="coerce").to_numpy(float)
    bad=~(np.isfinite(ra)&np.isfinite(dec))
    # Native TNS CSV uses sexagesimal hours for RA and sexagesimal degrees for Dec.
    if bad.any():
        ix=np.flatnonzero(bad)
        try:
            sky=SkyCoord(df.ra.iloc[ix].astype(str).to_numpy(),df.dec.iloc[ix].astype(str).to_numpy(),unit=(u.hourangle,u.deg))
            ra[ix],dec[ix]=sky.ra.deg,sky.dec.deg
        except (ValueError,TypeError):
            for i in ix:
                try:
                    c=SkyCoord(str(df.ra.iloc[i]),str(df.dec.iloc[i]),unit=(u.hourangle,u.deg))
                    ra[i],dec[i]=c.ra.deg,c.dec.deg
                except (ValueError,TypeError):
                    pass
    return ra,dec

def benchmark_names(names,ra,dec):
    """Check named independent control transients and sky locations in a candidate release."""
    outcomes={}
    for name,(x,y) in CONTROLS.items():
        ii=np.where(names==name)[0]
        if not len(ii):
            outcomes[name]=False
            continue
        sky=SkyCoord(ra[ii]*u.deg,dec[ii]*u.deg)
        outcomes[name]=bool((sky.separation(SkyCoord(x*u.deg,y*u.deg)).arcsec<=10.).any())
    return outcomes

def vet_snapshot(audit,src,points):
    frame,cutoff,expected=fetch_table(audit,src)
    ra,dec=coordinates(frame)
    raw=frame.name.fillna("").astype(str).str.strip()
    normalized=raw.str.lower().str.replace(r"^(?:at|sn)\s*","",regex=True).str.replace(" ","",regex=False)
    valid=(ra>=0)&(ra<360)&(np.abs(dec)<=90)&normalized.str.match(r"^\d{4}[a-z]+$").to_numpy()
    frame=frame.loc[valid].reset_index(drop=True)
    ra,dec,names=ra[valid],dec[valid],normalized.to_numpy()[valid]
    controls=benchmark_names(names,ra,dec)
    dated=pd.to_datetime(frame.get("date",pd.Series(dtype=str)),errors="coerce",utc=True)
    latest=dated.max() if len(dated) else pd.NaT
    if "prefix" in frame:
        prefixes=frame.prefix.fillna("").astype(str).str.strip().str.upper()
    else:
        prefixes=frame.name.fillna("").astype(str).str.extract(r"^\s*(AT|SN)",expand=False,flags=re.I).str.upper()
    kind=frame.get("type",pd.Series("",index=frame.index)).fillna("").astype(str).str.strip().str.lower()
    unnamed_class=(prefixes.eq("AT")&kind.isin(["","--","unknown","nan","unclassified","none"])).sum()
    census_utc=pd.to_datetime(cutoff,errors="coerce",utc=True) if cutoff else pd.NaT
    now=pd.Timestamp.now(tz="UTC")
    recent=bool(pd.notna(census_utc) and pd.Timedelta(0)<=now-census_utc<=pd.Timedelta(days=2))
    good_scope=completeness_claim(src["description"])
    # A snapshot that contains famous events may still contain exclusively
    # classified supernovae. Also require unclassified AT entries and a
    # trustworthy accounting of all rows before treating any empty cone as checked.
    unique=len(set(names))
    eligible_controls=[name for name,day in CONTROL_DATES.items() if
                       pd.notna(census_utc) and census_utc>=pd.Timestamp(day,tz="UTC")]
    covered=bool(good_scope and pd.notna(census_utc) and census_utc<=now and
                "2018cow" in eligible_controls and len(eligible_controls)>=3 and
                all(controls[name] for name in eligible_controls) and valid.all() and
                unique==len(frame) and unnamed_class>0 and pd.notna(latest) and
                latest>=census_utc-pd.Timedelta(days=10) and latest<=census_utc+pd.Timedelta(days=1) and
                (expected==len(frame) or (src["kind"]=="download" and bool(cutoff))))
    status=("ok" if recent else "complete_through_"+census_utc.strftime("%Y-%m-%d")) if covered else None
    info(audit,catalogue=src["label"],endpoint=src["url"],scope=src["description"][:1400],
         coverage_cutoff=cutoff,latest_discovery=str(latest),provider_count=expected,rows=len(valid),parsed=len(frame),
         unique_names=unique,unclassified_AT=int(unnamed_class),controls=controls,controls_before_cutoff=eligible_controls,
         verified_whole_census=covered,complete_through=cutoff if covered else None,current_within_48h=recent,vetting_status=status)
    # For positive-only records demand at least one accurate independent position
    # and positional validation of the other controls when present.
    if not any(controls.values()) or not points or not len(frame):
        return {},None
    p=SkyCoord([v["ra"] for v in points]*u.deg,[v["dec"] for v in points]*u.deg)
    c=SkyCoord(ra*u.deg,dec*u.deg)
    ii,jj,ss,_=search_around_sky(p,c,35*u.arcsec)
    hits=defaultdict(list)
    for i,j,sep in zip(ii,jj,ss.arcsec):
        ev=frame.iloc[j]
        hits[points[i]["k"]].append({"name":str(ev["name"]),"type":str(ev.get("type","")),
            "date":str(ev.get("date","")),"ra":ra[j],"dec":dec[j],"sep_arcsec":float(sep),
            "source_snapshot":src["url"],"provider":src["label"],"complete_through":cutoff})
    return hits,status

def tns_lookup(points,audit):
    """Return per-field historical associations and coverage after official TNS fallback."""
    candidates=[]
    with ThreadPoolExecutor(max_workers=5) as pool:
        futures=[pool.submit(inventory_tap,audit,i) for i in TAP.items()]
        futures+= [pool.submit(inventory_docs,audit,i) for i in PUBLIC_DOCS.items()]
        futures.extend([pool.submit(inventory_zenodo,audit),pool.submit(inventory_registry,audit)])
        for f in futures:
            candidates.extend(f.result())
    seen=set()
    out=defaultdict(list)
    checked_through=None
    candidates.sort(key=lambda c:not completeness_claim(c["description"]))
    for src in candidates[:10]:
        identity=(src["url"],src.get("table",""))
        if identity in seen:
            continue
        seen.add(identity)
        try:
            hits,coverage_status=vet_snapshot(audit,src,points)
            for k,h in hits.items():
                out[k].extend(h)
            if coverage_status=="ok":
                checked_through="ok"
                break
            if coverage_status and (checked_through is None or coverage_status>checked_through):
                checked_through=coverage_status
        except Exception as exc:
            info(audit,catalogue=src["label"],error=short(exc))
    status={p["k"]:(checked_through if checked_through else
                  "partial: public TNS archives contain a positional match, full current census unverified"
                  if out.get(p["k"]) else "incomplete: no verified full current public TNS replica") for p in points}
    info(audit,candidate_distributions=len(candidates),tested_distributions=len(seen),
         fields_with_matches=len(out),coverage_status=checked_through or "partial or unavailable")
    return out,status

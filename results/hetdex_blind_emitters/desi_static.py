"""Targeted DESI DR1 Redrock redshift lookup in the public iron HEALPix archive."""
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
import io
from html.parser import HTMLParser
import re
from urllib.parse import urljoin

from astropy import units as u
from astropy.io import fits
from astropy_healpix import HEALPix
import numpy as np
import requests

ROOT="https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/"
PORTAL="https://portal.nersc.gov/cfs/desi/public/dr1/spectro/redux/iron/healpix/"
RADIUS=1.25
PROGRAMS=("dark","bright","backup","other")
SURVEYS=("main","sv1","sv2","sv3","special","cmx")
TIMEOUT=(4,16)
ASU=("https://vizier.cds.unistra.fr/viz-bin/asu-tsv","https://vizier.cfa.harvard.edu/viz-bin/asu-tsv")


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.folders=set()

    def handle_starttag(self,tag,attrs):
        value=dict(attrs).get("href","")
        if tag=="a" and re.fullmatch(r"[A-Za-z0-9_-]+/",value):
            self.folders.add(value[:-1])


def listing(url):
    result=dict(url=url,attempts=[])
    for endpoint in (url,url.replace(ROOT,PORTAL)):
        try:
            r=requests.get(endpoint,timeout=TIMEOUT)
            item=dict(url=endpoint,http=r.status_code)
            result["attempts"].append(item)
            if r.status_code==404:
                return None,dict(result,status="404")
            r.raise_for_status()
            html=Links()
            html.feed(r.text)
            if not html.folders:
                raise ValueError("No recognized directory entries")
            return html.folders,dict(result,status="listed",sha256=sha256(r.content).hexdigest(),
                                    directories=sorted(html.folders))
        except (OSError,ValueError) as exc:
            result["attempts"].append(dict(url=endpoint,error=str(exc)[:160]))
    return None,dict(result,status="failed")


def inventory(pixels):
    """Prune only from positively identified Apache directory indexes; probe other paths directly."""
    names,root=listing(ROOT)
    surveys=sorted(names or SURVEYS)
    receipt=[root]
    complete=names is not None
    groupings=defaultdict(set)
    for pixel in pixels:
        groupings[pixel//100].add(pixel)
    todo=[]
    def children(s):
        return s,*listing(urljoin(ROOT,s+"/"))
    with ThreadPoolExecutor(max_workers=6) as pool:
        for survey,programs,meta in pool.map(children,surveys):
            receipt.append(meta)
            if meta["status"]=="404":
                continue
            if not programs:
                complete=False
            todo.extend((survey,q) for q in sorted(programs or PROGRAMS))
    def program_index(t):
        s,p=t
        return s,p,*listing(urljoin(ROOT,f"{s}/{p}/"))
    target=[]
    outcomes=defaultdict(Counter)
    with ThreadPoolExecutor(max_workers=7) as pool:
        for s,p,groups,meta in pool.map(program_index,todo):
            receipt.append(meta)
            for group,ids in groupings.items():
                if meta["status"]=="404" or (groups and str(group) not in groups):
                    for pix in ids:
                        outcomes[pix]["no_product_in_directory"]+=1
                else:
                    target.append((s,p,group))
    def pix_index(t):
        s,p,g=t
        return s,p,g,*listing(urljoin(ROOT,f"{s}/{p}/{g}/"))
    files=[]
    with ThreadPoolExecutor(max_workers=8) as pool:
        for s,p,g,ids,meta in pool.map(pix_index,target):
            receipt.append(meta)
            for pix in groupings[g]:
                if meta["status"]=="404" or (ids and str(pix) not in ids):
                    outcomes[pix]["no_product_in_directory"]+=1
                else:
                    files.append((s,p,pix))
    return files,outcomes,dict(complete_survey_program_enumeration=complete,receipts=receipt)


def redrock(file,points):
    from focus import RangeFile
    from vetting import separation
    s,p,pix=file
    url=urljoin(ROOT,f"{s}/{p}/{pix//100}/{pix}/redrock-{s}-{p}-{pix}.fits")
    receipt=dict(url=url,survey=s,program=p,healpix=pix,attempts=[])
    hits=defaultdict(list)
    # Most redrock files are small, consisting only of fit and fibre tables.
    # Use a single bounded range where possible, otherwise exact per-HDU ranges.
    size_limit=12*1024*1024
    for endpoint in (url,url.replace(ROOT,PORTAL)):
        try:
            with requests.get(endpoint,headers={"Range":f"bytes=0-{size_limit-1}","Accept-Encoding":"identity"},
                              stream=True,timeout=TIMEOUT) as result:
                if result.status_code==404:
                    receipt["attempts"].append(dict(url=endpoint,http=404))
                    return {},dict(receipt,status="404_no_published_file")
                result.raise_for_status()
                total=int(result.headers.get("Content-Length",0))
                span=result.headers.get("Content-Range","")
                if result.status_code==206:
                    match=re.fullmatch(r"bytes 0-([0-9]+)/([0-9]+)",span)
                    if not match:
                        raise ValueError("Invalid file range")
                    total=int(match[2])
                elif result.status_code!=200:
                    raise ValueError(f"Unexpected HTTP {result.status_code}")
                if total and total<=size_limit:
                    data=result.content
                    if len(data)!=total or not data.startswith(b"SIMPLE  ="):
                        raise ValueError("Invalid or truncated FITS payload")
                    src=io.BytesIO(data)
                    receipt["attempts"].append(dict(url=endpoint,http=result.status_code,bytes=len(data),sha256=sha256(data).hexdigest()))
                else:
                    src=None
            if src is None:
                src=RangeFile(endpoint)
            with src, fits.open(src,memmap=False) as hd:
                fm=hd["FIBERMAP"].data
                ras,decs=np.asarray(fm["TARGET_RA"],float),np.asarray(fm["TARGET_DEC"],float)
                decs_ok=np.isfinite(ras)&np.isfinite(decs)&(abs(decs)<=90)
                n_compared=int(decs_ok.sum())
                if not n_compared:
                    raise ValueError("No finite target positions in FIBERMAP")
                coords=[(point,np.asarray(separation(*point,ras,decs))) for point in points]
                matches=[(point,ds,np.flatnonzero(decs_ok & (ds<=RADIUS))) for point,ds in coords]
                if any(len(ix) for _,_,ix in matches):
                    rz=hd["REDSHIFTS"].data
                    rzid=np.asarray(rz["TARGETID"]).astype(np.int64)
                    for point,ds,ix in matches:
                        for j in ix:
                            objtype=str(fm["OBJTYPE"][j]).strip().upper() if "OBJTYPE" in fm.names else ""
                            if objtype in ("SKY","NON","BAD"):
                                continue
                            tid=int(fm["TARGETID"][j])
                            indices=np.flatnonzero(rzid==tid)
                            if len(indices)!=1:
                                raise ValueError("Nonunique redshift row for a spatial match")
                            rr=rz[indices[0]]
                            obj=dict(ra=str(ras[j]),dec=str(decs[j]),sep_arcsec=float(ds[j]),
                                 survey=s,prog=p,healpix=str(pix),specid=str(tid),targetid=str(tid),
                                 data_release="DESI-DR1",spectype=str(rr["SPECTYPE"]).strip(),
                                 redshift=str(float(rr["Z"])),redshift_warning=str(int(rr["ZWARN"])),
                                 deltachi2=float(rr["DELTACHI2"]),objtype=objtype,
                                 coadd_fiberstatus=int(fm["COADD_FIBERSTATUS"][j]) if "COADD_FIBERSTATUS" in fm.names else None,
                                 redrock_url=url,catalog="DESI_DR1_static")
                            if all(k in fm.names for k in ("MEAN_FIBER_RA","MEAN_FIBER_DEC")):
                                obj["mean_fibre_separation_arcsec"]=float(separation(*point,
                                      float(fm["MEAN_FIBER_RA"][j]),float(fm["MEAN_FIBER_DEC"][j])))
                            hits[point].append(obj)
                if hasattr(src,"receipts"):
                    receipt["attempts"].append(dict(url=endpoint,http=206,ranges=src.receipts))
            return hits,dict(receipt,status="searched",finite_fibermap_positions=n_compared)
        except (OSError,ValueError,KeyError,TypeError,IndexError) as exc:
            hits.clear()
            receipt["attempts"].append(dict(url=endpoint,error=str(exc)[:240]))
    return {},dict(receipt,status="unavailable")


def asu(point):
    """Single-cone CDS ASU interface; an unrecognized response cannot establish an empty cone."""
    ra,dec=point
    columns=("TargetID","RAICRS","DEICRS","Survey","Prog","Healpix","z","ZWARN","OType","delChi2")
    params={"-source":"V/161/zcatdr1","-c":f"{ra:.9f} {dec:+.9f}","-c.r":str(RADIUS),"-c.u":"arcsec",
            "-out":",".join(columns),"-out.max":"500"}
    trials=[]
    for url in ASU:
        try:
            resp=requests.get(url,params=params,timeout=TIMEOUT)
            resp.raise_for_status()
            data=resp.content.decode("utf-8-sig")
            body=[l for l in data.splitlines() if l.strip() and not l.startswith("#")]
            if not body or not set(columns).issubset(set(x.strip() for x in body[0].split("\t"))):
                # ASU can omit even the header on an empty query. Only accept an
                # explicit QUERY_STATUS OK with a real zero-row count or message.
                if ("V/161" in data and "QUERY_STATUS" in data and
                    re.search(r"(?:0\s+(?:records|rows|sources)|no\s+(?:records|rows|sources)\s+(?:found|selected))",data,re.I)):
                    body=["\t".join(columns)]
                else:
                    raise ValueError("Unrecognized CDS table or empty-cone response")
            result=[]
            from vetting import separation
            for m in csv.DictReader(body,delimiter="\t"):
                m={k.strip():v.strip() for k,v in m.items() if k and v is not None}
                if not re.fullmatch(r"-?\d+",m.get("TargetID","")):
                    continue
                ds=float(separation(ra,dec,float(m["RAICRS"]),float(m["DEICRS"])))
                if ds<=RADIUS:
                    result.append(dict(ra=m["RAICRS"],dec=m["DEICRS"],specid=m["TargetID"],targetid=m["TargetID"],
                              survey=m["Survey"],prog=m["Prog"],healpix=m["Healpix"],redshift=m["z"],
                              redshift_warning=m["ZWARN"],spectype=m["OType"],deltachi2=m["delChi2"],
                              data_release="DESI-DR1",sep_arcsec=ds,catalog="CDS_ASU_DESI_DR1"))
            trials.append(dict(url=resp.url,sha256=sha256(resp.content).hexdigest(),bytes=len(resp.content),status="queried"))
            return result,dict(status="queried",attempts=trials)
        except (OSError,ValueError,KeyError,TypeError) as exc:
            trials.append(dict(url=url,error=str(exc)[:200]))
    return [],dict(status="unavailable",attempts=trials)


def search(positions):
    from focus import coadd_url
    positions=sorted(set(positions))
    if not positions:
        return {},dict(catalog=[],redrock=[],asu=[])
    hp=HEALPix(nside=64,order="nested")
    per_pixel=defaultdict(list)
    for ra,dec in positions:
        for pix in hp.cone_search_lonlat(ra*u.deg,dec*u.deg,RADIUS*u.arcsec):
            per_pixel[int(pix)].append((ra,dec))
    files,outcomes,tree=inventory(per_pixel)
    bypoint={p:dict(matches=[],files={},expected_HEALPix=[]) for p in positions}
    receipts=[]
    def get(f):
        return f,redrock(f,per_pixel[f[2]])
    with ThreadPoolExecutor(max_workers=7) as pool:
        for (_,_,pix),(hits,receipt) in pool.map(get,files):
            receipts.append(receipt)
            outcomes[pix][receipt["status"]]+=1
            for q in per_pixel[pix]:
                bypoint[q]["matches"].extend(hits.get(q,[]))
    for pix,points in per_pixel.items():
        for q in points:
            bypoint[q]["expected_HEALPix"].append(pix)
            bypoint[q]["files"]=dict(Counter(bypoint[q]["files"])+outcomes[pix])
    asu_receipts=[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        for q,(rows,receipt) in zip(positions,pool.map(asu,positions)):
            rec=bypoint[q]
            files=rec["files"]
            complete=tree["complete_survey_program_enumeration"] and not files.get("unavailable")
            if complete:
                static_status="searched" if files.get("searched") else "no_published_pixel_product"
            else:
                static_status="partial_or_unavailable"
            rec["static_status"]=static_status
            rec["ASU_status"]=receipt["status"]
            asu_receipts.append(dict(ra=q[0],dec=q[1],**receipt))
            # Prefer native FITS metadata when both catalogues report the same
            # target/coadd. Preserve both query traces for independent checking.
            matches={}
            for m in rows+rec["matches"]:
                key=(coadd_url(m),m["specid"])
                matches[key]=m
            rec["matches"]=list(matches.values())
            rec["status"]=("matched" if matches else "no_desi_dr1_counterpart" if
                    (complete or receipt["status"]=="queried") else "partial_or_unavailable")
    return bypoint,dict(catalog=tree,redrock=receipts,asu=asu_receipts)

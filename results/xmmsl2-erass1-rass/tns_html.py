"""Normal, anonymous public TNS HTML search, with catalogue and cone controls."""
from collections import defaultdict
from datetime import datetime, timezone
import re
import time
from urllib.parse import parse_qs, urljoin, urlparse

import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from bs4 import BeautifulSoup

BASE = "https://www.wis-tns.org/search"
RADIUS = 35.
PAGE_SIZE = 50
CALIBRATORS = [("2018cow",244.000925,22.268008), ("2017gfo",197.4504,-23.3815),
               ("2023ixf",210.9107,54.31165)]
ZERO = re.compile(r"\b(?:no\s+(?:(?:matching|search)\s+)?(?:results?|objects?|records?|entries|items?|transients?)|"
                  r"(?:results?|objects?|records?)\s*(?:found)?\s*[:=]\s*0\b|"
                  r"\b0\s+(?:results?|objects?|records?|matches)\b|no data available in table)\b",re.I)
OBJECT = re.compile(r"/object/((?:1|2)\d{3}[a-z]+)(?:\b|/)",re.I)


def data_coordinate(ra,dec):
    ra,dec=str(ra).strip(),str(dec).strip().replace("−","-")
    try:
        r,d=float(ra),float(dec)
        if 0<=r<360 and -90<=d<=90:
            return SkyCoord(r*u.deg,d*u.deg)
    except ValueError:
        pass
    return SkyCoord(ra,dec,unit=(u.hourangle,u.deg))


def shortened(err):
    return " ".join(str(err).split())[:370]


def text(node):
    return " ".join(node.stripped_strings)


def header_index(headers, kind):
    for i,s in enumerate(headers):
        word=re.sub(r"[^a-z0-9]","",s.lower())
        if kind=="name" and word in ("name","tnsname","iauname","objectname"):
            return i
        if kind=="ra" and (word in ("ra","raj2000","raicrs","radeg") or word.startswith("rightascension")):
            return i
        if kind=="dec" and (word in ("dec","decl","decj2000","declj2000","decicrs","decdeg","declination") or word.startswith("declination")):
            return i
        if kind=="type" and word in ("type","objtype","objecttype","classification","spectype"):
            return i
        if kind=="date" and word.startswith(("discoverydate","discoverytime","discoveryutc")):
            return i
        if kind=="reported" and word.startswith(("reportingdate","reportdate","publictimestamp","reportingtime")):
            return i
    return None


def parse_page(raw):
    doc=BeautifulSoup(raw,"html.parser")
    if not doc.title or "transient name server" not in text(doc.title).lower():
        raise ValueError("Response is not a public Transient Name Server page")
    table=cols=headers=None
    for t in doc.select("table"):
        head=t.select_one("thead tr") or next((tr for tr in t.select("tr") if tr.find("th")),None)
        if not head:
            continue
        hdr=[text(c) for c in head.find_all(["th","td"],recursive=False)]
        idx={key:header_index(hdr,key) for key in ("name","ra","dec","type","date","reported")}
        if all(idx[k] is not None for k in ("name","ra","dec")):
            table,cols,headers=t,idx,hdr
            break
    if table is None:
        raise ValueError("No genuine object-results table with name, RA and Dec headings")
    rows=[]
    for tr in table.find_all("tr"):
        tds=tr.find_all("td",recursive=False)
        if not tds:
            continue
        if len(tds)==1 and ZERO.search(text(tds[0])):
            continue
        if len(tds)!=len(headers):
            raise ValueError("Unrecognized or incomplete public result-table row")
        named=tds[cols["name"]]
        links=[(a,OBJECT.search(urlparse(a.get("href","")).path)) for a in named.select("a[href]")]
        matches=[(a,m) for a,m in links if m]
        if len(matches)!=1:
            raise ValueError("Result row lacks a unique official object link")
        a,m=matches[0]
        pos=data_coordinate(text(tds[cols["ra"]]),text(tds[cols["dec"]]))
        rows.append({"name":text(a),"tns_slug":m[1].lower(),"ra":pos.ra.deg,"dec":pos.dec.deg,
                     "type":text(tds[cols["type"]]) if cols["type"] is not None else "",
                     "date":text(tds[cols["date"]]) if cols["date"] is not None else "",
                     "report_date":text(tds[cols["reported"]]) if cols["reported"] is not None else "",
                     "url":urljoin(BASE,a["href"])})
    # No-records copy may appear immediately above the table.
    container=table.find_parent("main") or table.parent
    empty_notice=bool(ZERO.search(text(container))) if not rows else False
    if not rows and not empty_notice:
        raise ValueError("Empty table without an explicit zero-result message")
    nexts=[]
    for a in doc.select('a[rel~="next"], li.pager__item--next a, li.pager-next a, a[title*="next page" i], ul.pager a, nav.pager a, .pagination a, a[href*="page="]'):
        href=a.get("href","")
        if not href:
            continue
        label=" ".join((text(a),str(a.get("title","")),str(a.get("aria-label","")))).casefold()
        classes=" ".join(a.parent.get("class",[])+a.get("class",[])).lower()
        if "last" in (label+" "+classes):
            continue
        if "next" in (label+" "+classes) or "next" in a.get("rel",[]) or text(a) in (">","›","»"):
            url=urljoin(BASE,href)
            p=urlparse(url)
            if p.scheme!="https" or p.hostname!="www.wis-tns.org" or p.path.rstrip("/")!="/search":
                raise ValueError("Foreign/invalid next-page URL")
            nexts.append(url)
    if len(set(nexts))>1:
        raise ValueError("Conflicting next-page links")
    # Detect pages whose pager shape is unsupported rather than dropping rows.
    pager=doc.select_one("ul.pager, nav.pager, .pager__items, .pagination")
    numeric_targets=set()
    if pager:
        for a in pager.select("a[href]"):
            q=parse_qs(urlparse(a["href"]).query)
            for v in q.get("page",[]):
                if v.isdigit():
                    numeric_targets.add(int(v))
    return rows,(nexts[0] if nexts else None),empty_notice,numeric_targets


def tns_html_lookup(points,audit):
    """Return hit and status dictionaries keyed with the calling vetter's 'k'."""
    hits,status=defaultdict(list),{}
    last_fetch=0.
    def announce(**kw):
        audit.record({"service":"TNS-public-HTML-audit",**kw})
    def request(params):
        nonlocal last_fetch
        time.sleep(max(0.,.5-(time.monotonic()-last_fetch)))
        last_fetch=time.monotonic()
        # Normal submitted search fields; no format/export, browser emulation,
        # alternate identity, authentication or challenge-handling mechanism.
        return parse_page(audit.get("TNS-public-HTML",BASE,params=params))
    def cone(x,y):
        pos=SkyCoord(x*u.deg,y*u.deg)
        params={"ra":pos.ra.to_string(unit=u.hourangle,sep=":",precision=3,pad=True),
                "decl":pos.dec.to_string(sep=":",alwayssign=True,precision=2,pad=True),
                "radius":str(RADIUS),"coords_unit":"arcsec","num_page":str(PAGE_SIZE)}
        known={}
        page=0
        while True:
            rows,link,empty,visible=request(dict(params,page=page))
            for r in rows:
                sep=pos.separation(SkyCoord(r["ra"]*u.deg,r["dec"]*u.deg)).arcsec
                if sep>RADIUS+.2:
                    raise ValueError("An official search result falls outside the requested cone")
                if r["tns_slug"] in known:
                    raise ValueError("Object repeated across pages of one cone")
                known[r["tns_slug"]]=dict(r,sep_arcsec=float(sep),source="official TNS public HTML search")
            if not link:
                if any(n>page for n in visible) or len(rows)>=PAGE_SIZE:
                    raise ValueError("Additional or potentially truncated result page not validated")
                break
            nextnum=parse_qs(urlparse(link).query).get("page",[""])[0]
            if nextnum!=str(page+1) or page>=99:
                raise ValueError("Unverified TNS cone pagination")
            # Submit original verified coordinate filters on each numbered page.
            page+=1
        return list(known.values()),page+1,empty
    try:
        # The public listing validates page numbering, live discovery/report dates,
        # classified and unclassified content, and at least one recent field.
        first,link,_,_=request({"num_page":str(PAGE_SIZE)})
        if not first or not link or parse_qs(urlparse(link).query).get("page")!=["1"]:
            raise ValueError("Public catalogue first page or pagination not verified")
        second,_,_,_=request({"num_page":str(PAGE_SIZE),"page":"1"})
        if not second or not ({r["tns_slug"] for r in second}-{r["tns_slug"] for r in first}):
            raise ValueError("Public catalogue second page did not advance")
        dated=[]
        all_rows=first+second
        for r in all_rows:
            for tag in ("report_date","date"):
                stamp=pd.to_datetime(r[tag],utc=True,errors="coerce")
                if pd.notna(stamp):
                    dated.append((stamp,r))
        if not dated:
            raise ValueError("No dated public listing records")
        most_recent,newest=max(dated,key=lambda item:item[0])
        age=pd.Timestamp.now(tz="UTC")-most_recent
        if not -pd.Timedelta(days=1)<=age<=pd.Timedelta(days=30):
            raise ValueError("Public listing did not establish recent coverage")
        # Prefer an explicitly unclassified AT from the two listing pages.
        def unclassified(r):
            return (r["name"].strip().upper().startswith("AT") and
                    r["type"].strip().lower() in ("","-","--","unknown","unclassified","n/a"))
        recent_example=next((r for r in first if unclassified(r)),newest)
        verified=[]
        for label,x,y in [*CALIBRATORS,(recent_example["tns_slug"],recent_example["ra"],recent_example["dec"])]:
            rows,_,_=cone(x,y)
            if not any(r["tns_slug"]==label and r["sep_arcsec"]<10 for r in rows):
                raise ValueError("Current public search failed named control "+label)
            verified.append(label)
        off=SkyCoord(CALIBRATORS[0][1]*u.deg,CALIBRATORS[0][2]*u.deg).directional_offset_by(0*u.deg,105*u.arcsec)
        rows,_,empty=cone(off.ra.deg,off.dec.deg)
        if any(r["tns_slug"]=="2018cow" for r in rows):
            raise ValueError("Out-of-cone AT 2018cow was returned")
        announce(utc=datetime.now(timezone.utc).isoformat(),archive_url=BASE,
                 positive_controls=verified,most_recent_public_record=most_recent.isoformat(),
                 recent_unclassified_AT=recent_example["tns_slug"] if unclassified(recent_example) else None,
                 out_of_cone_control={"offset_arcsec":105,"matches":len(rows),"explicit_zero":empty},
                 first_page=len(first),second_page=len(second))
    except Exception as exc:
        state="unverified official HTML search: "+shortened(exc)
        announce(error=state)
        return hits,{p["k"]:state for p in points}
    locations=defaultdict(list)
    for p in points:
        locations[(p["ra"],p["dec"])].append(p["k"])
    stop=""
    repeated_errors=0
    for (x,y),ks in locations.items():
        if stop:
            for k in ks:
                status[k]=stop
            continue
        try:
            found,pages,empty=cone(x,y)
            state="ok"
            repeated_errors=0
            announce(ra=x,dec=y,radius_arcsec=RADIUS,pages=pages,objects=len(found),
                     explicit_zero=empty,utc=datetime.now(timezone.utc).isoformat())
        except Exception as exc:
            found=[]
            state="public HTML cone unchecked: "+shortened(exc)
            repeated_errors+=1
            if isinstance(exc,requests.HTTPError) and exc.response.status_code in (401,403,429):
                stop="unqueried after official public search access/rate block"
            elif repeated_errors>=3:
                stop="unqueried after repeated invalid/unavailable official pages"
            announce(ra=x,dec=y,error=state)
        for k in ks:
            hits[k].extend(found)
            status[k]=state
    return hits,status

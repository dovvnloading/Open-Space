"""Cross-check public TDE compilations and preserve per-host bibliographic searches."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile
import zipfile

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord, search_around_sky

AGENT={"User-Agent":"Open-Space archival X-ray literature checks"}
RELEVANT=re.compile(r"X.ray|tidal|TDE|transien|flar|variab|outburst|switch|erupt|changing.look|high.state|low.state|turn.off|turn.on",re.I)


def get(url,**kwargs):
    headers=dict(AGENT)
    if url.startswith("https://api.github.com/") and os.environ.get("GH_TOKEN"):
        headers["Authorization"]="Bearer "+os.environ["GH_TOKEN"]
    r=requests.get(url,headers=headers,timeout=(12,60),**kwargs)
    r.raise_for_status()
    return r


def digest(b):
    return hashlib.sha256(b).hexdigest()


def designation(s):
    return re.sub(r"[^a-z0-9]","",str(s).lower())


def coord(ra,dec,hour=False):
    a,d=str(ra),str(dec)
    if hour or re.search(r"[:h\s]",a.strip()):
        c=SkyCoord(a,d,unit=(u.hourangle,u.deg))
        return c.ra.deg,c.dec.deg
    return float(a),float(d)


def read_event(data,key):
    record=data.get(key,data)
    if isinstance(record,list):
        return []
    if "coordinate" in record:
        # OTTER's published JSON schema.
        alias_block=record.get("name",{})
        name=alias_block.get("default_name",key) if isinstance(alias_block,dict) else str(alias_block)
        result=[]
        for c in record["coordinate"]:
            if "ra" not in c or "dec" not in c:
                continue
            x,y=coord(c["ra"],c["dec"],"hour" in str(c.get("ra_units","")))
            result.append((name,x,y))
        return result
    if not isinstance(record.get("ra"),list) or not isinstance(record.get("dec"),list):
        return []
    av=record["ra"][0]
    dv=record["dec"][0]
    x,y=coord(av["value"],dv["value"],"hour" in str(av.get("u_value","")))
    return [(record.get("name",key),x,y)]


def catalog_snapshots(out):
    objects=[]
    receipts=[]
    for org in ("astrocatalogs","astro-otter"):
        try:
            repos=[]
            for page in (1,2,3):
                r=get(f"https://api.github.com/orgs/{org}/repos",params={"per_page":100,"page":page})
                part=r.json()
                repos.extend(part)
                if len(part)<100:
                    break
            selected=[p for p in repos if re.fullmatch(r"tde-\d+[-\d]*|tidaldisruptions-\d+[-\d]*|otter-data|otterdb",p["name"],re.I)]
            receipts.append(dict(organization=org,matched_repositories=[p["full_name"] for p in selected]))
        except Exception as e:
            receipts.append(dict(organization=org,error=str(e)))
            continue
        for p in selected:
            rec=dict(repository=p["full_name"])
            try:
                sha=get(p["url"]+"/commits/"+p["default_branch"]).json()["sha"]
                url="https://codeload.github.com/"+p["full_name"]+"/zip/"+sha
                rec.update(url=url,sha=sha)
                with get(url,stream=True) as response,tempfile.SpooledTemporaryFile(max_size=8<<20) as fp:
                    h=hashlib.sha256()
                    for b in response.iter_content(1<<20):
                        fp.write(b)
                        h.update(b)
                    rec.update(bytes=fp.tell(),sha256=h.hexdigest())
                    fp.seek(0)
                    ngood,nschema,nfail=0,0,0
                    with zipfile.ZipFile(fp) as z:
                        for f in z.infolist():
                            if not f.filename.lower().endswith((".json",".json.gz")):
                                continue
                            try:
                                b=z.read(f)
                                raw=gzip.decompress(b) if f.filename.endswith(".gz") else b
                                data=json.loads(raw)
                                if not isinstance(data,dict):
                                    continue
                                key=Path(f.filename).name.split(".json")[0]
                                rows=read_event(data,key)
                                if not rows:
                                    nschema+=1
                                for name,x,y in rows:
                                    if np.isfinite([x,y]).all() and 0<=x<360 and -90<=y<=90:
                                        objects.append({"name":name,"ra":x,"dec":y,"repository":p["full_name"],"member":f.filename})
                                        ngood+=1
                            except (KeyError,ValueError,TypeError,IndexError):
                                nfail+=1
                    rec.update(positions=ngood,unrecognized_members=nschema,unparsed_members=nfail)
            except Exception as e:
                rec["error"]=str(e)
            receipts.append(rec)
    (out/"tde_inputs.json").write_text(json.dumps(receipts,indent=2)+"\n")
    pd.DataFrame(objects).to_csv(out/"tde_positions.csv.gz",index=False)
    return objects


def prior_matches(leads,events,out):
    columns=["detection_uid","source_name","published_name","repository","member","separation_arcsec"]
    rows=[]
    if len(leads) and len(events):
        src=SkyCoord(leads.ra.to_numpy()*u.deg,leads.dec.to_numpy()*u.deg)
        targets=SkyCoord([r["ra"] for r in events]*u.deg,[r["dec"] for r in events]*u.deg)
        si,ei,sep,_=search_around_sky(src,targets,45*u.arcsec)
        for i,j,d in zip(si,ei,sep.arcsec):
            event=events[j]
            rows.append([leads.iloc[i].detection_uid,leads.iloc[i]["name"],event["name"],event["repository"],event["member"],d])
    pd.DataFrame(rows,columns=columns).to_csv(out/"published_tde_matches.csv",index=False)


def find_papers(alias):
    # An exact alias match in the title/abstract is evidence. A high-ranking
    # search result that merely mentions the same instrument is insufficient.
    url="https://api.crossref.org/works"
    ident=designation(alias)
    rec={"alias":alias,"url":url,"query":"exact alias in title/abstract","hits":[]}
    try:
        r=get(url,params={"query.bibliographic":alias,"rows":30,"select":"DOI,title,abstract,published"})
        rec.update(sha256=digest(r.content),effective_url=r.url)
        items=r.json()["message"]["items"]
        for item in items:
            title=" ".join(item.get("title",[]))
            abstract=re.sub(r"<[^>]+>"," ",item.get("abstract",""))
            if ident in designation(title+" "+abstract):
                rec["hits"].append({"doi":item.get("DOI"),"title":title,
                                    "date":item.get("published"),"event_terms":bool(RELEVANT.search(title+" "+abstract))})
        rec["returned"]=len(items)
    except Exception as e:
        rec["error"]=str(e)
    return rec


def main(out):
    dest=out/"event_history"
    dest.mkdir(exist_ok=True)
    leads=pd.read_csv(out/"input.csv",low_memory=False).fillna("")
    vet=pd.read_csv(out/"vetting/vetting.csv",low_memory=False).fillna("")
    events=catalog_snapshots(dest)
    prior_matches(leads,events,dest)
    joined=leads.merge(vet,on="id",suffixes=("","_vet"),validate="one_to_one")
    alias_ids={}
    for _,row in joined.iterrows():
        if row.vet_status=="known_or_stellar":
            continue
        # Concentrate bibliographic requests on named extragalactic/compact
        # counterparts plus every source that passed the main identification.
        types=str(row.get("xmm_xtyp",""))+" "+str(row.get("host_type",""))
        if not (row.vet_status=="candidate" or re.search(r"gal|AGN|Sy[12]|QSO|LIN|BL.?Lac|XRB|HXB|LXB|ULX|^G\b",types,re.I)):
            continue
        for key in ("host_id","xmm_counterpart","name"):
            alias=str(row.get(key,"")).strip().strip("\0")
            if len(alias)<4 or alias.isdigit() or re.search(r"^Gaia|^WISE",alias,re.I):
                continue
            alias_ids.setdefault(alias,set()).add(row.id)
    with (dest/"source_bibliography.jsonl").open("w") as log,ThreadPoolExecutor(max_workers=2) as pool:
        for r in pool.map(find_papers,alias_ids):
            r["detection_uids"]=sorted(alias_ids[r["alias"]])
            log.write(json.dumps(r)+"\n")
    print("History check:",len(events),"published TDE coordinates and",len(alias_ids),"source aliases",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("out",type=Path)
    main(parser.parse_args().out)

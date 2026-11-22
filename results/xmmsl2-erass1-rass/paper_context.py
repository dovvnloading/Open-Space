"""Recover original electronic rows and documentation behind the final literature flags."""
import gzip
import hashlib
import io
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astroquery.vizier import VizierClass
from bs4 import BeautifulSoup

from vet import hash_table, parse_coord, positional_columns

LI = "J/MNRAS/512/3858"
ERO = "J/A+A/712/A171"
SEYFERT = "XMMSL2-total-69754"
PRIOR = {"erass3-h":253,"dr2hc":250,"dr2hg":9330,"df2hl":123}
FOCUS = {"XMMSL2-total-69754","XMMSL2-total-10038","XMMSL2-total-17748","XMMSL2-total-25333"}


def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str)+"\n")


def readme(cat,dest):
    urls=["https://cdsarc.cds.unistra.fr/ftp/"+cat+"/ReadMe",
          "https://vizier.cfa.harvard.edu/viz-bin/ReadMe/"+cat+"?format=html",
          "https://cdsarc.u-strasbg.fr/ftp/"+cat+"/ReadMe"]
    logs=[]
    for url in urls:
        try:
            r=requests.get(url,timeout=(8,40))
            r.raise_for_status()
            raw=r.content
            txt=raw.decode("utf-8","replace")
            if "<html" in txt[:1000].lower():
                txt=BeautifulSoup(raw,"html.parser").get_text("\n",strip=True)
            if "Byte-by-byte" not in txt and "Description of file" not in txt and "File Summary" not in txt:
                raise ValueError("no CDS catalogue documentation")
            dest.write_text(txt)
            return {"url":url,"sha256":hashlib.sha256(raw).hexdigest(),"saved":dest.name}
        except Exception as exc:
            logs.append([url,str(exc)[:280]])
    return {"failed":logs}


def get_table(name,center=None):
    errors=[]
    for host in ("vizier.cfa.harvard.edu","vizier.cds.unistra.fr"):
        try:
            viz=VizierClass(columns=["**","_RAJ2000","_DEJ2000"],row_limit=2000,vizier_server=host)
            viz.TIMEOUT=45
            res=(viz.query_region(center,radius=120*u.arcsec,catalog=name) if center is not None else viz.get_catalogs(name))
            if len(res):
                tab=res[0]
                if len(tab)>=2000:
                    raise ValueError("query truncated")
                return tab,host
            if center is None:
                raise ValueError("empty full table")
            # Retry the alternate host before reporting an empty cone.
            errors.append(host+" returned no table")
        except Exception as exc:
            errors.append(host+": "+str(exc)[:250])
    return None," | ".join(errors)


def row_values(tab,row):
    return {k:str(row[k]).strip() for k in tab.colnames if not np.ma.is_masked(row[k])}


def metadata(tab):
    return {k:{"description":str(tab[k].description or ""),"unit":str(tab[k].unit or "")} for k in tab.colnames}


def position(tab,row):
    ra,dec=positional_columns(tab)
    return SkyCoord(*parse_coord(row[ra],row[dec]),unit="deg")


def fink_snapshot(events):
    # Reopen the exact bytes already attempted by the primary run. A latest
    # event date in a small subset cannot establish a TNS completeness date.
    records=[]
    for e in events:
        if e.get("service")!="TNS-mirror" or e.get("status")!=200 or "parquet" not in e.get("url",""):
            continue
        item={"url":e["url"],"recorded_sha256":e["sha256"],"TNS_completeness_date":None}
        try:
            r=requests.get(e["url"],timeout=(10,60)); r.raise_for_status()
            digest=hashlib.sha256(r.content).hexdigest()
            if digest!=e["sha256"]:
                raise ValueError("recorded Fink file checksum differs")
            tab=pd.read_parquet(io.BytesIO(r.content))
            item.update(rows=len(tab),columns=list(tab),sha256=digest)
            names=[k for k in tab if re.search(r"name|object|designation",str(k),re.I)]
            if names:
                norm=tab[names].astype(str).apply(lambda s:s.str.lower().str.replace(" ","",regex=False))
                item["named_control_occurrences"]={k:int(norm.apply(lambda s:s.str.contains(k,regex=False)).any(axis=1).sum())
                                                for k in ("2018cow","2017gfo","2016gkg","2023ixf")}
            date_bounds={}
            for k in tab:
                if "date" in str(k).lower() and not pd.api.types.is_numeric_dtype(tab[k]):
                    dates=pd.to_datetime(tab[k],errors="coerce",utc=True)
                    if dates.notna().any():
                        date_bounds[k]=[str(dates.min()),str(dates.max())]
            item["event_or_record_date_ranges"]=date_bounds
        except Exception as exc:
            item["error"]=str(exc)[:320]
        records.append(item)
    return records


def run(out):
    dest=out/"paper_context"
    dest.mkdir(exist_ok=True,parents=True)
    src=pd.read_csv(out/"input.csv",low_memory=False)
    vet=pd.read_csv(out/"vetting/vetting.csv",low_memory=False).set_index("id")
    row=src.set_index("id").loc[SEYFERT]
    centers={k:SkyCoord(float(row[k+"_ra"])*u.deg,float(row[k+"_dec"])*u.deg)
             for k in ("xmm","ero","rass") if pd.notna(row[k+"_ra"])}
    summary={"position_degrees":{k:[v.ra.deg,v.dec.deg] for k,v in centers.items()}}
    read_for=[(LI,dest/"li_readme.txt"),(ERO,dest/"erass_readme.txt")]
    with ThreadPoolExecutor(max_workers=4) as pool:
        docs={cat:pool.submit(readme,cat,path) for cat,path in read_for}
        tabs={name:pool.submit(get_table,ERO+"/"+name,centers["ero"]) for name in
              ("erass3-h","erass3-m","dr2hc","dr2hg","df2hl","dr2mc","dr2mg","dr2ml")}
        li=pool.submit(get_table,LI+"/table")
        summary["documentation"]={cat:f.result() for cat,f in docs.items()}
        for name,task in tabs.items():
            tab,server=task.result()
            item={"catalog":ERO,"table":name,"server_or_status":server,"query_radius_arcsec":120,
                  "previous_returned_row":PRIOR.get(name),"rows":[]}
            if tab is not None:
                item.update(columns=metadata(tab),rows_sha256=hash_table(tab),cone_rows=len(tab))
                for r in tab:
                    pos=position(tab,r)
                    d={k:round(float(pos.separation(v).arcsec),3) for k,v in centers.items()}
                    if min(d.values())<=35:
                        item["rows"].append({"values":row_values(tab,r),"separation_arcsec":d})
            dump(dest/(name+".json"),item)
        tab,server=li.result()
    gates=src.id.map(vet.match_radius_arcsec)
    coords=SkyCoord(src.ra.to_numpy()*u.deg,src.dec.to_numpy()*u.deg)
    matched=[]
    detailed=[]
    if tab is not None:
        catalogue=SkyCoord([position(tab,r) for r in tab])
        a,b,sep,_=search_around_sky(coords,catalogue,35*u.arcsec)
        for i,j,s in zip(a,b,sep.arcsec):
            if s>gates.iloc[i]:
                continue
            det=src.iloc[i]
            r=tab[j]
            rec={"detection_uid":det.id,"source_group":det.source_group,"slew_name":det["name"],
                 "li_row":int(j+1),"li_ra":catalogue[j].ra.deg,"li_dec":catalogue[j].dec.deg,
                 "sep_arcsec":round(float(s),3),"gate_arcsec":gates.iloc[i]}
            matched.append(rec)
            if det.id in FOCUS:
                detailed.append(dict(rec,original_fields=row_values(tab,r)))
        summary["Li"]={"returned_rows":len(tab),"published_population":265,
            "searched_leads":len(src),"matching_detection_uids":len({r["detection_uid"] for r in matched}),
            "matching_slew_source_groups":len({r["source_group"] for r in matched}),
            "matching_distinct_Li_rows":len({r["li_row"] for r in matched}),
            "rows_sha256":hash_table(tab),"server":server,"columns":metadata(tab),"focal_rows":detailed}
    else:
        summary["Li"]={"status":server}
    pd.DataFrame(matched,columns=["detection_uid","source_group","slew_name","li_row","li_ra","li_dec",
                                 "sep_arcsec","gate_arcsec"]).to_csv(dest/"li_overlap.csv",index=False)
    with gzip.open(out/"vetting/vetting_provenance.jsonl.gz","rt") as f:
        events=[json.loads(line) for line in f]
    summary["TNS"]=[a for a in events if a.get("service") in {"TNS-broker-coverage","TNS-public-HTML-audit",
        "TNS-replica-audit","TNS","TNS-cone","TNS-broker-snapshot","TNS-mirror-inventory","TNS-mirror",
        "TNS-secondary","TNS-mirror-LFS"}]
    summary["Fink_snapshot_diagnostics"]=fink_snapshot(events)
    summary["TNS_search_outcomes"]=vet.service_outcomes.map(lambda v:json.loads(v).get("TNS")).value_counts().to_dict()
    dump(dest/"summary.json",summary)
    print("Original publication context",summary["Li"].get("returned_rows"),"Li sources;",
          summary["Li"].get("matching_detection_uids"),"slew detections in the prior sample",flush=True)


if __name__=="__main__":
    run(Path(sys.argv[1]))

"""Search public ANTARES alert loci with the NOIRLab client. Only alert-conditioned coverage."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from importlib import metadata
import hashlib
import inspect
import json
from pathlib import Path
import re
import subprocess
import sys

import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time

RADIUS = 35
CATALOG = "tns_public_objects"
CALIBRATORS = {"AT2018cow":(244.000925,22.268008),"SN2023ixf":(210.9107,54.31165),
               "AT2017gfo":(197.4504,-23.3815),"SN2016gkg":(23.5603,-29.4403)}


def safe(v):
    if hasattr(v,"to_dict"):
        v=v.to_dict()
    elif hasattr(v,"__dict__") and not isinstance(v,type):
        v={k:w for k,w in vars(v).items() if not k.startswith("_")}
    return json.loads(json.dumps(v,default=str))


def time_range(properties):
    times={}
    for name,value in properties.items():
        if "time" in name.lower() or "mjd" in name.lower() or "date" in name.lower():
            if value is None:
                continue
            if isinstance(value,(int,float)) and 30000<float(value)<90000:
                times[name]={"mjd":value,"utc":Time(value,format="mjd",scale="utc").isot}
            elif isinstance(value,(float,int,str)):
                times[name]=value
    return times


def one(field):
    from antares_client.search import cone_search
    sky=SkyCoord(field["ra"]*u.deg,field["dec"]*u.deg)
    records=[]
    cursor=cone_search(sky,Angle(RADIUS,unit=u.arcsec))
    for locus in cursor:
        ra,dec=float(locus.ra),float(locus.dec)
        offset=sky.separation(SkyCoord(ra*u.deg,dec*u.deg)).arcsec
        if offset>RADIUS+.2:
            raise ValueError("ANTARES returned a locus outside the requested angular cone")
        props=getattr(locus,"properties",{}) or {}
        cats=getattr(locus,"catalogs",[]) or []
        objects=getattr(locus,"catalog_objects",{}) or {}
        if not isinstance(objects,dict):
            raise ValueError("Unexpected ANTARES catalog_objects representation: "+type(objects).__name__)
        registered=sorted(set([str(x) for x in cats]+list(objects.keys())))
        tns=[(c,v) for c,v in objects.items() if c.lower()==CATALOG]
        tns_objects=[]
        for c,entries in tns:
            if not isinstance(entries,list):
                entries=[entries]
            tns_objects.extend({"catalog":c,"published_record":safe(entry)} for entry in entries)
        records.append({"locus_id":str(getattr(locus,"locus_id",getattr(locus,"id",""))),
                        "ra":ra,"dec":dec,"sep_arcsec":round(float(offset),3),
                        "ztf_object_id":props.get("ztf_object_id"),"catalog_names":registered,
                        "association_marked":CATALOG in registered,"tns_records":tns_objects,
                        "alert_counts":{k:v for k,v in props.items() if k in ("num_alerts","num_mag_values","num_filters")},
                        "broker_time_fields":time_range(props)})
    return {"status":"ok","loci":records}


def bounded(field):
    try:
        p=subprocess.run([sys.executable,__file__,"--field",json.dumps(field)],text=True,
                         capture_output=True,timeout=65)
        if p.returncode:
            raise RuntimeError((p.stdout+" "+p.stderr).strip()[-600:])
        return json.loads(p.stdout.splitlines()[-1])
    except Exception as exc:
        return {"status":"unavailable","reason":str(exc)[:700],"loci":[]}


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    bundle={"queried_utc":datetime.now(timezone.utc).isoformat(),"radius_arcsec":RADIUS,
        "service":"NOIRLab ANTARES","catalogue":CATALOG,
        "scope":"TNS associations of returned alert loci only. Absence is conditional on a recorded locus and ANTARES catalogue matching. ZTF sky and time coverage is incomplete.",
        "fields":[]}
    try:
        from antares_client.search import cone_search
        bundle["client"]={"package":"antares-client","version":metadata.version("antares-client"),
            "function":cone_search.__module__+"."+cone_search.__name__,"signature":str(inspect.signature(cone_search)),
            "function_source_sha256":hashlib.sha256(inspect.getsource(cone_search).encode()).hexdigest()}
    except Exception as exc:
        bundle["client"]={"error":str(exc)}
    try:
        url="https://api.antares.noirlab.edu/v1/catalogs"
        resp=requests.get(url,timeout=(6,20))
        resp.raise_for_status()
        catalogues=resp.json()
        bundle["catalogue_registration"]={"url":url,"sha256":hashlib.sha256(resp.content).hexdigest(),
            "rows":[x for x in catalogues.get("data",[]) if CATALOG in json.dumps(x).lower()]}
    except Exception as exc:
        bundle["catalogue_registration"]={"error":str(exc)[:450]}
    fields=[{"name":n,"role":"named TNS control","ra":ra,"dec":dec} for n,(ra,dec) in CALIBRATORS.items()]
    reg=pd.read_csv(out/"review_register.csv.gz",low_memory=False).fillna("")
    for _,r in reg.iterrows():
        forced=r.detection_uid in ("XMMSL2-total-69754","XMMSL2-total-10038")
        galaxy=bool(re.search(r"^(?:G|Galaxy|AGN|QSO|Sy[12]|LINER|BLLac|EXP|DEV|SER|REX)$",str(r.host_type),re.I))
        survives=str(r.disposition) in ("candidate","inconclusive") and pd.to_numeric(
            r.get("robust_model_floor_ratio",""),errors="coerce")>=10
        if forced or (galaxy and survives):
            fields.append({"name":r["name"],"detection_uid":r.detection_uid,"role":"galaxy field",
                           "ra":float(r.ra),"dec":float(r.dec)})
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f,reply in zip(fields,pool.map(bounded,fields)):
            record=dict(f,**reply)
            if f["role"]=="named TNS control":
                expected=re.sub(r"^(at|sn)","",f["name"].lower())
                matched=[]
                for l in reply.get("loci",[]):
                    payload=json.dumps(l["tns_records"]).lower().replace(" ","")
                    if expected in payload:
                        matched.append(l["locus_id"])
                record["named_control_recovered"]=bool(matched)
                record["matching_locus_ids"]=matched
            bundle["fields"].append(record)
    bundle["controls_recovered"]=[r["name"] for r in bundle["fields"] if r.get("named_control_recovered")]
    bundle["successful_galaxy_cones"]=sum(r["role"]=="galaxy field" and r["status"]=="ok" for r in bundle["fields"])
    (out/"broker_cones.json").write_text(json.dumps(bundle,indent=2,default=str)+"\n")
    print(f"ANTARES: {bundle['successful_galaxy_cones']} measured galaxy cones; "
          f"{len(bundle['controls_recovered'])}/{len(CALIBRATORS)} benchmark names recovered")


if __name__=="__main__":
    if sys.argv[1]=="--field":
        try:
            print(json.dumps(one(json.loads(sys.argv[2])),default=str))
        except Exception as exc:
            print(json.dumps({"status":"unavailable","reason":type(exc).__name__+": "+str(exc)[:650],"loci":[]}))
    else:
        run(Path(sys.argv[1]))

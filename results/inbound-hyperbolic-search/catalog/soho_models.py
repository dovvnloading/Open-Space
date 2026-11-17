"""Test a bound-energy boundary model against the spacecraft astrometry of two Kreutz objects."""
import argparse
from datetime import date
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import Request,urlopen

import numpy as np
from scipy.optimize import least_squares
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

HERE=Path(__file__).parent
OUT=HERE/"soho"
OUT.mkdir(exist_ok=True)
UA="Open-Space-SOHO/1.0 (+https://github.com/dovvnloading/Open-Space)"
AU=149597870.7
MU=.0002959122082855911
C=173.144632674240
receipts=[]
cat=json.loads(gzip.decompress((HERE/"sbdb.json.gz").read_bytes()))
rows=[dict(zip(cat["fields"],r)) for r in cat["data"]]

def fetch(url,label,body=None):
    r={"url":url,"file":label}
    if body:r.update(query=body,method="GET")
    receipts.append(r)
    try:
        with urlopen(Request(url,data=json.dumps(body).encode() if body else None,method="GET",
                             headers={"User-Agent":UA,"Content-Type":"application/json"}),timeout=90) as response:
            raw=response.read()
            r.update(status=response.status,date=response.headers.get("Date"),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        return raw
    except Exception as e:
        r["error"]=str(e)
        raise

def observations(identifier,packed,tag):
    rawfile=OUT/(tag+".mpc")
    if rawfile.exists():return rawfile.read_text().splitlines()
    for d in (identifier,identifier[2:],packed):
        try:
            raw=fetch("https://data.minorplanetcenter.net/api/get-obs",tag+"-obs",{"desigs":[d],"output_format":["OBS80"]})
            obs=json.loads(raw)[0]["OBS80"].splitlines()
            records=[l[:80] for l in obs if len(l)>=80 and l[5:12]==packed and l[77:80]=="249"]
            if records:break
        except Exception:pass
    else:
        raw=fetch("https://www.minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz",tag+"-bulk")
        records=[l[:80].decode() for l in gzip.decompress(raw).splitlines() if l[5:12].decode()==packed and l[77:80]==b"249"]
    if not records:raise ValueError("No matching MPC observations")
    rawfile.write_text("\n".join(records)+"\n")
    return records

def parse(lines):
    observations=[]
    sat={}
    for l in lines:
        yr,mo,dd=l[15:32].split()
        jd=date(int(yr),int(mo),1).toordinal()+1721424.5+float(dd)-1
        if l[14]=="s":
            # The geocentric observer coordinates are reported in the equatorial J2000 frame.
            # Satellite record columns 33 and 35:69 contain a unit code and three signed values.
            coords=re.findall(r"[+-]\s*\d+(?:\.\d+)?",l[34:70])
            if len(coords)!=3:raise ValueError(f"Bad spacecraft record {l!r}")
            scale=1/AU if l[32]=="1" else 1 if l[32]=="2" else None
            if scale is None:raise ValueError("Unknown satellite units "+l[32])
            sat[jd]=np.array([float(x.replace(" ",""))*scale for x in coords])
            continue
        if l[14]!="S":raise ValueError(f"Unexpected observer type {l!r}")
        h,m,s=map(float,l[32:44].split())
        d,dm,ds=map(float,l[44:56].split())
        ra=15*(h+m/60+s/3600)
        dec=(-1 if l[44:56].strip().startswith("-") else 1)*(abs(d)+dm/60+ds/3600)
        observations.append((jd,ra,dec))
    observations.sort()
    j,ra,dec=np.asarray(observations).T
    if any(t not in sat for t in j):raise ValueError("Missing SOHO observer vector")
    off=np.array([sat[t] for t in j])
    if not np.all((np.linalg.norm(off,axis=1)>.005)&(np.linalg.norm(off,axis=1)<.02)):
        raise ValueError("SOHO offsets not compatible with heliocentric spacecraft at L1")
    t=Time(j,format="jd",scale="utc")
    earth=np.asarray((get_body_barycentric("earth",t)-get_body_barycentric("sun",t)).xyz.to_value("AU")).T
    return t.tdb.jd,np.radians(ra),np.radians(dec),earth+off,off

def position(t,el):
    q,tp,inc,node,arg=el
    inc,node,arg=np.radians([inc,node,arg])
    B=(t-tp)*np.sqrt(MU/(2*q**3))
    d=2*np.sinh(np.arcsinh(1.5*B)/3)
    xp,yp=q*(1-d*d),2*q*d
    x1,y1=np.cos(arg)*xp-np.sin(arg)*yp,np.sin(arg)*xp+np.cos(arg)*yp
    x=np.cos(node)*x1-np.sin(node)*np.cos(inc)*y1
    y=np.sin(node)*x1+np.cos(node)*np.cos(inc)*y1
    z=np.sin(inc)*y1
    eps=np.radians(23.439291111)
    return np.stack([x,np.cos(eps)*y-np.sin(eps)*z,np.sin(eps)*y+np.cos(eps)*z],axis=1)

def analyze(target):
    row=next(r for r in rows if r["full_name"].strip().startswith(target+" "))
    pdes=row["pdes"]
    yr,des=pdes.split()
    packed="K"+yr[2:]+des[0]+f"{int(des[1:]):02d}0"
    tag=packed
    lines=observations(target,packed,tag)
    t,ra,dec,obs,off=parse(lines)
    seed=np.array([float(row[k]) for k in ("q","tp","i","om","w")])
    def res(x):
        el=seed+x*np.array([.001,.01,1,1,1])
        emitted=t
        for _ in range(5):
            dr=position(emitted,el)-obs
            emitted=t-np.linalg.norm(dr,axis=1)/C
        dr=position(emitted,el)-obs
        a=np.arctan2(dr[:,1],dr[:,0]);de=np.arctan2(dr[:,2],np.hypot(dr[:,0],dr[:,1]))
        return np.column_stack([np.arctan2(np.sin(ra-a),np.cos(ra-a))*np.cos(dec),dec-de]).ravel()*206264.806247
    fit=least_squares(res,np.zeros(5),bounds=([(.0002-seed[0])/.001,-150,-90,-180,-180],[(.08-seed[0])/.001,150,90,180,180]),
                      x_scale="jac",max_nfev=1500,ftol=1e-10,xtol=1e-10,gtol=1e-10,diff_step=1e-5)
    def metrics(v):
        r=v.reshape(-1,2)
        return {"coordinate_rms_arcsec":float(np.sqrt(np.mean(r*r))),"two_coordinate_radial_rms_arcsec":float(np.sqrt(np.mean(np.sum(r*r,axis=1)))),
                "median_radial_arcsec":float(np.median(np.linalg.norm(r,axis=1)))}
    fitted=seed+fit.x*np.array([.001,.01,1,1,1])
    return {"target":target,"packed":packed,"positions":len(t),"continuations":sum(l[14]=="s" for l in lines),
            "site":"249 SOHO; heliocentric observer = Earth ERFA epv00 + published satellite offset",
            "observer_geocentric_range_au":[float(np.linalg.norm(off,axis=1).min()),float(np.linalg.norm(off,axis=1).max())],
            "arc_days":float(t.max()-t.min()),"span_TDB":[float(t.min()),float(t.max())],"eccentricity_fixed":1.0,
            "parameters":["q_au","tp_jd_tdb","i_deg","ascending_node_deg","argument_perihelion_deg"],
            "published":seed.tolist(),"published_residuals":metrics(res(np.zeros(5))),
            "optimized":fitted.tolist(),"optimized_residuals":metrics(res(fit.x)),"converged":bool(fit.success),
            "fit_message":fit.message,"all_east_north_residual_arcsec":res(fit.x).reshape(-1,2).tolist(),
            "obs_sha256":hashlib.sha256((OUT/(tag+".mpc")).read_bytes()).hexdigest()}
p=argparse.ArgumentParser()
p.add_argument("targets",nargs="*",default=["C/2005 D3","C/2006 S8"])
answers=[]
for target in p.parse_args().targets:
    try: result=analyze(target)
    except Exception as e: result={"target":target,"error":str(e)}
    answers.append(result)
    print("SOHO",target,{k:v for k,v in result.items() if "all_east" not in k},flush=True)
(OUT/"results.json").write_text(json.dumps({"targets":answers,"receipts":receipts},indent=2,allow_nan=False)+"\n")

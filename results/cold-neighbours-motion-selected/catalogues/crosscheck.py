"""Motion-aware literature association screen. A positional alert requires inspection."""
import csv
import hashlib
import io
import json
import math
import re
import sys
import zipfile
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
def acquire():
    path=ROOT/"inputs"
    path.mkdir(exist_ok=True)
    manifest=json.loads((ROOT/"inputs.json").read_text())
    for key,record in manifest.items():
        dest=path/key
        if not dest.exists():
            with urlopen(record["url"], timeout=180) as f:b=f.read()
            if b.startswith(b"#"):
                stable=b"\n".join(l for l in b.splitlines() if l.strip() and not l.startswith(b"#"))
            else:stable=b
            if hashlib.sha256(stable).hexdigest()!=record["payload_sha256"]:
                raise ValueError(f"Data changed since this analysis: {key}")
            dest.write_bytes(b)
    return path

def tsv(path,table):
    parts=path.read_text().split("#Name: ")
    s=next(t for t in parts if t.split("\n",1)[0].endswith("/"+table))
    lines=[l for l in s.splitlines()[1:] if l.strip() and not l.startswith("#")]
    return csv.DictReader([lines[0],*lines[3:]],delimiter="\t")

def sex(s):
    a,b,c=s.replace("+","").lstrip("-").split(":")
    return (-1 if s.startswith("-") else 1)*(float(a)+float(b)/60+float(c)/3600)
def uvec(ra,dec):
    a,d=np.radians([ra,dec])
    return np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)]).T
def catalogue():
    path=acquire(); refs=[]
    def entry(key,name,ra,dec,epoch=None):
        try:a,d=float(ra),float(dec)
        except (ValueError,TypeError):return
        if np.isfinite(a+d): refs.append((key,name.strip(),a,d,epoch))
    def rows(key):return csv.DictReader((path/key).open(encoding="utf-8-sig"))
    for r in rows("UCS"):entry("UCS",r["name"],r["ra_j2000_formula"],r["dec_j2000_formula"],2000.0)
    for r in rows("20pc"):entry("20pc",r["DefaultName"],r["RA (deg)"],r["Dec (deg)"],float(r["Epoch"]) if r["Epoch"].strip() else None)
    for l in (path/"byw26").read_text().splitlines():
        m=re.match(r'^J(\d\d)(\d\d)(\d\d\.\d+)([+-]\d\d)(\d\d)(\d\d\.\d+)',l)
        if m:
            h,mi,s,d,dm,ds=m.groups()
            entry("byw26",m[0],15*sex(":".join((h,mi,s))),sex(":".join((d,dm,ds))))
    with zipfile.ZipFile(path/"spiff") as z:
        k=next(n for n in z.namelist() if n.endswith("/spherex_autotype_datasets_fulltable.csv") and not n.startswith("__"))
        for r in csv.DictReader(io.TextIOWrapper(z.open(k))):
            entry("spiff",r["designation"],15*sex(r["ra_j2000"]),sex(r["dec_j2000"]))
    specs=[
        ("karpc","table2","Candid","RA_ICRS","DE_ICRS",2016.0),
        ("euclid_spec","catalog","EUCLID","RA_ICRS","DE_ICRS",2024.25),
        ("euclid_phot","tablea1","ObjectId","RAJ2000","DEJ2000",2024.25),
        ("euclid_phot","tablea3","ObjectId","RAJ2000","DEJ2000",2024.25),
        ("byw_comp","table1","ID","_RA","_DE",None),
        ("byw20","table6","Source","RA_ICRS","DE_ICRS","MJD"),
        ("catwise19","bdcand","Name","RA_ICRS","DE_ICRS","MJD"),
    ]
    for k,t,idcol,ra,dec,ep in specs:
        for r in tsv(path/k,t):
            epoch=(2000+(float(r["MJD"])-51544.5)/365.25 if r["MJD"].strip() else None) if ep=="MJD" else ep
            entry(k,r[idcol],r[ra],r[dec],epoch)
    return refs

def run(input_path,output_path):
    refs=catalogue()
    pts=cKDTree(uvec([r[2] for r in refs],[r[3] for r in refs]))
    found=[]
    for r in csv.DictReader(Path(input_path).open()):
        a,d,mx,my=map(float,[r[k] for k in ("ra2015_deg","dec2015_deg","pmra","pmdec")])
        mu=math.hypot(mx,my)
        # Covers 1985–2027 including FAST movers. Wrong-epoch matches are retained as alerts.
        radius=(31*mu+5)*math.pi/(180*3600)
        for idx in pts.query_ball_point(uvec(a,d),radius):
            key,name,ra,dec,epoch=refs[idx]
            dx=((ra-a+180)%360-180)*math.cos(math.radians(d))*3600;dy=(dec-d)*3600
            t=(np.clip((dx*mx+dy*my)/(mu*mu),-30.405,11.595)+2015.405) if epoch is None else epoch
            sep=math.hypot(dx-(t-2015.405)*mx,dy-(t-2015.405)*my)
            if sep<=4.0:
                found.append(dict(source_name=r["source_name"],catalogue=key,matched_label=name,separation_arcsec=round(sep,3),epoch=round(t,4),epoch_mode="closest_approach" if epoch is None else "catalogue_or_release_epoch"))
    with Path(output_path).open("w") as f:
        w=csv.DictWriter(f,fieldnames=["source_name","catalogue","matched_label","separation_arcsec","epoch","epoch_mode"]);w.writeheader();w.writerows(found)
    from collections import Counter
    print(dict(Counter(r[0] for r in refs)),len(found),len(set(r["source_name"] for r in found)))
if __name__=="__main__":run(*sys.argv[1:])

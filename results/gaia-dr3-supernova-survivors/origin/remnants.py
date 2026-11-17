"""Recent angular traceback against Green's October 2024 catalogue, with perspective."""
import csv
import hashlib
import json
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
URL = "https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html"
SHA = "0dcffa17be68c3f2467d5e80eae40c6dddba68195b5a4765556c6b4e9f966396"
KMS_PC_PER_YEAR = 1.022712165e-6
MAS_RAD = np.pi / 648000000


def axes(ra, dec):
    a,b=np.deg2rad([ra,dec])
    return (np.array([np.cos(a)*np.cos(b),np.sin(a)*np.cos(b),np.sin(b)]),
            np.array([-np.sin(a),np.cos(a),0.0]),
            np.array([-np.cos(a)*np.sin(b),-np.sin(a)*np.sin(b),np.cos(b)]))


def main():
    data=urlopen(URL,timeout=50).read()
    if hashlib.sha256(data).hexdigest()!=SHA:
        raise ValueError("Green's catalogue content changed")
    entries=[]
    for line in data.decode().splitlines():
        m=re.search(r'HREF="snrs\.(G[0-9.+-]+)\.html".*?</A>\s+([\d ]+?)\s+([-+]\d\d) (\d\d)\s+(\S+)',line)
        if not m:
            continue
        name,hms,deg,minute,size=m.groups()
        h,mi,sec=map(float,hms.split())
        ra=15*(h+mi/60+sec/3600)
        dec=(-1 if deg.startswith("-") else 1)*(abs(int(deg))+float(minute)/60)
        lengths=re.findall(r'\d+(?:\.\d+)?',size)
        if lengths:
            entries.append((name,axes(ra,dec)[0],max(map(float,lengths))/120))
    if len(entries)<300:
        raise ValueError(f"incomplete parse ({len(entries)})")
    def keyed(p):
        with p.open() as f:return {r["source_id"]:r for r in csv.DictReader(f)}
    bj=keyed(HERE.parent/"control"/"bj.csv")
    kin=keyed(HERE.parent/"control"/"published_kinematics.csv")
    stars=keyed(HERE.parent/"control"/"controls.csv")
    g={sid:row for path in sorted(HERE.glob("rows_*.csv")) for sid,row in keyed(path).items()}
    t=np.arange(0.,200001.,250.)
    unit=np.array([e[1] for e in entries])
    radii=np.array([e[2] for e in entries]);out=[]
    for sid,row in g.items():
        n,e,north=axes(float(row["ra"]),float(row["dec"]))
        d=8500.0 if stars[sid]["name"]=="US 708" else float(bj[sid]["r_med_geo"])
        v=float(kin[sid]["rv_kms"])*KMS_PC_PER_YEAR*n + d*MAS_RAD*(float(row["pmra"])*e+float(row["pmdec"])*north)
        track=d*n-t[:,None]*v
        track/=np.linalg.norm(track,axis=1)[:,None]
        angle=np.rad2deg(np.arccos(np.clip(track@unit.T,-1,1)))
        tmin=angle.argmin(axis=0)
        sep=angle[tmin,np.arange(len(entries))]
        for j in range(len(entries)):
            if sep[j]<=radii[j] or j==np.argmin(sep/radii):
                out.append(dict(source_id=sid,name=stars[sid]["name"],remnant=entries[j][0],
                    hit=bool(sep[j]<=radii[j]),time_yr=t[tmin[j]],separation_deg=round(sep[j],5),
                    max_radius_deg=round(radii[j],5)))
    with (HERE/"remnant_checks.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    (HERE/"remnant_receipt.json").write_text(json.dumps(dict(url=URL,sha256=SHA,bytes=len(data),
        entries_in_source=data.count(b'HREF="snrs.G'),sources_with_parsed_size=len(entries),
        stars=len(g),dt_yr=250,lookback_yr=200000),indent=2)+"\n")


if __name__=="__main__":
    main()

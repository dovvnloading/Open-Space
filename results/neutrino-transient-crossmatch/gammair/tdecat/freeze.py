"""Freeze W2 visit measurements for the entire published Langis IR-flare list."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import statistics as st
from urllib.parse import quote
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
REV="bca013a266cd5d2a120aa4c6b34a2dc54bee1e5e"
BASE=f"https://raw.githubusercontent.com/dlangis/TDECat/{REV}/"
CAT="Section 1: TDECat/TDE_catalogue_all.csv"
FOLDER="Section 2: Photometry/Optical and Infrared/"
NAMES="AT2023ugy AT2023cvb AT2022upj AT2022fpx AT2022agi AT2022dyt AT2022aee AT2021uqv AT2020afhd AT2020ksf AT2020mot AT2020nov AT2020pj AT2019qiz AT2019dsg AT2019azh AT2017gge ASASSN-14li".split()
BASE_MJD=dt.datetime(1858,11,17)
def sexagesimal(s):
    a,b,c=map(float,s.lstrip("+-").split(":"))
    return (-1 if s.startswith("-") else 1)*(a+b/60+c/3600)
def write_csv(file,rows):
    with (HERE/file).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
def visit_table(data):
    obs=[]
    for r in csv.DictReader(io.StringIO(data.decode()),delimiter=";"):
        if r["Facility"] not in {"NEOWISE","ALLWISE"} or r["Filter"]!="WISE(W2)":
            continue
        try:
            t,m,e=map(float,[r["MJD"],r["Magnitude"],r["Error"]])
        except ValueError:
            continue
        if 0<e<=1.085736/3 and all(math.isfinite(x) for x in (t,m,e)):
            obs.append((t,m,e))
    groups=[]
    for row in sorted(obs):
        if not groups or row[0]-groups[-1][-1][0]>30:
            groups.append([])
        groups[-1].append(row)
    return [dict(mjd=st.median(r[0] for r in g),w2_mag=st.median(r[1] for r in g),
                n=len(g),mag_scatter_mad=st.median(abs(r[1]-st.median(x[1] for x in g)) for r in g))
            for g in groups if len(g)>=3]
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cache",type=Path,help="persistent raw files downloaded from this pinned revision")
    args=p.parse_args()
    receipts=[];visits=[];peaks=[]
    def read(path,local):
        url=BASE+quote(path)
        data=(args.cache/local).read_bytes() if args.cache else urlopen(url,timeout=60).read()
        receipts.append(dict(url=url,revision=REV,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
        return data
    sources=list(csv.DictReader(io.StringIO(read(CAT,"tdecat.csv").decode())))
    for name in NAMES:
        r=next(v for v in sources if name in [v["AT name"].replace(" ",""),v["Alternative name"]])
        disc=(dt.datetime.fromisoformat(r["Discovery date (UT)"])-BASE_MJD).total_seconds()/86400
        source_visits=visit_table(read(FOLDER+f"target_{name}_photometry.csv",f"tde/{name}.csv"))
        visits.extend(dict(source_id=name,**v) for v in source_visits)
        after=[v for v in source_visits if v["mjd"]>=disc]
        peak=min(after,key=lambda v:v["w2_mag"])
        k=source_visits.index(peak)
        peaks.append(dict(source_id=name,ra_deg=round(15*sexagesimal(r["RA"]),9),dec_deg=round(sexagesimal(r["DEC"]),9),
                peak_mjd=peak["mjd"],peak_utc=(BASE_MJD+dt.timedelta(days=peak["mjd"])).strftime("%Y-%m-%dT%H:%M:%SZ"),
                w2_mag=peak["w2_mag"],peak_n=peak["n"],previous_mjd=source_visits[k-1]["mjd"] if k else "",
                next_mjd=source_visits[k+1]["mjd"] if k+1<len(source_visits) else "",
                right_censored=int(k+1==len(source_visits)),anchor_kind="brightest observed post-discovery W2 visit"))
    write_csv("w2_visits.csv",visits)
    write_csv("observed_maxima.csv",peaks)
    (HERE/"receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")
if __name__=="__main__":
    main()

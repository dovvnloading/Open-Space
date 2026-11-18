"""Build both MIR tables without reference to alerts; from repo root give path to checkout."""
import argparse
import csv
from collections import defaultdict
import datetime as dt
import gzip
import hashlib
import io
import json
from math import asin,cos,radians,sin,sqrt
from pathlib import Path
import subprocess
import urllib.request

HERE=Path(__file__).resolve().parent
URL="https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/695/A228&-out.max=unlimited&-out.all"
REV="7d266947bd5c64e945ab040173b5ffb660c38ca3"
BASE="results/neowise-nuclear-flares/data/"
FIELDS=["source_id","ra_deg","dec_deg","peak_mjd","peak_utc","population","source_class","published_alias","flaires_host_within_5arcsec"]

def instant(mjd):
    return (dt.datetime(1858,11,17)+dt.timedelta(days=float(mjd))).isoformat(timespec="seconds")+"Z"

def dist_arcsec(a,b,c,d):
    x=sin(radians(d-b)/2)**2+cos(radians(b))*cos(radians(d))*sin(radians(c-a)/2)**2
    return 7200*asin(sqrt(min(x,1)))*180/3.141592653589793

def write(name,rows):
    with (HERE/name).open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("checkout",type=Path)
    p.add_argument("--download",action="store_true")
    opts=p.parse_args()
    pth=HERE/"flaires_raw.tsv.gz"
    if opts.download:
        with urllib.request.urlopen(URL,timeout=60) as r:
            pth.write_bytes(gzip.compress(r.read(),mtime=0))
    raw=gzip.decompress(pth.read_bytes())
    lines=[s for s in raw.decode().splitlines() if s.strip() and not s.startswith("#")]
    lines=lines[:1]+lines[3:]  # omit ASU unit and ruler rows
    main_rows=[]
    for row in csv.DictReader(lines,delimiter="\t"):
        name=row["Name"].strip()
        mjd=row["RefTime"].strip()
        main_rows.append(dict(source_id=name,ra_deg=float(row["RAJ2000"]),dec_deg=float(row["DEJ2000"]),
                 peak_mjd=mjd,peak_utc=instant(mjd),population="Flaires",source_class="dust-echo-like IR flare",
                 published_alias=row["TNSname"].strip(),flaires_host_within_5arcsec=""))
    assert len(main_rows)==823
    write("flaires_peaks.csv",main_rows)
    rec=[dict(dataset="Flaires",url=URL,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
              response_date=next(s.split("=")[1].split("\t")[0] for s in raw.decode().splitlines() if "request_date=" in s))]
    def git_csv(name):
        obj=BASE+name
        b=subprocess.check_output(["git","-C",str(opts.checkout),"show",f"{REV}:{obj}"])
        rec.append(dict(dataset=name,url=f"https://github.com/dovvnloading/Open-Space/blob/{REV}/{obj}",bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
        return list(csv.DictReader(io.StringIO(b.decode())))
    visits=defaultdict(list)
    for row in git_csv("flare-light-curves.csv"):
        visits[row["id"]].append(float(row["mjd"]))
    candidates=[]
    for row in git_csv("flares.csv"):
        if row["class"]=="off-nuclear" or row["vetted"]!="True":
            continue
        # Identify the unique visit within one day of the reported rounded peak.
        approx=51544.5+(float(row["peak_year"])-2000)*365.25
        dates=[v for v in visits[row["id"]] if abs(v-approx)<=1.0]
        if len(dates)!=1:
            raise ValueError(f"ambiguous visit for {row['id']}")
        ra,dec=float(row["ra"]),float(row["dec"])
        overlap=[r["source_id"] for r in main_rows if dist_arcsec(ra,dec,r["ra_deg"],r["dec_deg"])<=5]
        candidates.append(dict(source_id=row["id"],ra_deg=ra,dec_deg=dec,peak_mjd=dates[0],
                    peak_utc=instant(dates[0]),population="Open-Space 2MRS",source_class=row["class"],
                    published_alias=row["published_as"],flaires_host_within_5arcsec="|".join(overlap)))
    write("supplement_2mrs_peaks.csv",candidates)
    (HERE/"mir_receipts.json").write_text(json.dumps(rec,indent=2)+"\n")

if __name__=="__main__":
    main()

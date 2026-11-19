#!/usr/bin/env python3
"""Independent NASA Swift GRB-table archive read: fixed +/-1d / +/-1h windows.

This is a Swift GRB table (not the all-mission GCN/GBM population). Source
positions outside this table and spacecraft non-coverage remain untested.
"""
from datetime import datetime, timedelta, timezone
from html.parser import HTMLParser
from hashlib import sha256
from pathlib import Path
from urllib.request import Request,urlopen
from math import cos,sin,acos,radians,degrees
import csv,json,random,re

ROOT=Path(__file__).resolve().parent
URL="https://swift.gsfc.nasa.gov/archive/grb_table/fullview/"
T0=datetime(1858,11,17,tzinfo=timezone.utc)+timedelta(days=59988.0533299)
ra0,dec0=94.3,-7.8

class T(HTMLParser):
    def __init__(self):
        super().__init__();self.in_row=False;self.in_cell=False;self.cells=[];self.rows=[];self.cell=""
    def handle_starttag(self,tag,attrs):
        if tag=="tr":self.in_row=True;self.cells=[]
        if tag in ("td","th") and self.in_row:self.in_cell=True;self.cell=""
    def handle_data(self,d):
        if self.in_cell:self.cell+=d
    def handle_endtag(self,tag):
        if tag in ("td","th") and self.in_cell:
            self.cells.append(self.cell.strip());self.in_cell=False
        if tag=="tr" and self.in_row:
            if self.cells:self.rows.append(self.cells)
            self.in_row=False
def sep(ra,dec):
    s=sin((radians(dec-dec0))/2)**2+cos(radians(dec))*cos(radians(dec0))*sin(radians(ra-ra0)/2)**2
    return degrees(2*__import__("math").asin(min(1,s**.5)))
def f(x):
    try:return float(re.match(r"\s*([+-]?\d+(?:\.\d+)?)",str(x)).group(1))
    except (AttributeError,ValueError):return None

record={"source":URL,"collected_utc":datetime.now(timezone.utc).isoformat(),
        "time_window":"2023-02-12 01:16:47.7 to 2023-02-14 01:16:47.7 UTC",
        "caveat":"Swift archive includes selected Swift-localized GRBs, NOT all GCN/Fermi events"}
matches=[]
try:
    with urlopen(Request(URL,headers={"User-Agent":"Open-Space-KM3-source-catalog-audit/1.0"}),timeout=70) as u:
        buf=u.read(40_000_001)
        if len(buf)>40_000_000:raise ValueError("HTML exceeds 40MB safety cap")
        record.update({"http_status":u.status,"bytes":len(buf),"sha256":sha256(buf).hexdigest()})
    html=buf.decode("utf-8",errors="replace")
    t=T();t.feed(html)
    burst=[];seen=set()
    for row in t.rows:
        if len(row)<5:continue
        m=re.match(r"^(?:GRB\s*)?(\d{6}[A-Z])",row[0].strip())
        if not m:continue
        name=m.group(1)
        if name in seen:continue
        seen.add(name)
        year,month,day=2000+int(name[:2]),int(name[2:4]),int(name[4:6])
        try:
            day=datetime(year,month,day,tzinfo=timezone.utc)
        except ValueError:continue
        timestr=row[1].strip()
        try:
            hh,mm,ss=[float(v) for v in timestr.split(":")[:3]]
            when=day+timedelta(hours=hh,minutes=mm,seconds=ss)
        except (ValueError,IndexError):continue
        dt=(when-T0).total_seconds()/86400
        ra,dec=f(row[3]),f(row[4])
        if ra is None or dec is None:continue
        burst.append({"name":"GRB "+name,"ra":ra,"dec":dec,"dt":dt,"date_utc":when.isoformat()})
    if len(burst)<1000:raise ValueError("Unexpectedly few parseable Swift GRBs: "+str(len(burst)))
    filt=[x for x in burst if abs(x["dt"])<=1]
    matches=[]
    for v in filt:
        d=sep(v["ra"],v["dec"])
        if d<=3:
            matches.append({"name":v["name"],"ra_deg":v["ra"],"dec_deg":v["dec"],
                "delta_days":round(v["dt"],8),"angular_sep_deg":round(d,6),
                "inside_90pct":d<=2.2,"within_1h":abs(v["dt"])<=1/24,
                "status":"UNVERIFIED_SWIFT_GRB_COINCIDENCE","physical_chance_p":"NOT_ESTIMABLE_WITHOUT_ALL_MISSION_EXPOSURE"})
    shifts=[random.Random(230213+i).uniform(6.5,353.5) for i in range(999)]
    nullcounts=[]
    for s in shifts:
        ra_ref=(ra0+s)%360
        z=sum(degrees(2*__import__("math").asin(min(1,(sin(radians((v["dec"]-dec0)/2)**2+
           cos(radians(v["dec"]))*cos(radians(dec0))*sin(radians((v["ra"]-ra_ref)/2)**2)**.5)))<=3 for v in filt)
        nullcounts.append(z)
    obs=len(matches)
    record.update({"status":"COMPLETE_SWIFT_TABLE_ONLY","parsed_swift_grbs":len(burst),
          "within_1day":len(filt),"within_3deg_and_1day":len(matches),
          "within_1h_and_3deg":sum(r["within_1h"] for r in matches),
          "n_ra_shifts":len(nullcounts),"p_shift_count_ge_observed":round((1+sum(c>=obs for c in nullcounts))/1000,6),
          "chance_scope":"conditional only, not all-satellite GCN or time-dependent exposure",
          "uniform_all_sky_expected_3deg":len(filt)*28.267/41252.96125})
except Exception as ex:
    record.update({"status":"UNAVAILABLE_NOT_NEGATIVE","error":str(ex)[:300]})
with (ROOT/"swift_grb_matches.csv").open("w",newline="") as f:
    fields=["name","ra_deg","dec_deg","delta_days","angular_sep_deg","inside_90pct","within_1h","status","physical_chance_p"]
    w=csv.DictWriter(f,fields);w.writeheader();w.writerows(matches)
(ROOT/"swift_grb_receipt.json").write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
print(json.dumps({k:record[k] for k in ("status", "parsed_swift_grbs","within_1day","within_3deg_and_1day") if k in record}))

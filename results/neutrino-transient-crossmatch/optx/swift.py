"""Normalize the Swift discovery-observation and pointing projections."""
import csv
import gzip
import json
import re
import sys
from math import asin,sin,cos,radians,degrees,sqrt
from pathlib import Path

DIR=Path(__file__).resolve().parent
# Suppress externally triggered transient fields and explicit neutrino/GW triggers.
TRIGGER=re.compile(r"ICECUBE|NEUTRINO|AMON|^TOO\b|^HOST-SN|^(?:ZTF|WTP|TCP|LSXPS)|^IC[-_ ]?\d{6}|^S\d{6}[A-Za-z]|^(?:GRB|EP|AT|SN|ZTF|TCP|WTP)[-_ ]?(?:19|20|\d{6})|^(?:19|20)\d{2}[a-z]+$",re.I)


def distance(ra,dec,target_ra,target_dec):
    ra,dec,tra,tdec=map(radians,map(float,(ra,dec,target_ra,target_dec)))
    return degrees(2*asin(min(1,sqrt(sin((dec-tdec)/2)**2+cos(dec)*cos(tdec)*sin((ra-tra)/2)**2))))


def project(cache):
    from optical import write
    ids={}
    for tid,epoch,obs,sha,url,detail in json.loads((cache/"swift_details.json").read_text()):
        match=re.search(r"Obs\s+(\d{11})",detail)
        ids[tid]=(epoch,match[1] if match else obs,sha,url)
    with gzip.open(cache/"swiftmastr.gz","rt") as f:
        for line in f:
            if line.startswith("line["):
                keys=line.split("=")[1].split()
            if line.startswith("<DATA>"):break
        wanted={x[1] for x in ids.values()}
        pointings={}
        ix=keys.index("obsid")
        for line in f:
            parts=line.rstrip().split("|")
            if len(parts)>ix and parts[ix] in wanted:
                pointings[parts[ix]]={k:parts[keys.index(k)] for k in ["name","ra","dec","obsid"]}
    with gzip.open(cache/"swift11.gz","rt") as f:
        rows=[]
        for r in csv.DictReader(f):
            if r["TransientID"] not in ids:continue
            epoch,obs,sha,url=ids[r["TransientID"]]
            p=pointings.get(obs,{})
            rows.append(dict(source_id=r["IAUName"],transient_id=r["TransientID"],ra_deg=r["RA"],dec_deg=r["Decl"],
                             utc=epoch.replace(" ","T")+"Z",classification=r["Classification"],obsid=obs,target=p.get("name",""),
                             pointing_ra=p.get("ra",""),pointing_dec=p.get("dec",""),detail_sha256=sha,detail_url=url))
        write("data/swift_join.csv",rows)


def main():
    from optical import read, write
    if len(sys.argv)>1:project(Path(sys.argv[1]))
    rows=[];audit=[]
    for r in read("swift_join.csv"):
        reason=""
        if not r["target"]:reason="unverified discovery pointing"
        elif TRIGGER.search(r["target"]):reason="triggered transient field"
        elif distance(r["ra_deg"],r["dec_deg"],r["pointing_ra"],r["pointing_dec"])<2/60:
            reason="within 2 arcmin of pointing center"
        audit.append(dict(transient_id=r["transient_id"],source_id=r["source_id"],obsid=r["obsid"],target=r["target"],status=reason or "included"))
        if not reason:
            rows.append({**{k:r[k] for k in ["source_id","ra_deg","dec_deg","utc"]},"family":"xray",
                         "anchor":"first_XRT_detection","catalog":"Swift LSXPS","transient_id":r["transient_id"]})
    write("swift_events.csv",rows);write("swift_inventory.csv",audit)

if __name__=="__main__":
    main()

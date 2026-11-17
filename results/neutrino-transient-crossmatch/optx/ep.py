"""Extract autonomous EP-WXT coordinates and UTC from official GCN detection circulars.
Usage: python ep.py /path/to/archive.json.tar.gz
"""
import csv
import json
import re
import sys
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

DIR=Path(__file__).resolve().parent
CUTOFF=datetime.fromisoformat("2026-11-17T12:09:02+00:00").timestamp()*1000
EP=re.compile(r"\bEP\s*(\d{6})([a-z])\b",re.I)
ISSUE=re.compile(r"retract|correc|withdraw|spurious",re.I)
TITLE=re.compile(r"Einstein Probe.{0,32}detect|EP[- ]WXT.{0,32}detect|WXT.{0,20}detect",re.I)
WHEN=re.compile(r"(20\d\d[-/]\d\d[-/]\d\d)[ T]+(\d\d:\d\d:\d\d(?:\.\d+)?)",re.I)
POS=re.compile(r"\bR\.?\s*A\.?\s*[:=]\s*([\d.]+)\s*(?:deg(?:rees?)?\.?|°)?\s*[,;]?\s*(?:and\s+)?D(?:ec|\.E\.C)\.?\s*[:=]\s*([+−-]?[\d.]+)",re.I)


def dump(name,fields,rows):
    with (DIR/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fields,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    alerts=defaultdict(list)
    flags=defaultdict(list)
    with tarfile.open(sys.argv[1],"r:gz") as tar:
        for item in tar:
            if not item.name.endswith(".json"):
                continue
            rec=json.load(tar.extractfile(item))
            if not (1704067200000 <= rec.get("createdOn",0) <= CUTOFF):
                continue
            title=rec["subject"].strip()
            ids=EP.findall(title)
            if ISSUE.search(title):
                for date,letter in ids:
                    flags[date+letter.lower()].append(rec["circularId"])
            if "ep_ta@bao.ac.cn" not in rec.get("submitter","").lower():
                continue
            if not EP.match(title) or not TITLE.search(title):
                continue
            date,letter=ids[0]
            alerts[date+letter.lower()].append(rec)
    rows, audit=[],[]
    for name,group in sorted(alerts.items()):
        group.sort(key=lambda r:(r["createdOn"],r["circularId"]))
        first=group[0]
        reason=""
        if name in flags:
            reason="public correction/retraction requires adjudication"
        body=" ".join(first["body"].split())
        date=WHEN.search(body)
        point=POS.search(body)
        if not date or not point:
            reason=reason or "nonstandard timestamp or coordinate"
        if not reason:
            # Restrict the parse to the initial WXT detection, before any follow-up.
            if re.search(r"follow.up (?:observation|was performed)|(?:performed|conducted) (?:a|the) follow.up",body[:point.start()],re.I):
                reason="localization might describe a follow-up"
            elif re.search(r"neutrino|icecube",body[:max(date.end(),point.end())],re.I):
                reason="possible targeted follow-up"
            elif date[1].replace("-","").replace("/","")[2:]!=name[:6]:
                reason="event date differs from EP designation"
            elif not re.search(r"\bWXT\b",body[:point.end()],re.I):
                reason="instrument for first localization unclear"
        aid=str(first["circularId"])
        if not reason:
            ra=float(point[1].rstrip("."));dec=float(point[2].replace("−","-").rstrip("."))
            at=datetime.fromisoformat(date[1].replace("/","-")+"T"+date[2]).isoformat()+"Z"
            assert 0<=ra<360 and -90<=dec<=90
            rows.append(dict(source_id="EP"+name,ra_deg=ra,dec_deg=dec,utc=at,family="xray",anchor="WXT_detection",catalog="EP/GCN",gcn=aid))
        audit.append(dict(source_id="EP"+name,first_gcn=aid,matching_circulars="|".join(str(i["circularId"])for i in group),
                          flagged_gcn="|".join(str(x) for x in sorted(flags[name])),status=reason or "included"))
    dump("ep_events.csv",["source_id","ra_deg","dec_deg","utc","family","anchor","catalog","gcn"],rows)
    dump("ep_inventory.csv",["source_id","first_gcn","matching_circulars","flagged_gcn","status"],audit)

if __name__=="__main__":
    main()

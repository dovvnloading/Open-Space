"""Cross-check NASA Swift GRB classifications against autonomous GCN BAT position packets.
python bat.py <cache> fetch  downloads packets. python bat.py <cache>  projects.
"""
import csv
import hashlib
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from lxml import html
from optical import write

ROOT="https://gcn.gsfc.nasa.gov/"
def table(file):
    root=html.fromstring(file.read_bytes())
    return max(root.xpath("//table"),key=lambda t:len(t.xpath(".//tr"))).xpath(".//tr")

def text(node):
    return "".join(node.itertext()).strip()

def main():
    cache=Path(sys.argv[1]); notices={}
    for tr in table(cache/"gcnlatest"):
        cells=tr.xpath("td")
        if len(cells)<6 or not re.fullmatch(r"2[456]/\d\d/\d\d",text(cells[1])):continue
        notices[text(cells[0])]=(text(cells[1]),text(cells[2]),ROOT+cells[0].xpath(".//a/@href")[0])
    bursts=[]
    for tr in table(cache/"batfull"):
        c=tr.xpath("td")
        if len(c)<6 or not re.fullmatch(r"2[456]\d{4}[A-Z]?",text(c[0])):continue
        tid=text(c[2])
        if tid in notices:
            try:ra,dec=[float(next(el.itertext()).strip()) for el in c[3:5]]
            except ValueError:continue
            bursts.append((text(c[0]),tid,ra,dec,*notices[tid]))
    if "fetch" in sys.argv:
        def fetch(z):
            path=cache/("bat-"+z[1])
            if path.exists():return
            try:
                with urllib.request.urlopen(z[-1],timeout=40) as f: data=f.read()
                path.write_bytes(data)
            except OSError:pass
        with ThreadPoolExecutor(5) as pool:list(pool.map(fetch,bursts))
    rows=[];audit=[]
    for name,tid,ra,dec,date,clock,url in bursts:
        path=cache/("bat-"+tid);status="packet unavailable";digest=""
        if path.exists():
            raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
            packets=re.split(r"/{20,}",raw.decode(errors="replace"))
            pos=[p for p in packets if re.search(r"^NOTICE_TYPE:\s*Swift-BAT GRB Position\s*$",p,re.M)]
            if not pos:status="no flight GRB position notice"
            elif not any(re.search(r"This is (?:a|an) (?:rate|image) trigger",p) for p in pos):status="ambiguous flight trigger"
            elif any(re.search(r"this (?:is|was) (?:a |an )?(?:ground|slew|externally)",p,re.I) for p in pos):status="nonautonomous"
            else:status="included"
        if status=="included":
            rows.append(dict(source_id="GRB"+name,ra_deg=ra,dec_deg=dec,utc="20"+date.replace("/","-")+"T"+clock+"Z",family="xray",anchor="BAT_onboard_trigger",catalog="Swift BAT",trigger_id=tid))
        audit.append(dict(source_id="GRB"+name,trigger_id=tid,status=status,packet_url=url,sha256=digest))
    # Independent cross-check and gap fill from the mission's detailed BAT catalog.
    # This product stopped at 2025-06-05; never extrapolate its coverage.
    byid={r["trigger_id"]:r for r in audit}
    accepted={r["source_id"] for r in rows}
    url="https://swift.gsfc.nasa.gov/results/batgrbcat/summary_cflux/summary_general_info/summary_general.txt"
    for line in (cache/"batgeneral").read_text().splitlines():
        fields=[v.strip() for v in line.split("|")]
        if len(fields)<17 or not re.fullmatch(r"GRB2[456]\d{4}[A-Z]?", fields[0]):continue
        name,tid,_,ut,ra,dec=fields[:6]
        if name in accepted or fields[15] not in {"rate trigger","image trigger"}:continue
        try:
            ra,dec=float(ra),float(dec)
        except ValueError:continue
        rows.append(dict(source_id=name,ra_deg=ra,dec_deg=dec,utc=ut+"Z",family="xray",
                         anchor="BAT_onboard_trigger",catalog="Swift BAT",trigger_id=tid))
        accepted.add(name)
        if tid in byid:byid[tid]["status"]="included via mission summary"
        else:audit.append(dict(source_id=name,trigger_id=tid,status="included via mission summary",packet_url=url,sha256=""))
    if rows:write("bat_events.csv",sorted(rows,key=lambda r:r["utc"]))
    write("bat_inventory.csv",audit)

if __name__=="__main__":main()

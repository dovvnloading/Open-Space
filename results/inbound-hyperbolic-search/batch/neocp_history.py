"""Search a dated sample of unconfirmed NEOCP departures for archived astrometry."""
import base64
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date
import gzip
import hashlib
import json
from pathlib import Path
import re
from time import sleep
from urllib.parse import quote
from urllib.request import Request, urlopen

from bs4 import BeautifulSoup
from refit import fit

ROOT=Path(__file__).resolve().parent
TIMELINE="https://get-cp-orbits.jostjahn.de/object-timeline.html"
INITIAL=ROOT/"evidence/neocp_history_initial.json.gz"
PRIOR={r["id"]:r for r in json.loads(gzip.decompress(INITIAL.read_bytes()))["results"]} if INITIAL.exists() else {}
PSEUDO="https://www.projectpluto.com/neocp2/mpecs/{}.htm"
START,END="2026-10-19","2026-11-16"


def fetch(url):
    with urlopen(Request(url,headers={"User-Agent":"Open-Space astronomy research/1"}),timeout=28) as resp:
        raw=resp.read()
        return raw,{"url":resp.url,"http_date":resp.headers.get("Date"),"http_status":resp.status,
                    "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}


def get_target(item):
    ident,dt=item
    result={"id":ident,"departure":dt,"url":PSEUDO.format(quote(ident))}
    try:
        raw,receipt=fetch(result["url"])
        result["receipt"]=receipt
        soup=BeautifulSoup(raw,"html.parser")
        block=next((p for p in soup.find_all("pre") if "Astrometry:" in p.get_text()),None)
        if block is None:
            raise ValueError("No astrometry block")
        page=block.get_text().replace("\xa0"," ")
        lines=[]
        for ln in page.splitlines():
            time_field=re.search(r"(?:19|20)\d{2} \d{2} \d{2}\.\d+",ln)
            if not time_field or not ln[:time_field.start()].strip().startswith(ident):
                continue
            # Short temporary IDs are left aligned in the seven-byte provisional field.
            # Extract observation flags from the three bytes immediately before YYYY.
            flags=ln[time_field.start()-3:time_field.start()]
            record="     "+ident.ljust(7)+flags+ln[time_field.start():].rstrip()
            if len(record)==80 and len(record[77:].strip())==3:
                lines.append(record)
        # Some pages put the orbital elements in their own preformatted block.
        data=soup.get_text()
        el=re.search(r"(?ms)^Orbital elements:[^\n]*\n(.*?)(?=\n\s*(?:Residuals|Ephemeris|Ephemerides)|\Z)",data)
        body=el[1] if el else ""
        ecc=re.search(r"(?m)^\s*e\s+([0-9.]+)(?:\s+\+/-\s+([0-9.eE+-]+))?",body)
        result.update(published_e=float(ecc[1]) if ecc else None,
                      published_sigma_e=float(ecc[2]) if ecc and ecc[2] else None,
                      explicit_geocentric=bool(re.search(r"(?i)geocentric orbit|planetocentric|perigee",body)),
                      observations=lines,n_dates=len(set(l[15:25] for l in lines)),
                      stations=sorted(set(l[77:80] for l in lines)))
        if len(lines)>=3 and (result["n_dates"]>=3 or (result["published_e"] or 0)>=1):
            row={"spkid":ident,"full_name":ident,"tp":date.fromisoformat(dt[:10]).toordinal()+1721424.5}
            previous=PRIOR.get(ident,{})
            if previous.get("observations")==lines and previous.get("fit"):
                result["fit"]=previous["fit"]
            else:
                result["fit"]=fit((row,[(l+"\n").encode("ascii") for l in lines]))
    except Exception as exc:
        result["error"]=str(exc)
    sleep(.2)
    return result


def main():
    raw,receipt=fetch(TIMELINE)
    rows=BeautifulSoup(raw,"html.parser").find_all("tr")
    events={}
    for row in rows:
        line=" ".join(row.stripped_strings)
        m=re.search(r"(\d{4}-\d\d-\d\d) (\d\d:\d\d:\d\d) UTC\s+Departure\s+(\S+)\s+Not confirmed\b",line,re.I)
        if m and START<=m[1]<=END:
            events[m[3]]=m[1]+" "+m[2]+" UTC"
    if not events:
        raise ValueError("No unconfirmed timeline rows parsed: "+repr([" ".join(tr.stripped_strings) for tr in rows[:3]]))
    print("TARGETS",len(events),flush=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(get_target,sorted(events.items())))
    report={"selection":"Get-CP timeline, unconfirmed departures in inclusive UTC window",
            "start":START,"end":END,"timeline":receipt,"results":results}
    out=ROOT/"evidence"/"neocp_history.json.gz"
    out.write_bytes(gzip.compress(json.dumps(report,indent=2).encode(),mtime=0))
    print("HISTORY",len(results),sum(bool(r.get("observations")) for r in results),
          Counter(r.get("fit",{}).get("status") for r in results),flush=True)
    print("ARCHIVE",out.name,base64.b64encode(out.read_bytes()).decode())

if __name__=="__main__":
    main()

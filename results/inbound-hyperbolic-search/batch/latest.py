"""Refresh recent JPL arcs against the MPC Observations API; independently fit the returned positions."""
import base64
import csv
import gzip
import hashlib
import json
import re
from datetime import date, datetime, timezone
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from time import sleep

from refit import OUT, ROOT, aliases, cat, fit

ENDPOINT = "https://data.minorplanetcenter.net/api/get-obs"


def main():
    catalogue, _, _ = cat()
    targets = [r for r in catalogue if (r["last_obs"] or "")>"2026-09-11" or r["kind"]=="au"]
    with (ROOT.parent/"comparison/evidence/short_arc_screen.csv").open() as f:
        comparison = [r for r in csv.DictReader(f) if r["last_obs"] >= "2025-01-01"]
    targets += comparison
    records=[]
    for row in targets:
        des=row["pdes"]
        name=row["full_name"].strip()
        prefix=re.search(r"([CAPDX]/\d{4} [A-Z][A-Z]?\d*(?:-[A-Z])?)", name)
        permanent={"50788063":"1I","1003639":"2I","1004083":"3I"}.get(str(row["spkid"]))
        candidates=list(dict.fromkeys(([permanent] if permanent else [])+([prefix.group(1)] if prefix else [])+[des]+sorted(aliases(row))))
        record={"name":name,"spkid":row["spkid"],"url":ENDPOINT,"attempts":[],
                "selection":"recent near-bound comparison" if row in comparison else "catalogue refresh"}
        try:
            lines=[]
            for identifier in candidates:
                payload={"desigs":[identifier],"output_format":["OBS80"]}
                body=json.dumps(payload).encode()
                attempt={"request":payload}
                record["attempts"].append(attempt)
                request=Request(ENDPOINT,method="GET",data=body,headers={"Content-Type":"application/json","User-Agent":"Open-Space orbit research/1"})
                try:
                    with urlopen(request, timeout=35) as response:
                        content=response.read()
                        attempt.update(status_http=response.status,date=response.headers.get("Date"),
                                       response_sha256=hashlib.sha256(content).hexdigest())
                    lines=json.loads(content)[0]["OBS80"].splitlines()
                    if lines:
                        break
                except HTTPError as exc:
                    attempt.update(status_http=exc.code,message=exc.read(240).decode(errors="replace"))
                except (OSError,KeyError,ValueError,TypeError,IndexError) as exc:
                    attempt["error"]=str(exc)
            if not lines:
                raise ValueError("No observing records returned after identifier retries")
            accepted=[]
            packed=set(aliases(row))
            extended=re.fullmatch(r"(20\d{2}) ([A-HJ-Y])([A-HJ-Z])(\d+)",des)
            if extended and int(extended[4])>=620:
                letters="ABCDEFGHJKLMNOPQRSTUVWXYZ"
                alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
                seq=(int(extended[4])-620)*25+letters.index(extended[3])
                digits=""
                for _ in range(4):
                    seq,n=divmod(seq,62)
                    digits=alphabet[n]+digits
                packed.add("_"+alphabet[int(extended[1])-2000]+extended[2]+digits)
            for line in lines:
                if len(line)>=80 and (line[:5].strip() in packed or line[5:12].strip() in packed):
                    accepted.append((line[:80]+"\n").encode("ascii"))
            if not accepted:
                raise ValueError("No matching MPC designation bytes")
            if not row.get("tp"):
                # Date anchor for the comparison run only. Selected short-arc rows have no Tp column.
                row["tp"]=date.fromisoformat(row["last_obs"]).toordinal()+1721424.5
            days=sorted(l[15:25].decode() for l in accepted)
            record.update(first_observed=days[0],last_observed=days[-1],query_rows=len(lines),
                          matched_lines=len(accepted),expected_last_obs=row["last_obs"],fit=fit((row,accepted)))
        except Exception as exc:
            record["error"]=str(exc)
        print("MPC_REFIT",des,record.get("last_observed"), record.get("fit",{}).get("status"),record.get("error",""),flush=True)
        records.append(record)
        sleep(.15)
    p=OUT/"latest_mpc_refits.json"
    p.write_text(json.dumps(records,indent=2)+"\n")
    print("ARCHIVE",p.name+".gz",base64.b64encode(gzip.compress(p.read_bytes(),mtime=0)).decode())

if __name__=="__main__":
    main()

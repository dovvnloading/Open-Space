"""Retry current MPC arcs that timed out using a reproducible station/night sample."""
import base64
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import date
import gzip
import hashlib
import json
import os
from urllib.request import Request, urlopen

from refit import ROOT, aliases, cat, fit


def run(item):
    record,row=item
    key=str(row["spkid"])
    permanent={"50788063":"1I","1004083":"3I","1003639":"2I"}.get(key)
    successful=[a["request"]["desigs"][0] for a in record["attempts"] if a.get("status_http")==200]
    ident=permanent or successful[-1]
    req={"desigs":[ident],"output_format":["OBS80"]}
    base={"spkid":key,"name":record["name"],"selection":record["selection"],"url":record["url"],"request":req}
    try:
        with urlopen(Request(record["url"],method="GET",data=json.dumps(req).encode(),
                             headers={"Content-Type":"application/json","User-Agent":"Open-Space astronomy research"}),timeout=45) as response:
            raw=response.read()
            base.update(date=response.headers.get("Date"),status_http=response.status,
                        response_sha256=hashlib.sha256(raw).hexdigest())
        lines=json.loads(raw)[0]["OBS80"].splitlines()
        alias=aliases(row)
        groups={}
        accepted=0
        for line in lines:
            if len(line)<80 or not (line[:5].strip() in alias or line[5:12].strip() in alias):
                continue
            accepted+=1
            if line[14] not in {"C","B","A","P"," "} or not line[15:19].isdigit():
                continue
            groups.setdefault((line[15:25],line[77:80]),[]).append(line[:80])
        if len(groups)<3:
            raise ValueError("Under three usable station/nights")
        # Order by UTC and station and retain the median-timestamp record in each group.
        reps=[sorted(g,key=lambda l:l[15:32])[(len(g)-1)//2] for k,g in sorted(groups.items())]
        count=min(180,len(reps))
        inds=[round(j*(len(reps)-1)/(count-1)) for j in range(count)]
        selected=[reps[i] for i in inds]
        if not row.get("tp"):
            row["tp"]=date.fromisoformat(row["last_obs"]).toordinal()+1721424.5
        base.update(total_api_lines=len(lines),matched_lines=accepted,station_nights=len(groups),
                    selected_mpc80=selected,fit=fit((row,[(line+"\n").encode("ascii") for line in selected])))
    except Exception as exc:
        base["error"]=str(exc)
    print("STRATIFIED",base["name"],base.get("fit",{}).get("status"),base.get("error"),flush=True)
    return base


def main():
    os.environ["FO_UNPENALIZED"]="1"
    allrows,_,_=cat()
    with (ROOT.parent/"comparison/evidence/short_arc_screen.csv").open() as f:
        allrows+=list(csv.DictReader(f))
    by_id={str(row["spkid"]):row for row in allrows}
    full=json.loads(gzip.decompress((ROOT/"evidence/latest_mpc_refits.json.gz").read_bytes()))
    targets=[(r,by_id[str(r["spkid"])]) for r in full if r.get("fit",{}).get("status")=="fit timeout"]
    with ThreadPoolExecutor(max_workers=3) as pool:
        output=list(pool.map(run,targets))
    path=ROOT/"evidence"/"recent_stratified.json.gz"
    path.write_bytes(gzip.compress(json.dumps(output,indent=2).encode(),mtime=0))
    data=base64.b64encode(path.read_bytes()).decode()
    for i in range(0,len(data),4000):
        print("ARCHPART",path.name,i,data[i:i+4000])

if __name__=="__main__":
    main()

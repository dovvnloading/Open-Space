"""Inventory and inspect the newest complete UTC ingestion day of Rubin MPC submissions.

This feed mirrors MPC obs_sbn. The UTC ingestion date (object path) can differ
from the actual observation timestamp. It supplements the MPC ITF snapshot.
"""
import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib.parse import urlencode, quote
from urllib.request import urlopen

CUTOFF = "2026-11-17"
ROOT = "https://storage.googleapis.com"
BUCKET = "asteroid-institute-public"
PREFIX = "production/rubin/mpc/obs_sbn/daily/"


def main(out):
    query = {"prefix": PREFIX, "maxResults": 1000, "fields": "nextPageToken,items(name,size,generation,updated)"}
    pages = []
    dates = collections.defaultdict(list)
    while True:
        url = ROOT+"/storage/v1/b/"+BUCKET+"/o?"+urlencode(query)
        with urlopen(url, timeout=90) as response:
            raw = response.read()
        pages.append({"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()})
        obj = json.loads(raw)
        for item in obj.get("items", []):
            if not item["name"].endswith(".csv.gz"):
                continue
            hits = re.findall(r"20\d\d[-/]\d\d[-/]\d\d",item["name"])
            if hits:
                date = hits[-1].replace("/","-")
                if date < CUTOFF:
                    dates[date].append(item)
        if not obj.get("nextPageToken"):
            break
        query["pageToken"] = obj["nextPageToken"]
    result = dict(source="https://b612.ai/rubin-mpc-downloads/",cutoff_exclusive=CUTOFF,
                  observed_at=dt.datetime.now(dt.timezone.utc).isoformat(), listing_pages=pages,
                  n_dates=len(dates), first=min(dates,default=None), latest=max(dates,default=None),
                  ingestion_dates_in_primary_window=[x for x in sorted(dates) if "2026-10-19" <= x < CUTOFF],
                  retrieved=[])
    status, obs_dates = collections.Counter(), collections.Counter()
    inspected, superseded, feed_itf = {}, {}, []
    seen, duplicates = set(), 0
    raw_itf_in_window, exported_itf = 0, 0
    for date in sorted(dates,reverse=True):
        if date < "2026-10-19" and seen:
            break
        families = collections.defaultdict(list)
        for item in dates[date]:
            families[item["name"].rsplit("/",1)[0]].append(item)
        latest = max(families,key=lambda fam:max(i["updated"] for i in families[fam]))
        superseded[date] = sorted(k for k in families if k!=latest)
        inspected[date] = 0
        for item in sorted(families[latest],key=lambda i:i["name"]):
            url=ROOT+"/"+BUCKET+"/"+quote(item["name"])+ "?generation="+item["generation"]
            with tempfile.TemporaryFile() as f, urlopen(url,timeout=150) as response:
                digest=hashlib.sha256()
                n=0
                while block:=response.read(2**20):
                    n+=len(block)
                    digest.update(block)
                    f.write(block)
                if n!=int(item["size"]):
                    raise IOError("Truncated CSV gzip object")
                f.seek(0)
                with io.TextIOWrapper(gzip.GzipFile(fileobj=f)) as stream:
                    reader=csv.DictReader(stream)
                    result["fields"]=reader.fieldnames
                    for row in reader:
                        if row.get("stn")!="X05":
                            raise ValueError(f"Unexpected station {row.get('stn')}")
                        if row["id"] in seen:
                            duplicates+=1
                            continue
                        seen.add(row["id"])
                        inspected[date]+=1
                        obsdate=row["obstime"][:10]
                        obs_dates[obsdate]+=1
                        status[row["status"]]+=1
                        if row["status"]=="I" and "2026-10-19"<=obsdate<CUTOFF:
                            raw_itf_in_window+=1
                            line=row["obs80"].rstrip("\r\n")
                            if len(line)==80 and line[77:80]=="X05":
                                feed_itf.append(line)
                                exported_itf+=1
                result["retrieved"].append(dict(url=url,bytes=n,sha256=digest.hexdigest(),updated=item.get("updated")))
    feed_file=out/"x05-feed.mpc"
    feed_file.write_text("".join(line+"\n" for line in feed_itf))
    result.update(rows=len(seen),by_observation_date=dict(sorted(obs_dates.items())),by_status=dict(status),
                  ingestion_date_counts=inspected, older_csv_families={k:v for k,v in superseded.items() if v},
                  duplicate_ids=duplicates,primary_obstime_status_I=raw_itf_in_window, exported_obs80_rows=exported_itf,
                  exported_obs80_sha256=hashlib.sha256(feed_file.read_bytes()).hexdigest())
    (out/"rubin-feed.json").write_text(json.dumps(result,indent=2)+"\n")
    print("RUBIN_FEED_COUNTERS="+json.dumps({k:v for k,v in result.items() if k not in {"fields","retrieved","listing_pages"}}))

if __name__=="__main__":
    main(Path(sys.argv[1]))

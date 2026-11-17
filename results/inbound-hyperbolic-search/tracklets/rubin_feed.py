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
    status, rows = collections.Counter(), 0
    obs_dates = collections.Counter()
    if dates:
        for item in sorted(dates[max(dates)],key=lambda i:i["name"]):
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
                        rows+=1
                        obs_dates[row["obstime"][:10]]+=1
                        status[row.get("status","<absent>")]+=1
                result["retrieved"].append(dict(url=url,bytes=n,sha256=digest.hexdigest(),updated=item.get("updated")))
    result.update(rows=rows,by_observation_date=dict(sorted(obs_dates.items())),by_status=dict(status))
    (out/"rubin-feed.json").write_text(json.dumps(result,indent=2)+"\n")
    print("RUBIN_FEED="+json.dumps(result))

if __name__=="__main__":
    main(Path(sys.argv[1]))

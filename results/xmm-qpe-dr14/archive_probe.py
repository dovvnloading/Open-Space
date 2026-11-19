#!/usr/bin/env python3
"""Network/service probe on the GitHub Actions runner; NOT an archive screen."""
import datetime as dt
import hashlib
import io
import json
import time
import xml.sax.saxutils as sax
from pathlib import Path
import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)
ENDPOINT = "https://datalab.noirlab.edu/tap/sync"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "Open-Space-QPE-DR14-gated-host-selection/1.0"

def votable(rows):
    data = "".join(
        "<TR><TD>%s</TD><TD>%.8f</TD><TD>%.8f</TD></TR>" %
        (sax.escape(str(k)), ra, dec) for k, ra, dec in rows
    )
    return ("""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="1.3">
<RESOURCE type="results"><TABLE name="xray">
<FIELD name="detid" datatype="char" arraysize="*"/>
<FIELD name="ra" datatype="double" unit="deg"/>
<FIELD name="dec" datatype="double" unit="deg"/>
<DATA><TABLEDATA>""" + data +
            """</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>""").encode()

def request_probe(name, method, url, **kwargs):
    record = {"name": name, "url": url}
    start = time.monotonic()
    try:
        r = SESSION.request(method, url, timeout=(15, 60), **kwargs)
        body = r.content
        record.update(
            http_status=r.status_code, bytes=len(body),
            content_type=r.headers.get("content-type"),
            elapsed_s=round(time.monotonic()-start, 2),
            sha256=hashlib.sha256(body).hexdigest(),
            excerpt=body[:1600].decode("utf-8", "replace")
        )
        record["query_ok"] = bool(r.ok and "QUERY_STATUS\" value=\"ERROR" not in record["excerpt"]
            and not "<INFO" in record["excerpt"] and
            "error" not in record["excerpt"][:250].lower())
        if r.ok and "csv" in str(r.headers.get("content-type", "")).lower():
            record["csv_line_count"] = len(body.decode("utf-8","replace").splitlines())
    except Exception as exc:
        record.update(error=repr(exc), elapsed_s=round(time.monotonic()-start, 2), query_ok=False)
    return record

def main():
    rec={"utc": dt.datetime.now(dt.timezone.utc).isoformat(),
         "purpose":"Test anonymous LS DR10 bulk host crossmatch. No screen unless full audited cohort can be frozen.",
         "tests":[]}
    sample=[("GSN069", 19.786098, -34.191813),
            ("RXJ1301", 195.50054, 27.78267),
            ("eROQPE2", 38.704042, -44.325458)]
    query=("SELECT TOP 5 ls_id,ra,dec,type,brick_primary,flux_r "
           "FROM ls_dr10.tractor WHERE 't'=Q3C_RADIAL_QUERY("
           "ra,dec,19.786098,-34.191813,0.002)")
    rec["tests"].append(request_probe("tap_cone", "POST", ENDPOINT,
       data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":query}))
    uploaded=votable(sample)
    adql=("SELECT TOP 30 x.detid, x.ra AS xra, x.dec AS xdec,"
          " t.ls_id, t.ra AS gra, t.dec AS gdec, t.type, t.brick_primary, t.flux_r"
          " FROM TAP_UPLOAD.xray AS x, ls_dr10.tractor AS t"
          " WHERE 't'=Q3C_JOIN(x.ra,x.dec,t.ra,t.dec,0.000555555556)"
          " AND t.type IN ('EXP','DEV','SER','REX') AND t.brick_primary=1")
    rec["tests"].append(request_probe("tap_upload", "POST", ENDPOINT,
       data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv",
             "QUERY":adql,"UPLOAD":"xray,param:xray"},
       files={"xray":("xray.xml",uploaded,"application/x-votable+xml")}))
    sql=("WITH x(detid,xra,xdec) AS (VALUES "
         "('GSN069',19.786098,-34.191813),('RXJ1301',195.50054,27.78267)) "
         "SELECT x.detid,t.ls_id,t.ra,t.dec,t.type,t.brick_primary "
         "FROM x,ls_dr10.tractor t WHERE "
         "q3c_join(x.xra,x.xdec,t.ra,t.dec,0.000555555556) "
         "AND t.type IN ('EXP','DEV','SER','REX') LIMIT 20")
    rec["tests"].append(request_probe("query_sql_values", "GET",
        "https://datalab.noirlab.edu/query",
        params={"sql":sql,"ofmt":"csv"}))
    rec["all_sites_responded"] = all("http_status" in r for r in rec["tests"])
    (OUT/"archive_probe_01.json").write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps(rec,indent=2),flush=True)

if __name__=="__main__":
    main()

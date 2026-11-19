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
         "purpose":"Try independent SDSS/Legacy host services after LS DR10 TAP outages; no sky screen.",
         "tests":[]}
    sdss="https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
    single=("SELECT TOP 20 p.objID,p.ra,p.dec,p.type,p.mode,p.clean "
            "FROM PhotoPrimary p WHERE p.ra BETWEEN 195.499 AND 195.503 "
            "AND p.dec BETWEEN 27.780 AND 27.785 AND p.type=3")
    rec["tests"].append(request_probe("sdss_single", "GET", sdss,
       params={"cmd":single,"format":"csv"}))
    sql=("SELECT x.detid, p.objID, p.ra, p.dec, p.type, p.mode, p.clean "
         "FROM (VALUES ('RXJ1301',195.50054,27.78267),"
         "('GSN069',19.786098,-34.191813)) AS x(detid,ra,dec) "
         "CROSS APPLY dbo.fGetNearbyObjEq(x.ra,x.dec,0.0333333) AS n "
         "JOIN PhotoPrimary AS p ON n.objID=p.objID "
         "WHERE p.type=3")
    rec["tests"].append(request_probe("sdss_batch", "GET", sdss,
       params={"cmd":sql,"format":"csv"}))
    direct=("SELECT TOP 15 p.objID,p.ra,p.dec,p.type FROM PhotoPrimary p "
            "JOIN dbo.fGetNearbyObjEq(195.50054,27.78267,0.0333333) AS n "
            "ON p.objID=n.objID WHERE p.type=3")
    rec["tests"].append(request_probe("sdss_nearby", "GET", sdss,
       params={"cmd":direct,"format":"csv"}))
    beta="https://beta.datalab.noirlab.edu/tap/sync"
    rec["tests"].append(request_probe("noirlab_beta_cone", "GET", beta,
       params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv",
               "QUERY":"SELECT TOP 5 ls_id,ra,dec,type FROM ls_dr10.tractor WHERE 't'=Q3C_RADIAL_QUERY(ra,dec,19.786098,-34.191813,0.002)"}))
    (OUT/"archive_probe_03.json").write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps(rec,indent=2),flush=True)

if __name__=="__main__":
    main()

"""Legacy spectrum identifiers and SDSS model redshifts for flagged MWS records."""
import hashlib
import json
import urllib.parse
import urllib.request
from screen import HERE

sql = ("SELECT plate,mjd,fiberID,run2d,ra,dec,z,zErr,subClass FROM SpecObjAll WHERE "
       "(ra BETWEEN 210.508 AND 210.513 AND dec BETWEEN 2.900 AND 2.904) OR "
       "(ra BETWEEN 234.424 AND 234.428 AND dec BETWEEN 43.787 AND 43.791)")
url="https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?"+urllib.parse.urlencode(dict(cmd=sql,format="csv"))
with urllib.request.urlopen(url,timeout=90) as response:
    raw=response.read()
    if not raw.startswith(b"#Table1"):raise ValueError(raw[:200])
    report=dict(url=url,status=response.status,server_date=response.getheader("Date"),
                raw_sha256=hashlib.sha256(raw).hexdigest())
(HERE/"legacy_rv_hits.csv").write_bytes(raw.split(b"\n",1)[1])
(HERE/"legacy_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
